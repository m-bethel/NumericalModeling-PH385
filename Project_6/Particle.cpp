/*
Particle.cpp

Implements the Particle class. Handles the core physics of the simulation:
random directional stepping over steradians and reflection off container walls.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 02/18/2026
*/

#include "Particle.h"
#include <iostream>
#include <fstream>
#include <cmath>    // sin, cos, acos, sqrt, M_PI
#include <vector>

using namespace std;

// Constructor — initializes position to the origin and seeds the RNG
Particle::Particle()
    : m_position{0.0, 0.0, 0.0},       // all particles start at (0, 0, 0)
      m_rng(std::random_device{}()),    // seed with a non-deterministic random value
                                        // LIMITATION: std::random_device may fall back to a
                                        // pseudo-random seed on some systems (e.g. MinGW on Windows)
                                        // which could produce identical sequences across particles
      m_dist(0, 1.0)                    // uniform distribution over [0, 1]
{}

// Takes one step in a uniformly random direction over steradians
// Uses spherical coordinates to ensure equal probability in all directions
void Particle::step(double step_size)
{
    // Generate azimuthal angle phi uniformly over [0, 2π]
    double phi = 2 * M_PI * m_dist(m_rng);

    // Generate polar angle theta using inverse CDF to ensure
    // uniform distribution over the surface of a sphere
    // LIMITATION: acos() is defined over [-1, 1] — if floating point rounding
    // causes (1 - 2*random) to slightly exceed this range, acos() returns NaN
    // This is rare but theoretically possible
    double theta = acos(1 - 2 * m_dist(m_rng));

    // Convert spherical to cartesian and scale by step size
    double dx = sin(theta) * cos(phi) * step_size;
    double dy = sin(theta) * sin(phi) * step_size;
    double dz = cos(theta) * step_size;

    // Apply displacement to current position
    m_position[0] += dx;
    m_position[1] += dy;
    m_position[2] += dz;
}

// Returns the squared distance from the origin: x² + y² + z²
// No sqrt needed since we want r² directly for the MSD calculation
double Particle::rSquared()
{
    double r_squared = 0.0;

    // Sum the square of each coordinate
    for (int i = 0; i < 3; i++)
    {
        r_squared += m_position[i] * m_position[i];
    }
    return r_squared;
}

// Returns a read-only pointer to the internal position array
const double* Particle::getPosition()
{
    return m_position;
}

// Reflects the particle off the walls of the cubic container at ±0.5
// If a particle overshoots a wall, it bounces back by the overshoot amount
// LIMITATION: only handles a single reflection per step — if step_size is large
// enough for a particle to overshoot a wall and hit the opposite wall in the
// same step, the reflection will be incorrect. With step_size = 0.025 and
// a container of width 1.0 this is not an issue, but larger step sizes could
// cause particles to escape the container or behave unphysically
void Particle::reflect()
{
    for (int i = 0; i < 3; i++)
    {
        // Reflect off the positive wall (+0.5)
        if (m_position[i] > 0.5)
        {
            m_position[i] = 1.0 - m_position[i];
        }

        // Reflect off the negative wall (-0.5)
        if (m_position[i] < -0.5)
        {
            m_position[i] = -1.0 - m_position[i];
        }
    }
}