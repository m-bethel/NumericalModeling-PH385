/*
Particle.cpp

Implements the Particle class. Handles the core physics of the simulation:
Paritcle movement and periodic motion of particles in grid

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 02/18/2026
*/

#include "Particle.h"
#include <iostream>
#include <fstream>
#include <cmath>    // sin, cos, acos, sqrt, M_PI
#include <vector>
#include <cstdlib>

using namespace std;

// Constructor — initializes position to the origin and seeds the RNG
Particle::Particle(int gridSize, double x, double y)
    : m_position{x, y},   // all particles start at (10, 10)
      m_gridSize(gridSize),
      m_velocity{((double)rand()/RAND_MAX)*0.4 - 0.2,
                 ((double)rand()/RAND_MAX)*0.4 - 0.2},
      m_force{0.0,0.0},
      m_forcePrev{0.0,0.0}
{}

double Particle::kineticEnergy()
{
    double energy = 0.0;
    energy = (m_velocity[0]* m_velocity[0]+m_velocity[1]* m_velocity[1])/2;
    return energy;
}
const double* Particle::getVelocity()
{
    return m_velocity;
}
void Particle::capVelocity(double maxSpeed)
{
    for (int i = 0; i < 2; i++)
    {
        if (m_velocity[i] > maxSpeed)  m_velocity[i] = maxSpeed;
        if (m_velocity[i] < -maxSpeed) m_velocity[i] = -maxSpeed;
    }
}
void Particle::savePrevForce()
{
    copy(begin(m_force), end(m_force), begin(m_forcePrev));
}
void Particle::updatePosition(double dt)
{
    for (int i=0; i<2; i++)
    {
        m_position[i] = m_position[i] + m_velocity[i]*dt+0.5*m_force[i]*dt*dt;
    }

}
void Particle::updateVelocity(double dt)
{
    for (int i=0;i<2;i++)
    {
        m_velocity[i] = m_velocity[i] + 0.5*(m_force[i] + m_forcePrev[i])*dt;
    }
}
const double* Particle::getForce()
{
    return m_force;
}
void Particle::addForce(double fx, double fy)
{
    m_force[0] += fx;
    m_force[1] += fy;
}
void Particle::setForce(double fx, double fy)
{
    m_force[0] = fx;
    m_force[1] = fy;
}
void Particle::scaleVelocity(double factor)
{
    for (int i = 0; i < 2; i++)
    {
        m_velocity[i] *= factor;
    }
}
const double* Particle::getPosition()
{
    return m_position;
}
void Particle::periodic()
{
    for (int i = 0; i < 2; i++)
    {
        // If the particle has moved past the right/top boundary
        if (m_position[i] >= m_gridSize) 
        {
            m_position[i] -= m_gridSize;
        }
        // If the particle has moved past the left/bottom boundary
        else if (m_position[i] < 0) 
        {
            m_position[i] += m_gridSize;
        }
    }
}