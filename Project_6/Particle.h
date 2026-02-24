/*
Particle.h

Declares the Particle class, which represents a single particle in the
diffusion simulation. Each particle tracks its own 3D position and contains
its own random number generator for independent stepping.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 02/18/2026
*/

#pragma once
#include <random>   // mt19937 random number engine and uniform distribution

class Particle
{
    private:
        // Current x, y, z position of the particle in 3D space
        // LIMITATION: fixed at 3 dimensions — class cannot be used for 2D or nD simulations
        // without rewriting the array size and all associated loops
        double m_position[3];

        // Mersenne Twister random number engine, seeded uniquely per particle
        // LIMITATION: each particle owns its own RNG, which increases memory usage
        // significantly when simulating large numbers of particles (e.g. 2500+)
        // A shared RNG passed in from Diffusion would be more memory efficient
        std::mt19937 m_rng;

        // Uniform distribution over [0, 1] used to generate random angles
        std::uniform_real_distribution<double> m_dist;

    public:
        // Constructor — initializes particle at the origin (0, 0, 0)
        // LIMITATION: starting position is hardcoded to the origin
        // A more flexible design would accept a starting position as a parameter
        Particle();

        // Takes one random step of a given size in a uniformly random
        // direction over steradians (spherical coordinates)
        void step(double step_size = 1.0);

        // Reflects the particle off the container walls at ±0.5
        // LIMITATION: container boundaries are hardcoded to ±0.5
        // A more flexible design would accept boundary size as a parameter
        void reflect();

        // Returns the squared distance from the origin (x² + y² + z²)
        // Used to calculate mean squared displacement without needing sqrt
        // LIMITATION: always measures displacement from the origin (0,0,0)
        // If the starting position were ever changed, this would give incorrect results
        double rSquared();

        // Returns a read-only pointer to the internal position array
        // LIMITATION: returns a raw pointer — caller must know the array is size 3
        // A safer alternative would be std::array<double,3> which carries its size
        const double* getPosition();
};