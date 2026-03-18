/*
Particle.h

Declares the Particle class, which represents a single particle in the
molecular dynamics simulation. Each particle tracks its own 3D position.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/2/2026
*/

#pragma once

class Particle
{
    private:
        // Current x, y, z position of the particle in 3D space
        // LIMITATION: fixed at 3 dimensions — class cannot be used for 2D or nD simulations
        // without rewriting the array size and all associated loops
        double m_position[2];
        double m_velocity[2];
        double m_force[2];
        double m_forcePrev[2];
        int m_gridSize;

    public:
        // Constructor — initializes particle at the origin (0, 0, 0)
        // LIMITATION: starting position is hardcoded to the origin
        // A more flexible design would accept a starting position as a parameter
        Particle(int gridSize, double x, double y);

        const double* getVelocity();
        const double* getForce();
        void capVelocity(double maxSpeed);
        void addForce(double fx, double fy);
        void setForce(double fx, double fy);
        void scaleVelocity(double factor);
        void savePrevForce();
        void updatePosition(double dt);
        void updateVelocity(double dt);
    

        double kineticEnergy();

        // Tracks the particle as it exits the container, but reinserts it at the "other side"
        // of the grid with the same condition. Periodic Boundary conditions
        void periodic();

        // Returns a read-only pointer to the internal position array
        // LIMITATION: returns a raw pointer — caller must know the array is size 3
        // A safer alternative would be std::array<double,3> which carries its size
        const double* getPosition();
};