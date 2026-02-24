/*
Diffusion.h

Declares the Diffusion class, which manages the full diffusion simulation.
Owns a collection of Particle objects, runs the simulation loop, records
mean squared displacement at each step, and handles all file output.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 02/18/2026
*/

#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <string>       // std::string for file name storage
#include <vector>       // std::vector for particle collection and MSD storage
#include "Particle.h"   // Particle class used in the simulation

class Diffusion
{
    private:
        // Total number of steps to simulate
        int m_steps;

        // Total number of particles in the simulation
        int m_numParticles;

        // Collection of all particles, each starting at the origin
        // LIMITATION: all particles are stored in memory simultaneously
        // For very large particle counts this could become a memory concern
        std::vector<Particle> m_particles;

        // Stores the mean squared displacement <r²> recorded at each step
        // LIMITATION: all MSD values are stored in memory for the full simulation
        // For very long simulations this grows linearly with step count
        std::vector<double> m_msd;

        // Name of the output CSV file for MSD data
        // LIMITATION: hardcoded to "msd.csv" in the constructor — if multiple
        // simulations are run in the same directory, the file will be overwritten
        std::string m_fileName;

    public:
        // Constructor — creates numParticles particles and sets up the simulation
        Diffusion(int numParticles, int steps);

        // Runs the full simulation: steps and reflects all particles each
        // iteration and records the mean squared displacement at each step
        // LIMITATION: simulate() can only be called once — calling it again
        // without resetting would append to m_msd and give incorrect averages
        void simulate(double step_size = 1.0);

        // Writes MSD data to a CSV file and runs a gnuplot script to plot it
        // LIMITATION: output() depends on gnuplot being installed and accessible
        // in the system PATH — will silently fail if gnuplot is not available
        void output();
};

#endif // DIFFUSION_H