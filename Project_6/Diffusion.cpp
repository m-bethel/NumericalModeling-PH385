/*
Diffusion.cpp

Implements the Diffusion class. Manages the simulation loop across all
particles, records mean squared displacement at each step, exports results
to a CSV file, and generates a gnuplot script to visualize the output.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 02/18/2026
*/

#include "Diffusion.h"
#include <iostream>
#include <iomanip>
#include <fstream>  // ofstream for writing CSV and gnuplot script
#include <cmath>
#include <vector>

using namespace std;

// Constructor — initializes simulation parameters and creates all particles
// Each Particle is default constructed at the origin (0, 0, 0)
Diffusion::Diffusion(int numParticles, int steps)
    : m_numParticles(numParticles),     // store particle count
      m_steps(steps),                   // store step count
      m_particles(numParticles),        // create numParticles Particle objects
      m_fileName("msd.csv")             // default output file name
{}

// Runs the full diffusion simulation
// For each step: moves every particle, reflects off walls,
// then calculates and stores the mean squared displacement
// LIMITATION: particles are processed sequentially, one at a time
// A parallel implementation (e.g. using OpenMP) could significantly
// speed this up for large particle counts
void Diffusion::simulate(double step_size)
{
    // Outer loop — iterate over each time step
    for (int j = 0; j < m_steps; j++)
    {
        double sum = 0.0;   // accumulates r² across all particles for this step

        // Inner loop — update every particle for this step
        for (Particle& p : m_particles)
        {
            p.step(step_size);   // move particle in a random direction
            p.reflect();         // bounce off container walls if needed
            sum += p.rSquared(); // add this particle's r² to the running sum
        }

        // Calculate mean squared displacement for this step and store it
        // LIMITATION: only the average is stored — individual particle
        // trajectories are not recorded, so they cannot be analyzed later
        // or used to produce a 3D position animation
        m_msd.push_back(sum / m_numParticles);
    }
}

// Writes simulation results to a CSV file and generates a gnuplot plot
void Diffusion::output()
{
    // Open the output CSV file
    // LIMITATION: no error checking on whether the file opened successfully
    // If the file cannot be created (e.g. bad permissions), the program will
    // silently write nothing and produce an empty or missing plot
    ofstream outFile(m_fileName);

    // Write the header row
    outFile << "step,msd" << endl;

    // Write one row per step: step number and mean squared displacement
    for (int i = 0; i < m_steps; i++)
    {
        outFile << i << "," << m_msd[i] << endl;
    }

    outFile.close();

    cout << "\nGenerating 2D plot with gnuplot...\n";

    // Create a gnuplot script file to plot <r²> vs step
    ofstream gnuplotScript("plot_script.gnu");

    gnuplotScript << "set title 'Mean Squared Displacement vs Step'\n";
    gnuplotScript << "set xlabel 'Step'\n";                          // x axis label
    gnuplotScript << "set ylabel '<r^2>'\n";                         // y axis label
    gnuplotScript << "set datafile separator \",\"\n";               // tell gnuplot file is comma separated
    gnuplotScript << "set key noautotitle\n";                        // suppress automatic legend entry
    gnuplotScript << "plot '" << m_fileName << "' skip 1 using 1:2 with lines linewidth 2\n"; // plot columns 1 and 2, skip header
    gnuplotScript << "pause -1 'Press any key to close'\n";          // keep window open until user closes it

    gnuplotScript.close();

    // Execute the gnuplot script to display the plot
    // LIMITATION: system() is a blocking call on some platforms and non-blocking
    // on others — behavior may differ between Linux, macOS, and Windows
    // LIMITATION: system() has security implications if the filename were ever
    // derived from user input — safe here since it is hardcoded
    system("gnuplot plot_script.gnu");
}