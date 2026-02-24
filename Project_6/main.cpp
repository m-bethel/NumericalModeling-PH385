/*
Main File

Simulates the random walk diffusion of particles in a closed cubic container.
Prompts the user for simulation parameters, runs the diffusion model, and
outputs the mean squared displacement data with an automated gnuplot graph.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 02/18/2026
*/

#include <iostream>     // cout, cin
#include <limits>       // numeric_limits for input buffer clearing
#include "Diffusion.h"  // Diffusion class — manages the full simulation
#include "Particle.h"   // Particle class — represents a single diffusing particle

using namespace std;

int main()
{
    int steps;      // number of simulation steps
    int particles;  // number of particles to simulate
    char choice;    // user's menu selection

    // Prompt user to accept default conditions or enter their own
    cout << "Would you like to accept default initial conditions? [y,n]\n";
    cout << "Particles: 2500 \nSteps: 1000";

    cin >> choice;

    // Clear the input buffer to avoid leftover characters affecting later input
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    // Use default conditions or prompt user for custom values
    if (choice == 'y' || choice == 'Y')
    {
        cout << "Using default conditions" << endl;
        steps = 1000;       // default step count per project specification
        particles = 2500;   // default particle count per project specification
    }
    else
    {
        // LIMITATION: no input validation — if the user enters a negative number,
        // zero, or a non-integer, the simulation will behave incorrectly or crash
        cout << "Enter the number of particles you want to model: ";
        cin >> particles;
        cout << "Enter the number of steps to model: ";
        cin >> steps;
    }

    // Create the diffusion simulation with the specified number of
    // particles and steps — all particles begin at the origin (0, 0, 0)
    Diffusion simulation(particles, steps);

    // Run the simulation with a step size of 0.025 units as per project spec
    // LIMITATION: step size is hardcoded here — ideally this would also be
    // a user-configurable parameter, especially for the investigation portion
    // of the project where different step sizes are compared
    simulation.simulate(0.025);

    // Write MSD data to CSV and generate the gnuplot graph
    simulation.output();

    return 0;
}