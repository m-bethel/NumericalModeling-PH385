/*
main.cpp - Main File
This program models the second order phase transition of a 20x20x20 lattice
using Ising Model and Monte Carlo Methods

Simulation Parameters:
-20x20x20 lattice.
-Periodic boundary conditions.
-200 Monte Carlo sweeps.
-Temperatures ranging from 0.01 to 4.00, increasing by 0.01 (units are in J/k_b ).
-Plot the second order phase transition.
-Use the simulated data near (but not above) the critical temperature to find the coefficient \Beta

Author: Miles Bethel    (miles.d.bethel@gmail.com)
Date: 04/01/2026

================================================================================
Discussion on Computational Trade-offs & Limitations:
================================================================================
- Limitations of the Method: 
  The Monte Carlo (Metropolis) method can suffer from "critical slowing down" near the 
  phase transition, where the system takes a very long time to reach equilibrium. It is 
  also subject to finite-size effects; a 20^3 lattice is an approximation of a macroscopic 
  system. This finite size causes the theoretical sharp phase transition (at Tc) to smooth 
  out and slightly shift.

- Ideal Parameters: 
  For a high-quality, scientifically accurate curve, a lattice size of 20x20x20 or larger, 
  coupled with 200+ sweeps per temperature step (and potentially a separate "burn-in" 
  phase to reach thermal equilibrium before measuring), and a fine temperature step of 
  0.01 are recommended.

- Lower Limits & Saving Computation Time: 
  If you need to iterate quickly (e.g., debugging), you can scale down parameters. Because 
  this simulation uses a "warm start" (retaining the lattice state from the previous 
  temperature), small temperature steps (e.g., 0.01) mean the system doesn't have to adjust 
  much between steps. In these cases, you can "get away" with as few as 50 sweeps. However, 
  if you increase your step size (e.g., 0.05) or decrease your grid size (e.g., 10x10x10), 
  fluctuations become more extreme, and you may need to increase your sweeps back to 
  100-200 to average out the statistical noise properly.
================================================================================
*/

#include "Lattice.h" // Include the custom Lattice class definition
#include <vector>    // Required for using the std::vector container
#include <numeric>   // Required for math operations like std::accumulate
#include <random>    // Required for the random number generators
#include <iostream>  // Required for printing output to the console
#include <fstream>   // Required for writing data to an external file

using namespace std; // Allows using standard library components without the std:: prefix

int main() {

    int sweeps = 200;       // The number of times to iterate over the entire lattice per temperature
    int size = 20;          // The length of one side of the 3D cube (20x20x20 = 8000 spins)
    double H = 0;           // The external magnetic field strength (0 for standard Ising model)
    double start = 0.01;    // The starting temperature (near absolute zero)
    double stop = 4.00;     // The ending temperature (past the theoretical 3D phase transition)
    double step = 0.01;     // The increment to increase the temperature by each loop iteration

    double maxChi = 0.00;   // A variable to keep track of the highest observed magnetic susceptibility
    double Tc;              // A variable to store the estimated Critical Temperature

    Lattice lattice(size);  // Instantiate the lattice object with the specified dimensions

    // Choose one of the following, one starts with the lattice with 
    // randomized spin, the other will all spin +1.
    // lattice.initializeRandom(); // Starts in a high-energy disordered state
    
    // We use initializeAligned() for a "warm start" at low temperatures. 
    // This perfectly aligns the spins (ground state) preventing frozen multi-domain states at T=0.01.
    lattice.initializeAligned(); 

    ofstream outFile("ising_data.csv"); // Create and open a CSV file to save the data
    // Write the column headers to the first row of the CSV
    outFile << "Temperature,Magnetization,Chi" << endl; 

    // Loop through temperatures from start to stop, incrementing by step
    for (double T = start; T <= stop; T+= step) {
        double magSum = 0;      // Accumulator for the magnetization values across all sweeps
        double magSqSum = 0.0;  // Accumulator for the square of magnetization (needed for variance/Chi)
        
        // Perform the specified number of Monte Carlo sweeps for this temperature
        for ( int s = 0; s < sweeps; ++s ){
            lattice.performSweep(T, H); // Attempt to flip every spin in the lattice once
            double M = lattice.calculateMagnetization(); // Get the total magnetization of the lattice
            magSum += M;        // Add the current magnetization to our sum
            magSqSum += (M*M);  // Add the squared magnetization to our squared sum
        }
        
        // Calculate the average magnetization over all sweeps at this temperature
        double currentMag = (magSum / sweeps); 
        // Calculate the average of the squared magnetization
        double currentMag2 = (magSqSum / sweeps); 
        
        // Calculate Magnetic Susceptibility (Chi) using the fluctuation-dissipation theorem formula.
        // Variance = <M^2> - <M>^2. Multiplied by 100 to scale the output for easier plotting.
        double chi = 100*(currentMag2 - (currentMag*currentMag)) / T; 
        
        // Check if we have found a new peak in susceptibility. We only check above T=2.0 
        // to avoid any potential low-temperature mathematical artifacts/noise.
        if (T > 2.0 && chi > maxChi){
            maxChi = chi; // Update the new highest Chi value
            Tc = T;       // Record the temperature where this peak occurred (Critical Temp)
        }

        // Write the calculated averages for this temperature to the CSV file
        outFile << T << "," << currentMag << "," << chi << endl;
        
        // Print a progress update to the console roughly every 50 steps
        // (T * 100) casts the double to an int (e.g., 2.50 becomes 250), so % 50 works.
        if (int (T * 100) % 50 == 0) cout << "Running Simulation... T = " << T << endl;
    }
    
    outFile.close(); // Close the CSV file to ensure data is saved properly
    
    // Output the final estimated Critical Temperature based on the susceptibility peak
    cout << "Critical Temp: " << Tc << endl; 
    
    return 0; // End the program successfully
}