/*
Lattice.h - Header File
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
*/
#ifndef LATTICE_H  // Include guard: Prevents this file from being included multiple times
#define LATTICE_H

#include <vector>  // Required to use std::vector for our internal array

// The Lattice class defines the structure and behaviors of our 3D Ising model
class Lattice {
private:
    int L;        // The length of one dimension of the cube (e.g., 20)
    int N;        // The total number of spins in the lattice (L * L * L)
    std::vector<int> spins; // A flat 1D vector used to store the 3D lattice for memory speed
    
    // A private helper function to convert 3D coordinates (x,y,z) into a 1D index
    int get_idx(int x, int y, int z) const;

public:
    // Constructor: Initializes the lattice with a given dimension size
    Lattice(int size);

    // Fills the lattice with a random mix of +1 and -1 spins (High T state)
    void initializeRandom();
    
    // Fills the lattice with all +1 spins (Absolute Zero state)
    void initializeAligned();

    // Runs one complete Metropolis Monte Carlo sweep across the lattice
    void performSweep(double T, double H);

    // Calculates the sum of the 6 nearest neighbors for a specific spin using PBC
    int getNeighborSum(int x, int y, int z) const;

    // Multiplies the spin at a specific location by -1 (flips it)
    void flipSpin(int x, int y, int z);

    // Returns the value (+1 or -1) of a spin at a specific location
    int getSpin(int x, int y, int z) const;

    // Computes the absolute average magnetization of the entire lattice
    double calculateMagnetization() const;

    // A getter function to return the dimension size (L)
    int getSize() const { return L; } 
};

#endif // Ends the include guard