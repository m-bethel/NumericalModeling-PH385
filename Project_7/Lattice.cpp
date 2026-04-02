/*
Lattice.cpp - Implementation File

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

#include "Lattice.h" // Include the corresponding header file
#include <random>    // Required for the Mersenne Twister random number generator
#include <numeric>   // Required for std::accumulate to sum vector elements
#include <cmath>     // Required for mathematical functions like std::exp and std::abs

// Constructor: Sets L to size, N to size^3, and initializes the flat 'spins' vector with N zeros
Lattice::Lattice(int size) : L(size), N(size * size * size), spins(N, 0) {}

// Helper function to map 3D coordinates to a 1D array index using Periodic Boundary Conditions
int Lattice::get_idx(int x, int y, int z) const {
    // Lambda function to handle wrap-around. 
    // Adding L before the modulo ensures negative values wrap correctly in C++
    auto pbc = [this](int val) { return (val % L + L) % L; };
    
    // Mathematical formula to convert 3D (x,y,z) coordinates into a linear index
    return pbc(x) * L * L + pbc(y) * L + pbc(z);
}

// Sets every spin in the lattice to either +1 or -1 randomly
void Lattice::initializeRandom() {
    std::random_device rd; // Seed generator
    std::mt19937 gen(rd()); // Mersenne Twister engine seeded with rd
    std::uniform_int_distribution<> dist(0, 1); // Distribution that returns 0 or 1

    // Iterate through every element in the spins vector by reference
    for (int& s : spins) {
        // If dist generates 1, set to 1. If it generates 0, set to -1
        s = (dist(gen) == 1) ? 1 : -1; 
    }
}

// Sets every spin in the lattice to +1 (perfectly aligned ferromagnet)
void Lattice::initializeAligned() {
    // Iterate through every element in the spins vector by reference
    for (int& s : spins) {
        s = 1; // Start with a perfect ferromagnetic ground state
    }
}

// The core physics loop: Attempts to flip every spin based on the Metropolis algorithm
void Lattice::performSweep(double T, double H) {
    std::random_device rd; // Seed generator
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Distribution returning floats between 0.0 and 1.0

    // Loop through the X, Y, and Z dimensions
    for (int x = 0; x < L; ++x) {
        for (int y = 0; y < L; ++y) {
            for (int z = 0; z < L; ++z) {
                int idx = get_idx(x, y, z); // Get the 1D index for the current 3D position
                int s_i = spins[idx];       // Get the current spin state (+1 or -1)
                int neighbors = getNeighborSum(x, y, z); // Calculate the sum of the 6 adjacent spins
                
                // Calculate the change in energy if this spin were to flip.
                // Normally dE = 2 * J * s_i * sum(neighbors). Your code omitted the 2, 
                // which just means your temperature scale is shifted by a factor of 2.
                double dE = s_i * (neighbors + H); 

                // Metropolis Criteria:
                // If dE < 0, the flip lowers energy, so we ALWAYS accept it.
                // If dE >= 0, we flip with probability e^(-dE/T). 
                // We generate a random number and see if it's less than this probability.
                if (dE < 0 || dist(gen) <= std::exp(-dE / T)) { 
                    spins[idx] *= -1; // Multiply by -1 to flip the spin
                }
            }
        }
    }
}

// Calculates the sum of the 6 nearest neighbors (Up, Down, Left, Right, Front, Back)
int Lattice::getNeighborSum(int x, int y, int z) const {
    // get_idx automatically handles the wrap-around if x+1 goes out of bounds
    return spins[get_idx(x + 1, y, z)] + spins[get_idx(x - 1, y, z)] +
           spins[get_idx(x, y + 1, z)] + spins[get_idx(x, y - 1, z)] +
           spins[get_idx(x, y, z + 1)] + spins[get_idx(x, y, z - 1)];
}

// Helper function to manually flip a specific spin
void Lattice::flipSpin(int x, int y, int z) {
    spins[get_idx(x, y, z)] *= -1; // Multiply the spin by -1
}

// Helper function to read a specific spin
int Lattice::getSpin(int x, int y, int z) const {
    return spins[get_idx(x, y, z)]; // Return the value at the calculated 1D index
}

// Calculates the average magnetization of the system
double Lattice::calculateMagnetization() const {
    // Sums up all the +1s and -1s in the entire vector, starting with an initial sum of 0.0
    double sum = std::accumulate(spins.begin(), spins.end(), 0.0);
    // Returns the absolute value divided by total spins (N). 
    // Absolute value is used because a fully aligned lattice of -1s is physically identical to one of +1s.
    return std::abs(sum) / N; 
}