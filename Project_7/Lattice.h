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
*/
#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

class Lattice {
private:
    int L;        // Dimension length
    int N;        // Total spins (L^3)
    std::vector<int> spins;
    int get_idx(int x, int y, int z) const;

public:
    // Constructor
    Lattice(int size);

    // Initialize with random spins
    void initializeRandom();
    void initializeAligned();

    void performSweep(double T, double H);

    // Get the sum of neighbors for a specific spin
    int getNeighborSum(int x, int y, int z) const;

    // Flip a spin at a location
    void flipSpin(int x, int y, int z);

    // Get value of a spin
    int getSpin(int x, int y, int z) const;

    // Calculate total magnetization
    double calculateMagnetization() const;

    int getSize() const { return L; }
};

#endif