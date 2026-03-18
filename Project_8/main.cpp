/*
Main File

(README)
Models the molecular dynamics of particles as system energy is increased
and decreased. It will create an animation and graph of the particles moving and 
a Q-T graph, respectively

Output is a png that represents Tempuerature vs. Total Internal Energy 

E = KE + PE

To adjust parameters, adjust everything here, specifically, number of particles,
how much energy to add in each step and the time steps as well.

You can adjust step count and a few other paramenters in the Simulate.cpp file
But it I don't recommend it without understanding the code throughly.
Check the Simulation.cpp (README) for more paramter info.

If you adjust the parameters to the point where it's not very stable, you can always add more particles
The more particles, the more stable and nicer it is. You kinda just trade computational time for stability
Feel free to use whatever number, the program will place them at the center in a grid just fine.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/2/2026
*/
#include <iostream>     // cout, cin
#include <limits>       // numeric_limits for input buffer clearing
#include "Simulate.h"   // Simulate class - manages the full simulation and animation
#include "Particle.h"   // Particle class - represents a particle in the simulation
#include <ctime>

using namespace std;

int main()
{
    int particles = 400;
    double deltaEnergy = 0.0055;
    int gridSize = 50;
    srand(42); // time(0) change me back for random seeds

    Simulate simulation (particles, deltaEnergy, gridSize);

    simulation.simulate(0.001);
    simulation.output();
    
    return 0;
}