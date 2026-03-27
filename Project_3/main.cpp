/*
Main file

Simulates a three body orbit of the Solar System using verlet integration, originally
designed to do a 2 dimentional orbit simulation, but made for 3 dimensions.
Make sure to input/redefine the 3rd dimension if you would like that
 (set to zero currently)

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/

#include <cmath>
#include <iostream>
#include "system.h"
#include "object.h"

using namespace std;

int main() {
    // Simulation config
    double simulationTime = 11.2;        // Simulation duration in Earth years (roughly 1 Jupiter orbit)
    const double DAY = 1.0 / 365.25;     // One day represented as a fraction of a year
    double timeStep = 1 * DAY;           // 1-day increments for the solver

    // Experimental mass multipliers to test system stability
    double jupiterMassMult = 1000; 
    double earthMassMult = 1;
    double sunMassMult = 1;

    cout << "=== Solar System Simulation: Sun, Earth, Jupiter ===" << endl;

    System solarSystem;

    // SUN: Mass=1 (Solar Mass), Position=(0,0,0), Velocity=(0,0,0)
    Object sun("Sun", 1.0 * sunMassMult, Vector3D(0,0,0), Vector3D(0,0,0));
    solarSystem.addObject(sun);

    // EARTH: Mass~3e-6, 1 AU from Sun, Velocity=2*pi (1 orbit per year)
    Object earth("Earth", 3.0027e-6 * earthMassMult, Vector3D(1,0,0), Vector3D(0, 2*M_PI, 0));
    solarSystem.addObject(earth);

    // JUPITER: Mass~1e-3, ~5.2 AU from Sun, lower orbital velocity
    Object jupiter("Jupiter", 0.000954 * jupiterMassMult, Vector3D(-5.2, 0, 0), Vector3D(0, -0.876897*M_PI, 0));
    solarSystem.addObject(jupiter);

    cout << "Starting simulation for " << simulationTime << " years..." << endl;
    solarSystem.simulate(simulationTime, timeStep, "orbit_data.csv");
    cout << "Done. Results saved to orbit_data.csv" << endl;

    return 0;
}