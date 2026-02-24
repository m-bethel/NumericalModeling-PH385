/*
Main file

Simulates a three body orbit of the Solar System using verlet integration, originally
designed to do a 2 dimentional orbit simulation, but made for 3 dimensions.
Make sure to input/redefine the 3rd dimension if you would like that
 (set to zero currently)

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/

#include <iostream>
#include "system.h"
#include "object.h"

using namespace std;

int main()
{
    double mass;
    double mass_Sun;
    double mass_Jupiter;
    double mass_earth;
    string name;
    string sun = "Sun";
    string jupiter = "Jupiter";
    string earth = "Earth";
    double position;
    double x,y,z;
    double velocity;
    double vx,vy,vz;


    cout << "=== Three Body Simulation of Earth, Jupiter, and the Sun ===" << endl;

    // ===== Create Pendulum object =====
    // Constructor initializes the pendulum with the specified parameters
    Object sun("Sun", 1, Vector3D(0,0,0), Vector3D(0,0,0));
    solarSystem.addObject(sun);

    Object earth("Earth", 3.0027e-6, Vector3D(1,0,0), Vector3D(0,2*M_PI,0));
    solarSystem.addObject(earth);

    // Simulation parameters
    double timeStep = 1 * DAY;          // 1 day time steps
    double simulationTime = 365 * DAY;  // Simulate 1 year
    
    std::cout << "Starting Earth orbit simulation...\n";
    std::cout << "Time step: " << timeStep/DAY << " days\n";
    std::cout << "Duration: " << simulationTime/DAY << " days (1 year)\n\n";
    
    // Run simulation
    solarSystem.simulate(simulationTime, timeStep, "orbit_data.csv");
    
    std::cout << "\nYou can plot the orbit using Python:\n";
    std::cout << "import pandas as pd\n";
    std::cout << "import matplotlib.pyplot as plt\n";
    std::cout << "data = pd.read_csv('orbit_data.csv')\n";
    std::cout << "earth = data[data['object'] == 'Earth']\n";
    std::cout << "plt.plot(earth['x'], earth['y'])\n";
    std::cout << "plt.xlabel('x (AU)')\n";
    std::cout << "plt.ylabel('y (AU)')\n";
    std::cout << "plt.title('Earth Orbit')\n";
    std::cout << "plt.axis('equal')\n";
    std::cout << "plt.grid(True)\n";
    std::cout << "plt.show()\n";
    
    return 0;
    return 0;
}