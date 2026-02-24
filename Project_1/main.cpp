#include <iostream>
#include <limits>
#include "PingPong.h"

using namespace std;

/*
 * Ping Pong Ball Trajectory Simulator
 * 
 * VALID INPUT RANGES:
 * - Mass: 2.5-2.9 grams (regulation: 2.7g)
 * - Diameter: 3.8-4.2 cm (regulation: 4.0cm)
 * - Initial height: 0.1-10 meters (below 0.1m may have collision detection issues)
 * - Velocity: < 50 m/s (higher speeds exceed model validity)
 * - Spin: < 500 rad/s (model accuracy decreases above this)
 * 
 * LIMITATIONS:
 * - Does not simulate table/paddle collisions
 * - Assumes flight in still air at sea level
 * - Results not avalible for flight times > 5 seconds unless 
 * - changed manually in PingPong.cpp file variable totalTime
 */

int main() 
{
    // Initialize variable to pass to the header file
    double mass;
    double diameter;
    double density;
    double dragCoefficient;
    double x, y, z;
    double vx, vy, vz;
    double sx, sy, sz;
    char choice;

    // User prompts for initial conditions
    // This asks user if they would like to use default conditons provided
    cout << "Would you like to accept default initial conditions? [y,n]\n";
    cout << "Mass: 2.7 grams \nDiameter: 4.0 centimeters \nAir density: 1.27 kg/m^3 \nDrag Coefficient: 0.5 \nPosition: (0, 0, 5) m \nVelocity: (4 4 10) m/s \nSpin Vector (-50, -100, 100)" << endl;
    
    cin >> choice;

    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    
    // If statement to check if you would like to use default
    // conditions or input your own
    if (choice  == 'y' || choice == 'Y')
    {
        cout << "Using default conditions" << endl;
        mass = 0.0027, diameter = 0.04, dragCoefficient = 0.5,
        density = 1.27, x = 0, y = 0, z = 5, vx = 4, vy = 4,
        vz = 10, sx = -50, sy = -100, sz = 100;
    } else
    {
        cout << "Enter Initial Conditions: \n" << endl;
        cout << "Mass and Diameter (mass diameter): (meters)";
        cin >> mass >> diameter;
        cout << "Density and Drag Coefficient: (Density Drag Coefficent) ";
        cin >> density >> dragCoefficient;
        cout << "Position (x y z): (meters) ";
        cin >> x >> y >> z;
        cout << "Velocity (x y z): (meters/s) ";
        cin >> vx >> vy >> vz;
        cout << "Initial Spin Components (x y z):";
        cin >> sx >> sy >> sz;
    }

    // Creates object ball1
    PingPongBall ball1( mass, diameter, density,
                        dragCoefficient, x, y, z, 
                        vx, vy, vz, sx, sy, sz);

    // Checks for negative z value, cuz then it cause issues
    if (z < 0)
    {
        cout << "Error: z value cannot be negative. Using default height of 5.0 meters.\n";
        z= 5.0;
    }
    
    ball1.outputProperties();
    ball1.simulate();
    ball1.plotTrajectory();
    
    return 0;
}