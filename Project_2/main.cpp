/*
 * main.cpp
 * 
 * Main driver program for simulating a driven damped pendulum system.
 * This program demonstrates chaotic behavior and creates Poincaré sections
 * for phase space analysis, as described in Exercise 3.12.
 * 
 * The driven damped pendulum equation:
 *   d²θ/dt² = -(g/L)·sin(θ) - q·(dθ/dt) + F·cos(ΩD·t)
 * 
 * Where:
 *   θ = angular displacement from vertical
 *   g = gravitational acceleration
 *   L = pendulum length
 *   q = damping coefficient
 *   F = driving force amplitude
 *   ΩD = driving force angular frequency
 */

#include <iostream>
#include "program.h"

using namespace std;

int main()
{
    // Declare variables for pendulum parameters
    double mass;                // Mass of pendulum bob (kg)
    double length;              // Length of pendulum rod (m)
    double dampingCoefficient;  // Damping coefficient (dimensionless)
    double initialAngle;        // Initial angular displacement (radians)
    double initialOmega;        // Initial angular velocity (rad/s)
    double drivingForceMag;     // Driving force amplitude (dimensionless)
    double drivingForceFreq;    // Driving force angular frequency ΩD (rad/s)
    char choice;                // User's choice for default vs. custom parameters

    // ===== Display default parameters to user =====
    cout << "=== Driven Damped Pendulum Simulation ===" << endl;
    cout << "Would you like to accept default initial conditions? [y,n]\n" << endl;
    cout << "Default Parameters:" << endl;
    cout << "  Mass: 1.0 kg" << endl;
    cout << "  Length: 9.8 m" << endl;
    cout << "  Damping Coefficient: 0.5" << endl;
    cout << "  Initial Angle: 0.2 radians" << endl;
    cout << "  Initial Angular Velocity: 0.0 rad/s" << endl;
    cout << "  Driving Force Magnitude: 1.2" << endl;
    cout << "  Driving Force Angular Frequency: 2/3 rad/s" << endl;
    cout << "\nYour choice [y/n]: ";
    
    cin >> choice;
    
    // ===== Handle user's choice =====
    if (choice == 'y' || choice == 'Y')
    {
        // Use default parameters for chaotic behavior demonstration
        cout << "\nUsing default conditions..." << endl;
        
        mass = 1.0;                 // 1 kg bob
        length = 9.8;               // ~9.8 m length (convenient for g/L ≈ 1)
        dampingCoefficient = 0.5;   // Moderate damping
        initialAngle = 0.2;         // Small initial displacement (0.2 rad ≈ 11.5°)
        initialOmega = 0.0;         // Start from rest
        drivingForceMag = 1.2;      // Strong driving force (produces chaos)
        drivingForceFreq = 2.0/3.0; // Driving angular frequency ΩD = 2/3 rad/s
    }
    else
    {
        // Prompt user for custom parameters
        cout << "\n=== Enter Custom Initial Conditions ===" << endl;
        
        cout << "Mass (kg): ";
        cin >> mass;
        
        cout << "Length (meters): ";
        cin >> length;
        
        cout << "Damping Coefficient: ";
        cin >> dampingCoefficient;
        
        cout << "Initial Angle (radians): ";
        cin >> initialAngle;
        
        cout << "Initial Angular Velocity (rad/s): ";
        cin >> initialOmega;
        
        cout << "Driving Force Magnitude: ";
        cin >> drivingForceMag;
        
        cout << "Driving Force Angular Frequency ΩD (rad/s): ";
        cin >> drivingForceFreq;
        
        cout << "\nUsing custom parameters..." << endl;
    }

    // ===== Create Pendulum object =====
    // Constructor initializes the pendulum with the specified parameters
    Pendulum object(mass, length, dampingCoefficient,
                    initialAngle, initialOmega, 
                    drivingForceMag, drivingForceFreq);

    // ===== Run the simulation =====
    // Step 1: Display the system properties
    object.outputProperties();
    
    // Step 2: Run RK4 integration and collect data
    // This will take some time for long simulations (e.g., 72000 seconds)
    object.simulate();
    
    // Step 3: Write data files and generate plots
    object.plotTrajectory();
    
    cout << "\n=== Simulation Complete ===" << endl;
    cout << "Check the generated PNG files to see the results." << endl;
    cout << "The Poincaré section reveals the chaotic attractor structure." << endl;
    
    return 0;
}