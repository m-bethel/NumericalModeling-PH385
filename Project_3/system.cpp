/*
system.cpp-file

This file instantiates the solarSystem, and adds other bodies to it.
Simulates the 3-body problem using Verlet Integration.

VERLET INTEGRATION ANALYSIS:
----------------------------
BONUSES:
1. Symplectic/Energy Conserving: Unlike Euler, Verlet doesn't "drift" in energy. 
   Orbits remain stable and don't spiral inward or outward indefinitely.
2. Time-Reversible: The physics remains consistent if you run time backward.
3. Computational Efficiency: Only requires one force calculation per step.

LIMITATIONS:
1. Handling Velocity-Dependent Forces: Standard Verlet is difficult to use 
   if forces depend on velocity (like air resistance).
2. Constant Time-Step: It performs best with a fixed dt. Variable time-steps 
   break the stability and reversibility of the algorithm.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/


#include <cmath>
#include "system.h"
#include <iostream>
#include <fstream>

// G = 4*pi^2 is the constant when units are Solar Masses, AU, and Years
const double G = 4.0 * M_PI * M_PI; 

System::System() : time(0.0) {}

void System::addObject(const Object& obj) {
    objects.push_back(obj);
}

// Sums the gravitational pull from every other object in the system
Vector3D System::calculateAcceleration(size_t i) {
    Vector3D totalAccel(0, 0, 0);
    for (size_t j = 0; j < objects.size(); j++) {
        if (i == j) continue; // Don't calculate gravity on yourself
        
        Vector3D r = objects[j].getPosition() - objects[i].getPosition();
        double distance = r.magnitude();
        
        if (distance > 0) {
            // Newton's Universal Gravitation: a = G * m_other / r^2
            double forceMagnitude = G * objects[j].getMass() / (distance * distance);
            Vector3D forceDirection = r.normalized();
            totalAccel = totalAccel + forceDirection * forceMagnitude;
        }
    }
    return totalAccel;
}

// The Velocity Verlet sequence:
void System::step(double dt) {
    // 1. Update positions of all objects using current v and a
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i].updatePosition(dt);
    }
    
    // 2. Calculate new accelerations based on the new positions
    std::vector<Vector3D> newAccelerations;
    for (size_t i = 0; i < objects.size(); i++) {
        newAccelerations.push_back(calculateAcceleration(i));
    }
    
    // 3. Update acceleration states (moves current a to prev_a)
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i].setAcceleration(newAccelerations[i]);
    }
    
    // 4. Update velocities using the average of old and new acceleration
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i].updateVelocity(dt);
    }
    
    time += dt;
}

// Runs the main loop and logs data for Python visualization
void System::simulate(double duration, double dt, const std::string& outputFile) {
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Could not open file " << outputFile << "\n";
        return;
    }

    outFile << "time,object,x,y,z\n"; // Header for Pandas

    while (time < duration) {
        // Log current positions before taking the next step
        for (size_t i = 0; i < objects.size(); ++i) {
            Vector3D pos = objects[i].getPosition();
            outFile << time << "," << objects[i].getName() << "," 
                    << pos.x << "," << pos.y << "," << pos.z << "\n";    
        }
        step(dt);
    }
    outFile.close();
}