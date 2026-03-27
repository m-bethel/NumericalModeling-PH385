/*
system.h file

Header file for the system.cpp, providing the input for the actualy simulation

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/

#ifndef SYSTEM_h
#define SYSTEM_h

#include "object.h"
#include <string>
#include <vector>

// Manages the collection of objects and the simulation clock
class System {
private:
    std::vector<Object> objects; // List of all planets/stars in simulation
    double time;                 // Current simulation time
public:
    System();

    void addObject(const Object& obj);           // Add a body to the system
    Vector3D calculateAcceleration(size_t i);    // Calculate net gravitational pull on object i
    void step(double dt);                        // Advance system by one time step
    void simulate(double duration, double dt, const std::string& outputFile); // Run loop and save to CSV
    void printState();                           // Output positions to console
};

#endif