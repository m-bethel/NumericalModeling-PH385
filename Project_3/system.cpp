#include "system.h"
#include <iostream>
#include <fstream>

// Physical Constants
const double G = 6.67430e-11;
const double AU = 1.496e11;

System::System() : time(0.0) {}
void System::addObject(const Object& obj)
{
    objects.push_back(obj);
}

Vector3D System::calculateAcceleration(size_t i) {
    Vector3D totalAccel(0, 0, 0);
    
    for (size_t j = 0; j < objects.size(); j++) {
        if (i == j) continue;
        
        Vector3D r = objects[j].getPosition() - objects[i].getPosition();
        double distance = r.magnitude();
        
        if (distance > 0) {
            // F = G * m1 * m2 / r^2, a = F/m1 = G * m2 / r^2
            double forceMagnitude = G * objects[j].getMass() / (distance * distance);
            Vector3D forceDirection = r.normalized();
            totalAccel = totalAccel + forceDirection * forceMagnitude;
        }
    }
    
    return totalAccel;
}

void System::step(double dt) {
    // Calculate accelerations for all objects
    std::vector<Vector3D> newAccelerations;
    for (size_t i = 0; i < objects.size(); i++) {
        newAccelerations.push_back(calculateAcceleration(i));
    }
    
    // Update positions using current velocities and accelerations
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i].updatePosition(dt);
    }
    
    // Set new accelerations
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i].setAcceleration(newAccelerations[i]);
    }
    
    // Update velocities using average of old and new accelerations
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i].updateVelocity(dt);
    }
    
    time += dt;
}
