/*
Class Object

Create Celectial objects to use in the modeling of a Solar System
specifically for a 3 body problem

Also include getters and setters for the system

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/



#include "object.h"
#include <iostream>
#include <vector>

using namespace std;

// Initialize object state; accelerations start at zero until first calculateAcceleration call
Object::Object(const string& name, double mass, const Vector3D& pos, const Vector3D& vel)
    : name(name), mass(mass), position(pos), velocity(vel),
    acceleration(0,0,0), prevAcceleration(0,0,0) {}

string Object::getName() const { return name; }
double Object::getMass() const { return mass; }
Vector3D Object::getPosition() const { return position; }
Vector3D Object::getVelocity() const { return velocity; }
Vector3D Object::getAcceleration() const { return acceleration; }

// Before setting new acceleration, we archive the current one as 'prevAcceleration'
// Velocity Verlet requires both a(t) and a(t+dt) to update velocity
void Object::setAcceleration(const Vector3D& acc) { 
    prevAcceleration = acceleration;
    acceleration = acc; 
}

// Position update: r(t+dt) = r(t) + v(t)dt + 0.5 * a(t)dt^2
void Object::updatePosition(double dt) {
    position = position + velocity * dt + acceleration * (0.5 * dt * dt);
}

// Velocity update: v(t+dt) = v(t) + 0.5 * [a(t) + a(t+dt)]dt
void Object::updateVelocity(double dt) {
    velocity = velocity + (acceleration + prevAcceleration) * (0.5 * dt);
}