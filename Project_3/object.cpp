#include "object.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

Object::Object(const string& name, double mass, const Vector3D& pos, const Vector3D& vel)
    : name(name), mass(mass), position(pos), velocity(vel),
    acceleration(0,0,0), prevAcceleration(0,0,0) {}


    string Object::getName() const
    {
        return name;
    }
    double Object::getMass() const
    {
        return mass;
    }
    Vector3D Object::getPosition() const 
    { 
    return position; 
    }

    Vector3D Object::getVelocity() const 
    { 
        return velocity; 
    }

    Vector3D Object::getAcceleration() const 
    { 
        return acceleration; 
    }

    void Object::setAcceleration(const Vector3D& acc) 
    { 
        prevAcceleration = acceleration;
        acceleration = acc; 
    }

    void Object::updatePosition(double dt) 
    {
        position = position + velocity * dt + acceleration * (0.5 * dt * dt);
    }

    void Object::updateVelocity(double dt) 
    {
        velocity = velocity + (acceleration + prevAcceleration) * (0.5 * dt);
    }

