/*
Class Object

Create Celectial objects to use in the modeling of a Solar System
specifically for a 3 body problem

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

struct Vector3D {
    double x, y, z;
    
    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    
    Vector3D operator+(const Vector3D& v) const {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }
    
    Vector3D operator-(const Vector3D& v) const {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }
    
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    Vector3D operator/(double scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }
    
    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    
    Vector3D normalized() const {
        double mag = magnitude();
        return (mag > 0) ? (*this / mag) : Vector3D(0, 0, 0);
    }
};

#endif

#ifndef OBJECT_H
#define OBJECT_H

#include <string>

class Object
{
private:
    std::string name;
    double mass;
    Vector3D position;
    Vector3D velocity;
    Vector3D acceleration;
    Vector3D prevAcceleration;
    
public:
    Object(const std::string& name, double mass, const Vector3D& pos, const Vector3D& vel);
    
    // Getters
    std::string getName() const;
    double getMass() const;
    Vector3D getPosition() const;
    Vector3D getVelocity() const;
    Vector3D getAcceleration() const;
    
    // Setters
    void setAcceleration(const Vector3D& acc);
    
    // Verlet integration methods
    void updatePosition(double dt);
    void updateVelocity(double dt);
};

#endif //OBJECT_H