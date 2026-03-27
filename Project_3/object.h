/*
object.h file

Create Celectial objects to use in the modeling of a Solar System
specifically for a 3 body problem

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 1/26/2026
*/

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

// Helper structure to handle 3D vector math (Positions, Velocities, Accelerations)
struct Vector3D {
    double x, y, z;
    
    // Default constructor
    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    
    // Operator overloads to allow natural math: v1 + v2, v1 * scalar, etc.
    Vector3D operator+(const Vector3D& v) const { return Vector3D(x + v.x, y + v.y, z + v.z); }
    Vector3D operator-(const Vector3D& v) const { return Vector3D(x - v.x, y - v.y, z - v.z); }
    Vector3D operator*(double scalar) const { return Vector3D(x * scalar, y * scalar, z * scalar); }
    Vector3D operator/(double scalar) const { return Vector3D(x / scalar, y / scalar, z / scalar); }
    
    // Calculates the Euclidean distance from origin
    double magnitude() const { return std::sqrt(x*x + y*y + z*z); }
    
    // Returns a vector in the same direction but with a length of 1
    Vector3D normalized() const {
        double mag = magnitude();
        return (mag > 0) ? (*this / mag) : Vector3D(0, 0, 0);
    }
};

#endif

#ifndef OBJECT_H
#define OBJECT_H

#include <string>

// Represents a single celestial body (Planet/Star)
class Object {
private:
    std::string name;
    double mass;
    Vector3D position;
    Vector3D velocity;
    Vector3D acceleration;      // Current acceleration a(t)
    Vector3D prevAcceleration;  // Acceleration from the previous step a(t-dt)
    
public:
    Object(const std::string& name, double mass, const Vector3D& pos, const Vector3D& vel);

    // Getters for physical properties
    std::string getName() const;
    double getMass() const;
    Vector3D getPosition() const;
    Vector3D getVelocity() const;
    Vector3D getAcceleration() const;

    // Updates acceleration and stores the old one (Crucial for Velocity Verlet)
    void setAcceleration(const Vector3D& acc);
    
    // Integration steps for position and velocity
    void updatePosition(double dt);
    void updateVelocity(double dt);
};

#endif