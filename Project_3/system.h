#ifndef SYSTEM_h
#define SYSTEM_h

#include "object.h"
#include <string>
#include <vector>

class System
{
private:
    std::vector<Object> objects;
    double time;
public:
    System();

    void addObject(const Object& obj);
    Vector3D calculateAcceleration(size_t i);
    void step(double dt);
    void simulate(double duration, double dt, const std::string& outputFile);
    void printState();


};

#endif // SYSTEM_H