#ifndef PINGPONGBALL_H
#define PINGPONGBALL_H

#include <string>

class PingPongBall
{
private:
    // Physical properties
    double m_mass;
    double m_diameter;
    double m_dragCoefficient;
    double m_magnusCoefficient;
    double m_density;    // kg/m^3
    
    // Initial state
    double m_initialPosition[3];
    double m_initialVelocity[3];
    double m_initialSpin[3];
    
    // Current state
    double m_position[3];
    double m_velocity[3];
    double m_spin[3];
    
    // Output file
    std::string m_fileName;
    
    // Physical constants
    static constexpr double GRAVITY = 9.81;        // m/s^2
    static constexpr double TIME_STEP = 0.0001;    // seconds
    
    // Private helper methods
    void calculateAccelerationComponents(double accel[3]);
    void eulerStep();

public:
    // Constructor
    PingPongBall(double mass, double diameter, double density,
                 double dragCoefficient, double x, double y, double z, 
                 double vx, double vy, double vz, double sx, double sy, double sz);
    
    // Public methods
    void simulate();
    void plotTrajectory();
    void outputProperties() const;
    
    // Getters (optional, for accessing state if needed)
    const double* getPosition() const { return m_position; }
    const double* getVelocity() const { return m_velocity; }
};

#endif // PINGPONGBALL_H