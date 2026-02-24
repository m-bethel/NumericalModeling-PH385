#include "PingPong.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

// Constructor
// Add this to PingPongBall.cpp

PingPongBall::PingPongBall(double mass, double diameter, double density,
                 double dragCoefficient, double x, double y, double z, 
                 double vx, double vy, double vz, double sx, double sy, double sz)
    : m_magnusCoefficient(0.04 * mass),
      m_fileName("pingpong_trajectory.out")
{
    // Properties
    m_mass = mass;
    m_diameter = diameter;
    m_dragCoefficient = dragCoefficient;
    m_density = density;
    

    // Initial position
    m_initialPosition[0] = x;
    m_initialPosition[1] = y;
    m_initialPosition[2] = z;
    
    // Initial velocity
    m_initialVelocity[0] = vx;
    m_initialVelocity[1] = vy;
    m_initialVelocity[2] = vz;
    
    // Initial spin
    m_initialSpin[0] = sx;
    m_initialSpin[1] = sy;
    m_initialSpin[2] = sz;
    
    // Copy to current state
    for (int i = 0; i < 3; i++)
    {
        m_position[i] = m_initialPosition[i];
        m_velocity[i] = m_initialVelocity[i];
        m_spin[i] = m_initialSpin[i];
    }
}

void PingPongBall::calculateAccelerationComponents(double accel[3])
{
    
    // Calculate velocity magnitude
    double v = sqrt(m_velocity[0] * m_velocity[0] +
                    m_velocity[1] * m_velocity[1] + 
                    m_velocity[2] * m_velocity[2]);
    
    // Calculate drag force
    double dragMagnitude = 0.5 * m_density * v * v * m_dragCoefficient * 
                          (M_PI * m_diameter * m_diameter / 4.0);
    
    double dragForce[3] = {0.0, 0.0, 0.0};
    if (v > 0.001)
    {
        for (int i = 0; i < 3; i++)
        {
            dragForce[i] = -(dragMagnitude * m_velocity[i] / v);
        }
    }
    
    // Gravitational force
    double gravityForce[3] = {0.0, 0.0, -m_mass * GRAVITY};
    
    // Magnus force (spin effect)
    double magnusForce[3] = {0.0, 0.0, 0.0};
    if (v > 0.001)
    {
        // Cross product: spin × velocity
        double crossProduct[3];
        crossProduct[0] = m_spin[1] * m_velocity[2] - m_spin[2] * m_velocity[1];
        crossProduct[1] = m_spin[2] * m_velocity[0] - m_spin[0] * m_velocity[2];
        crossProduct[2] = m_spin[0] * m_velocity[1] - m_spin[1] * m_velocity[0];
        
        for (int i = 0; i < 3; i++)
        {
            magnusForce[i] = m_magnusCoefficient * crossProduct[i]; 
        }
    }
    
    // Total acceleration
    for (int i = 0; i < 3; i++)
    {
        accel[i] = (dragForce[i] + gravityForce[i] + magnusForce[i]) / m_mass;
    }  
}

void PingPongBall::eulerStep()
{
    double accel[3];
    calculateAccelerationComponents(accel);
    
    // Euler's method integration
    for (int i = 0; i < 3; i++)
    {
        m_velocity[i] += accel[i] * TIME_STEP;
        m_position[i] += m_velocity[i] * TIME_STEP;
    }
}

void PingPongBall::simulate()
{
    ofstream outFile(m_fileName);
    outFile << "# x y z\n";
    
    double totalTime = 5.0; // seconds
    int maxSteps = static_cast<int>(totalTime / TIME_STEP);
    
    cout << "Simulating ping pong ball trajectory...\n";
    
    for (int step = 0; step < maxSteps; step++)
    {
        // Save current position
        outFile << m_position[0] << " " << m_position[1] << " " << m_position[2] << "\n";
        
        // Update state using Euler's method
        eulerStep();
        
        // Stop if ball hits ground
        if (m_position[2] < 0.0)
        {
            cout << "Ball hit ground at t = " << step * TIME_STEP << " seconds\n";
            cout << "Impact position: (" << m_position[0] << ", " 
                 << m_position[1] << ", " << m_position[2] << ") meters\n";
            break;
        }
    }
    
    outFile.close();
    cout << "Trajectory saved to " << m_fileName << "\n";
}

void PingPongBall::plotTrajectory()
{
    cout << "\nGenerating 3D plot with gnuplot...\n";
    
    // Create a gnuplot script file
    ofstream gnuplotScript("plot_script.gnu");
    gnuplotScript << "set title 'Ping Pong Ball Trajectory'\n";
    gnuplotScript << "set xlabel 'X (meters)'\n";
    gnuplotScript << "set ylabel 'Y (meters)'\n";
    gnuplotScript << "set zlabel 'Z (meters)'\n";
    gnuplotScript << "set grid\n";
    gnuplotScript << "splot '" << m_fileName << "' with lines linewidth 2 title 'Trajectory'\n";
    gnuplotScript << "pause -1 'Press any key to close'\n";
    gnuplotScript.close();
    
    // Execute gnuplot with the script
    int result = system("gnuplot plot_script.gnu");
    
    if (result != 0)
    {
        cout << "Error: Could not execute gnuplot. Make sure gnuplot is installed.\n";
        cout << "You can manually run: gnuplot plot_script.gnu\n";
    }
}

void PingPongBall::outputProperties() const
{
    cout << "Ping Pong Ball Properties:\n";
    cout << "Mass: " << m_mass << " kg\n";
    cout << "Diameter: " << m_diameter << " m\n";
    cout << "Drag Coefficient: " << m_dragCoefficient << "\n";
    cout << "Initial Position: (" << m_initialPosition[0] << ", " 
         << m_initialPosition[1] << ", " << m_initialPosition[2] << ") m\n";
    cout << "Initial Velocity: (" << m_initialVelocity[0] << ", " 
         << m_initialVelocity[1] << ", " << m_initialVelocity[2] << ") m/s\n";
    cout << "Initial Spin: (" << m_initialSpin[0] << ", " 
         << m_initialSpin[1] << ", " << m_initialSpin[2] << ") rad/s\n";
    cout << "Time Step: " << TIME_STEP << " s\n\n";
}