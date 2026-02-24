/*
 * program.h
 * 
 * Header file for simulating a driven damped pendulum system using RK4 integration.
 * This simulation models chaotic behavior and generates Poincaré sections.
 * 
 * NUMERICAL LIMITATIONS AND ERROR ACCUMULATION:
 * 
 * 1. RK4 Integration Error:
 *    - Local truncation error: O(h^5) per step, where h = TIME_STEP
 *    - Global truncation error: O(h^4) over the entire simulation
 *    - For TIME_STEP = 0.01, errors accumulate over long simulations
 * 
 * 2. Floating-Point Precision:
 *    - Using double precision (64-bit) provides ~15-17 significant digits
 *    - Catastrophic cancellation can occur in trigonometric functions near π
 *    - Angle normalization reduces but doesn't eliminate floating-point drift
 * 
 * 3. Long-Term Simulation Issues:
 *    - For TIME = 72000s (20 hours), accumulated error can affect phase space structure
 *    - Chaotic systems are sensitive to initial conditions (butterfly effect)
 *    - Small numerical errors grow exponentially with Lyapunov exponent
 * 
 * 4. Poincaré Section Sampling:
 *    - Tolerance for time matching introduces systematic sampling error
 *    - Points may not be sampled at exactly t = 2πn/ΩD
 *    - Increasing TIME_STEP increases sampling uncertainty
 * 
 * 5. Recommendations:
 *    - Keep TIME_STEP ≤ 0.01 for accuracy
 *    - Verify results by comparing different TIME_STEP values
 *    - For very long simulations (>100000s), consider higher-order methods
 */

#ifndef PROGRAM_H
#define PROGRAM_H

#include <string>
#include <vector>

class Pendulum
{
private:
    // Physical Properties - Parameters defining the pendulum system
    double m_mass;                  // Mass of pendulum bob (kg)
    double m_length;                // Length of pendulum rod (m)
    double m_dampingCoefficient;    // Damping coefficient (dimensionless)
    double m_initialAngle;          // Initial angular displacement (radians)
    double m_initialOmega;          // Initial angular velocity (rad/s)
    double m_drivingForceMag;       // Amplitude of driving force (dimensionless)
    double m_drivingForceFreq;      // Angular frequency ΩD of driving force (rad/s)

    // Physical Constants
    static constexpr double GRAVITY = 9.81;      // Gravitational acceleration (m/s^2)
    static constexpr double TIME = 72000.0;      // Total simulation time (seconds)
    static constexpr double TIME_STEP = 0.01;    // Time step for RK4 integration (seconds)

    // Storage for full trajectory data
    std::vector<double> m_theta;    // Angular position at each time step (radians)
    std::vector<double> m_omega;    // Angular velocity at each time step (rad/s)
    std::vector<double> m_time;     // Time values (seconds)

    // Storage for Poincaré section points
    // These are sampled at t = 2πn/ΩD (when driving force is in phase)
    std::vector<double> m_poincare_theta;  // Theta values at Poincaré section
    std::vector<double> m_poincare_omega;  // Omega values at Poincaré section

    // Output filename for trajectory data
    std::string m_fileName;

    // Private helper methods
    
    /**
     * Calculate angular acceleration α for the driven damped pendulum
     * Equation: α = -(g/L)sin(θ) - q·ω + F·cos(ΩD·t)
     * 
     * @param theta Current angular position (radians)
     * @param omega Current angular velocity (rad/s)
     * @param t Current time (seconds)
     * @return Angular acceleration (rad/s^2)
     */
    double calculateAlpha(double theta, double omega, double t);
    
    /**
     * Perform one step of 4th-order Runge-Kutta integration
     * Updates theta and omega to their values at t + TIME_STEP
     * 
     * @param theta Reference to angular position (modified in place)
     * @param omega Reference to angular velocity (modified in place)
     * @param t Current time (seconds)
     */
    void RK4(double& theta, double& omega, double t);

public:
    /**
     * Constructor - Initialize pendulum with specified parameters
     * 
     * @param mass Mass of pendulum bob (kg)
     * @param length Length of pendulum rod (m)
     * @param dampingCoefficient Damping coefficient (dimensionless)
     * @param initialAngle Initial angular displacement (radians)
     * @param initialOmega Initial angular velocity (rad/s)
     * @param drivingForceMag Amplitude of driving force (dimensionless)
     * @param drivingForceFreq Angular frequency ΩD of driving force (rad/s)
     */
    Pendulum(double mass, double length,
             double dampingCoefficient, double initialAngle,
             double initialOmega, double drivingForceMag, 
             double drivingForceFreq);

    /**
     * Run the full simulation using RK4 integration
     * Stores full trajectory and Poincaré section points
     */
    void simulate();
    
    /**
     * Write trajectory data to file and generate plots using gnuplot
     * Creates three plots:
     *   - Poincaré section (omega vs theta at t = 2πn/ΩD)
     *   - Full phase space (complete trajectory)
     *   - Time series (theta and omega vs time)
     */
    void plotTrajectory();
    
    /**
     * Display pendulum properties and simulation parameters to console
     */
    void outputProperties() const;
};

#endif // PROGRAM_H