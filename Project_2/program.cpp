/*
 * program.cpp
 * 
 * Implementation of the Pendulum class for simulating a driven damped pendulum
 * using the 4th-order Runge-Kutta (RK4) method. This simulation models chaotic
 * behavior and generates Poincaré sections for phase space analysis.
 */

#include "program.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

/**
 * Constructor - Initialize the pendulum with given parameters
 * Reserves memory for trajectory vectors to improve performance
 */
Pendulum::Pendulum(double mass, double length,
                   double dampingCoefficient, double initialAngle,
                   double initialOmega, double drivingForceMag, double drivingForceFreq)
    : m_mass(mass),
      m_length(length),
      m_dampingCoefficient(dampingCoefficient),
      m_initialAngle(initialAngle),
      m_initialOmega(initialOmega),
      m_drivingForceMag(drivingForceMag),
      m_drivingForceFreq(drivingForceFreq),
      m_fileName("pendulum_trajectory.out")
{
    // Pre-allocate memory for trajectory vectors to avoid repeated reallocations
    int numSteps = static_cast<int>(TIME / TIME_STEP);
    m_theta.reserve(numSteps);
    m_omega.reserve(numSteps);
    m_time.reserve(numSteps);
    
    // Reserve space for Poincaré section (estimate ~1000 points)
    m_poincare_theta.reserve(1000);
    m_poincare_omega.reserve(1000);
}

/**
 * Calculate the angular acceleration for a driven damped pendulum
 * 
 * Equation of motion:
 * α = d²θ/dt² = -(g/L)·sin(θ) - q·(dθ/dt) + F·cos(ΩD·t)
 * 
 * Where:
 *   - First term: gravitational restoring torque
 *   - Second term: damping (friction) torque
 *   - Third term: external driving force
 * 
 * @param theta Current angular position (radians)
 * @param omega Current angular velocity (rad/s)
 * @param t Current time (seconds)
 * @return Angular acceleration α (rad/s²)
 */
double Pendulum::calculateAlpha(double theta, double omega, double t)
{
    // The driving force angular frequency ΩD is already in rad/s
    // So we use it directly without multiplying by 2π
    double angularFreq = m_drivingForceFreq;  // ΩD in rad/s
    
    // Gravitational torque term: -(g/L)·sin(θ)
    double gravityTerm = -(GRAVITY / m_length) * sin(theta);
    
    // Damping torque term: -q·ω
    double dampingTerm = -m_dampingCoefficient * omega;
    
    // Driving force term: F·cos(ΩD·t)
    double drivingTerm = m_drivingForceMag * cos(angularFreq * t);
    
    // Total angular acceleration
    return gravityTerm + dampingTerm + drivingTerm;
}

/**
 * Perform one step of 4th-order Runge-Kutta (RK4) integration
 * 
 * RK4 is a numerical method for solving ordinary differential equations (ODEs).
 * For our coupled system:
 *   dθ/dt = ω
 *   dω/dt = α(θ, ω, t)
 * 
 * RK4 computes four estimates (k1, k2, k3, k4) and takes a weighted average
 * to achieve O(h⁴) accuracy where h is the time step.
 * 
 * The method:
 *   k1 = f(t, y)              [slope at beginning]
 *   k2 = f(t + h/2, y + k1·h/2)  [slope at midpoint using k1]
 *   k3 = f(t + h/2, y + k2·h/2)  [slope at midpoint using k2]
 *   k4 = f(t + h, y + k3·h)      [slope at end using k3]
 *   y_new = y + (h/6)·(k1 + 2k2 + 2k3 + k4)
 * 
 * @param theta Reference to current angle (updated to new value)
 * @param omega Reference to current angular velocity (updated to new value)
 * @param t Current time
 */
void Pendulum::RK4(double& theta, double& omega, double t)
{
    double dt = TIME_STEP;

    // ===== K1: Evaluate at current point (t, θ, ω) =====
    double k1_theta = omega;  // dθ/dt = ω
    double k1_omega = calculateAlpha(theta, omega, t);  // dω/dt = α

    // ===== K2: Evaluate at midpoint using k1 (t + dt/2) =====
    double k2_theta = omega + 0.5 * dt * k1_omega;
    double k2_omega = calculateAlpha(theta + 0.5 * dt * k1_theta,
                                     omega + 0.5 * dt * k1_omega, 
                                     t + 0.5 * dt);

    // ===== K3: Evaluate at midpoint using k2 (t + dt/2) =====
    double k3_theta = omega + 0.5 * dt * k2_omega;
    double k3_omega = calculateAlpha(theta + 0.5 * dt * k2_theta,
                                     omega + 0.5 * dt * k2_omega, 
                                     t + 0.5 * dt);

    // ===== K4: Evaluate at endpoint using k3 (t + dt) =====
    // CRITICAL: Use full dt here, not 0.5*dt!
    double k4_theta = omega + dt * k3_omega;
    double k4_omega = calculateAlpha(theta + dt * k3_theta,
                                     omega + dt * k3_omega, 
                                     t + dt);

    // ===== Update theta and omega using weighted average =====
    // Weights: k1 and k4 get 1/6 each, k2 and k3 get 1/3 each
    theta += (dt / 6.0) * (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    omega += (dt / 6.0) * (k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega);
    
    // ===== Normalize theta to [-π, π] =====
    // This prevents theta from growing without bound and reduces
    // floating-point errors in trigonometric functions
    while (theta > M_PI) theta -= 2.0 * M_PI;
    while (theta < -M_PI) theta += 2.0 * M_PI;
}

/**
 * Run the complete pendulum simulation
 * 
 * This method:
 * 1. Initializes the system with given initial conditions
 * 2. Runs RK4 integration for the specified TIME duration
 * 3. Stores full trajectory data at every time step
 * 4. Collects Poincaré section points at t = 2πn/ΩD
 * 
 * The Poincaré section is created by sampling the phase space only when
 * the driving force is at the same phase, revealing the underlying structure
 * of the chaotic attractor.
 */
void Pendulum::simulate()
{
    cout << "\nRunning simulation..." << endl;
    
    // Initialize state variables
    double theta = m_initialAngle;
    double omega = m_initialOmega;
    double t = 0.0;
    
    // Calculate the period of the driving force: T = 2π/ΩD
    double angularFreq = m_drivingForceFreq;  // ΩD is already in rad/s
    double drivingPeriod = 2.0 * M_PI / angularFreq;

    // Time of next Poincaré section sample
    double nextSampleTime = drivingPeriod;
    
    // Track whether we've passed the sampling time in this iteration
    double prevTime = 0.0;

    // Store initial state
    m_theta.push_back(theta);
    m_omega.push_back(omega);
    m_time.push_back(t);
    
    // Display simulation parameters
    cout << "Driving period T = 2π/ΩD = " << drivingPeriod << " seconds" << endl;
    cout << "Expected Poincare points: ~" << (int)(TIME / drivingPeriod) << endl;
    
    // ===== Main simulation loop =====
    int numSteps = static_cast<int>(TIME / TIME_STEP);
    for (int i = 0; i < numSteps; ++i)
    {
        prevTime = t;
        
        // Advance system by one time step using RK4
        RK4(theta, omega, t);
        t += TIME_STEP;
        
        // Store current state for full trajectory
        m_theta.push_back(theta);
        m_omega.push_back(omega);
        m_time.push_back(t);
        
        // Check if we crossed a sampling time between prevTime and t
        // This catches the sample even if we step over it
        if (prevTime < nextSampleTime && t >= nextSampleTime)
        {
            // We crossed a sampling time - store this point
            m_poincare_theta.push_back(theta);
            m_poincare_omega.push_back(omega);
            
            // Advance to next sampling time
            nextSampleTime += drivingPeriod;
        }
    }
    
    // Display simulation results
    cout << "Simulation complete. " << m_theta.size() << " total data points." << endl;
    cout << "Poincare section contains " << m_poincare_theta.size() << " points." << endl;
}

/**
 * Write trajectory data to files and generate plots using gnuplot
 * 
 * Creates three output files:
 * 1. pendulum_trajectory.out - Full trajectory data (time, theta, omega)
 * 2. poincare_section.out - Poincaré section points (theta, omega)
 * 3. chaotic_pendulum.gnu - Gnuplot script for generating plots
 * 
 * Generates three PNG plots:
 * 1. poincare_section.png - Poincaré section showing strange attractor
 * 2. phase_space_full.png - Complete phase space trajectory
 * 3. time_series.png - Time evolution of theta and omega
 */
void Pendulum::plotTrajectory()
{
    // ===== Write full trajectory data file =====
    ofstream outFile(m_fileName);
    
    if (!outFile.is_open())
    {
        cerr << "Error: Could not open file " << m_fileName << endl;
        return;
    }
    
    // Write header comments
    outFile << "# Driven Damped Pendulum Trajectory (Full Simulation)" << endl;
    outFile << "# Columns: time(s)  theta(rad)  omega(rad/s)" << endl;
    outFile << fixed << setprecision(6);
    
    // Write all trajectory data
    for (size_t i = 0; i < 720/TIME_STEP; ++i)
    {
        outFile << m_time[i] << "\t" << m_theta[i] << "\t" << m_omega[i] << endl;
    }
    outFile.close();

    // ===== Write Poincaré section data file =====
    string poincareFile = "poincare_section.out";
    ofstream poincareOut(poincareFile);
    
    if (!poincareOut.is_open())
    {
        cerr << "Error: Could not open Poincaré file" << endl;
        return;
    }
    
    // Write header with sampling information
    poincareOut << "# Poincaré Section - Points sampled at t = 2πn/ΩD" << endl;
    poincareOut << "# where ΩD = " << m_drivingForceFreq << " rad/s" << endl;
    poincareOut << "# Sampling period = " << (2.0 * M_PI / m_drivingForceFreq) 
                << " seconds" << endl;
    poincareOut << "# Columns: theta(rad)  omega(rad/s)" << endl;
    poincareOut << fixed << setprecision(6);
    
    // Write Poincaré section points
    for (size_t i = 0; i < m_poincare_theta.size(); ++i)
    {
        poincareOut << m_poincare_theta[i] << "\t" 
                    << m_poincare_omega[i] << endl;
    }
    poincareOut.close();
    
    cout << "\nTrajectory data written to " << m_fileName << endl;
    cout << "Poincaré section data written to " << poincareFile << endl;

    // ===== Create gnuplot script =====
    string gnuplotScript = "chaotic_pendulum.gnu";
    ofstream gnuFile(gnuplotScript);

    if (!gnuFile.is_open())
    {
        cerr << "Error: Could not create gnuplot script" << endl;
        return;
    }

    // --- Plot 1: Poincaré Section ---
    // This shows the strange attractor structure by sampling at constant phase
    gnuFile << "set terminal png size 800,600" << endl;
    gnuFile << "set output 'poincare_section.png'" << endl;
    gnuFile << "set title 'Poincaré Section - Driven Damped Pendulum'" << endl;
    gnuFile << "set xlabel 'Theta (radians)'" << endl;
    gnuFile << "set ylabel 'Omega (rad/s)'" << endl;
    gnuFile << "set grid" << endl;
    // Use small points (ps 0.3) to see dense structure clearly
    gnuFile << "plot '" << poincareFile 
            << "' using 1:2 with points pt 7 ps 0.1 lc rgb 'blue' notitle" << endl;
    gnuFile << endl;
    
    // --- Plot 2: Full Phase Space ---
    // This shows the complete trajectory as a continuous line
    gnuFile << "set output 'phase_space_full.png'" << endl;
    gnuFile << "set title 'Full Phase Space - Driven Damped Pendulum'" << endl;
    gnuFile << "plot '" << m_fileName 
            << "' using 2:3 with points pt 7 ps 0.001 lc rgb 'red' title 'Full Trajectory'" << endl;
    gnuFile << endl;
        
    // --- Plot 3: Time Series ---
    // This shows how theta and omega evolve over time
    gnuFile << "set output 'time_series.png'" << endl;
    gnuFile << "set title 'Time Series - Driven Damped Pendulum'" << endl;
    gnuFile << "set xlabel 'Time (s)'" << endl;
    gnuFile << "set ylabel 'Angle/Velocity'" << endl;
    gnuFile << "plot '" << m_fileName 
            << "' using 1:2 with lines lw 1 title 'Theta (rad)', \\" << endl;
    gnuFile << "     '" << m_fileName 
            << "' using 1:3 with lines lw 1 title 'Omega (rad/s)'" << endl;
    
    gnuFile.close();
    
    // ===== Execute gnuplot to generate plots =====
    cout << "\nGenerating plots with gnuplot..." << endl;
    int result = system(("gnuplot " + gnuplotScript).c_str());
    
    if (result == 0)
    {
        cout << "Plots generated successfully:" << endl;
        cout << "  - poincare_section.png (Poincaré section at t = 2πn/ΩD)" << endl;
        cout << "  - phase_space_full.png (complete trajectory)" << endl;
        cout << "  - time_series.png (theta and omega vs time)" << endl;
    }
    else
    {
        cout << "Note: Gnuplot may not be installed or not in PATH." << endl;
        cout << "To plot manually, run: gnuplot " << gnuplotScript << endl;
    }
}

/**
 * Display pendulum properties and simulation parameters
 * Prints a formatted summary of all system parameters to the console
 */
void Pendulum::outputProperties() const
{
    cout << "\n=== Pendulum Properties ===" << endl;
    cout << "Mass: " << m_mass << " kg" << endl;
    cout << "Length: " << m_length << " m" << endl;
    cout << "Damping Coefficient: " << m_dampingCoefficient << endl;
    cout << "Initial Angle: " << m_initialAngle << " rad" << endl;
    cout << "Initial Angular Velocity: " << m_initialOmega << " rad/s" << endl;
    cout << "Driving Force Magnitude: " << m_drivingForceMag << endl;
    cout << "Driving Force Angular Frequency: " << m_drivingForceFreq << " rad/s" << endl;
    cout << "Simulation Time: " << TIME << " s" << endl;
    cout << "Time Step: " << TIME_STEP << " s" << endl;
    cout << "===========================\n" << endl;
}