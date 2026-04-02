#ifndef PROGRAM_H
#define PROGRAM_H

#include <vector>
#include <string>

// Physical constants for electrostatics
const double EPSILON_0 = 8.854e-12; // Vacuum permittivity (Farads/meter)
const double K = 8.99e9;            // Coulomb's constant (N*m^2/C^2), used for initial V = kq/r guesses
const double TOLERANCE = 1e-6;      // The convergence threshold; iteration stops when max change is below this
const int MAX_ITERATIONS = 10000;   // Failsafe to prevent infinite loops if the solver diverges (e.g., omega >= 2.0)

// Grid parameters defining the 3D physical space
const double GRID_MIN = -0.5;      // The start of the bounding box at -50 cm (in meters)
const double GRID_MAX = 0.5;       // The end of the bounding box at +50 cm (in meters)
const double GRID_SPACING = 0.025; // The physical distance between adjacent nodes (resolution)
// Calculates the total number of nodes needed along one axis to cover the space
const int GRID_SIZE = (int)((GRID_MAX - GRID_MIN) / GRID_SPACING) + 1;

// Disk parameters
const double DISK_RADIUS = 0.10;  // The radius of the central conducting disk (10 cm)
const double DISK_VOLTAGE = 0.25; // The fixed Dirichlet boundary condition for the disk (0.25 Volts)

// Point charge parameters
const double CHARGE = 1e-6; // The magnitude of the point charges (1 microCoulomb)

// A simple structure to hold the 3D coordinates and magnitude of a point charge
struct PointCharge {
    double x, y, z; // Physical position in meters
    double q;       // Charge magnitude in Coulombs
};

// The main class that handles the 3D Poisson equation solving
class PotentialSolver {
private:
    double omega; // The Successive Over-Relaxation (SOR) parameter
    // A 3D vector array storing the electric potential (Volts) at every (i, j, k) index
    std::vector<std::vector<std::vector<double>>> V; 
    // A 3D vector array acting as a mask; true = boundary condition (locked), false = free to update
    std::vector<std::vector<std::vector<bool>>> fixed; 
    std::vector<PointCharge> charges; // A list holding the 4 point charges
    
public:
    PotentialSolver(); // Constructor to set up initial memory and charge positions
    void setOmega(double w){ omega = w;} // Setter function to allow dynamic testing of different omega values

    void initializeConditions(); // Function to apply boundary conditions and initial guesses
    double iterate();            // Performs a single full sweep of the 3D grid using the SOR method
    int solve();                 // Loops the iterate() function until convergence is reached
    // Function to export a 2D cross-section of the 3D grid to a text file for plotting
    void outputSlice(const std::string& filename, int plane_index, char axis);
};

#endif // PROGRAM_H