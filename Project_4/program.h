#ifndef PROGRAM_H
#define PROGRAM_H

#include <vector>
#include <string>

// Physical constants
const double EPSILON_0 = 8.854e-12; // F/m
const double K = 8.99e9; // N*m^2/C^2
const double OMEGA = 1.5; // Over-relaxation parameter (1 < omega < 2)
const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 10000;

// Grid parameters
const double GRID_MIN = -0.5; // -50 cm in meters
const double GRID_MAX = 0.5;  // 50 cm in meters
const double GRID_SPACING = 0.025; // 2.5 cm in meters
const int GRID_SIZE = (int)((GRID_MAX - GRID_MIN) / GRID_SPACING) + 1;

// Disk parameters
const double DISK_RADIUS = 0.10; // 10 cm in meters
const double DISK_VOLTAGE = 0.25; // V

// Point charge parameters
const double CHARGE = 1e-6; // 1 uC in Coulombs

struct PointCharge {
    double x, y, z; // position in meters
    double q; // charge in Coulombs
};

class PotentialSolver {
private:
    std::vector<std::vector<std::vector<double>>> V; // Potential array
    std::vector<std::vector<std::vector<bool>>> fixed; // Fixed potential points
    std::vector<PointCharge> charges;
    
public:
    PotentialSolver();
    void initializeConditions();
    double iterate();
    void solve();
    void outputSlice(const std::string& filename, int plane_index, char axis);
    void outputVTK(const std::string& filename);
};

#endif // PROGRAM_H