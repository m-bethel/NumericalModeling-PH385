/*

Creating a 3-dimensional plot of charges and a disk held to 0.25V in a grid using 
the Jacobi method.  Iterates through multiple omega values to find the best one.

Author: Miles Bethel    (miles.d.bethel@gmail.com)
Date:03/30/2026
*/


#include <fstream>
#include <chrono> // Included for potential performance timing, though not actively used in this iteration
#include "program.h"
#include <iostream>

using namespace std;

int main() {
    // A predefined list of over-relaxation parameters to test. Approaching 2.0 gets faster, but 2.0 fails.
    vector<double> test_omegas = {1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.95, 1.98, 1.99, 2.0};
    
    // Open a CSV file to log the results of the parameter sweep
    ofstream results("convergence_data.csv");
    results << "omega,iterations" << endl; // Write the CSV header row

    // Loop through each omega value in the vector
    for (double w : test_omegas) {
        cout << "Testing Omega: " << w << "..." << flush; // flush forces the terminal to print immediately
        
        // Instantiate a fresh solver object so memory is completely reset for each test
        PotentialSolver solver;
        solver.setOmega(w); // Inject the current test omega
        solver.initializeConditions(); // Set up the disk and charges
        
        // Run the solver and capture how many iterations it took to reach tolerance
        int iters = solver.solve(); 
        
        // Report success to terminal and write the data point to the CSV file
        cout << " Done in " << iters << " iterations." << endl;
        results << w << "," << iters << endl;
        
        // Generate a dynamic filename based on the omega used (e.g., potential_w_1.900000.dat)
        string filename = "potential_w_" + to_string(w) + ".dat";
        
        // Export the z=0 plane (index 20 in a 41x41x41 grid) to that file for visualization
        solver.outputSlice(filename, 20, 'z'); 
    }
    
    results.close(); // Safely close the CSV file
    return 0; // Terminate execution successfully
}