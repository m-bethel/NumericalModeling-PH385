/*
================================================================================
FINAL PROJECT: 2D COMPUTATIONAL FLUID DYNAMICS (CFD) ENGINE
================================================================================
Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/23/2026

README: A CRASH COURSE IN FLUID SIMULATION FOR NEW PHYSICS STUDENTS
--------------------------------------------------------------------------------
Welcome to Computational Fluid Dynamics (CFD)! If you are new to physics programming, 
simulating water or air might seem like magic. It is actually just the rigorous 
application of two fundamental rules of the universe:

1. Conservation of Mass: Fluid cannot be created or destroyed. What goes in, must come out.
2. Conservation of Momentum (Newton's 2nd Law, F=ma): Fluid moves when pushed by forces.

In fluids, these rules are described by the "Navier-Stokes Equations". Because these 
equations are too complex to solve perfectly with algebra, we use computers to 
approximate them using a technique called "The Projection Method" (invented by Alexandre Chorin).

HOW THIS ALGORITHM WORKS (The 3-Step Dance):
Imagine the fluid domain as a giant grid of tiny squares (like a spreadsheet). 
Every time step, the computer does the following:

  STEP 1: THE PREDICTOR (Stokes.cpp)
  We calculate how the fluid *wants* to move based on its current momentum and viscosity 
  (thickness). However, this prediction is "flawed" because it ignores pressure. The 
  fluid might accidentally squish together or overlap, violating the Conservation of Mass.

  STEP 2: THE PRESSURE SOLVER (Poisson.cpp)
  We look at our flawed prediction and ask: "Where is the fluid bunching up too much, 
  and where is it spreading out?" We use a mathematical tool called a "Poisson Equation" 
  to calculate a Pressure Field. High pressure forms where fluid is bunching up; low 
  pressure forms where it is spreading out.

  STEP 3: THE CORRECTOR (Stokes.cpp)
  We use the pressure field from Step 2 to push back on the fluid. Fluid flows from 
  high pressure to low pressure. This final "correction" fixes the flawed prediction, 
  ensuring the fluid is perfectly incompressible.

THE SCENARIO:
This code simulates fluid flowing through a wind tunnel (left to right) past a 
circular obstacle. At the right speeds, this creates a beautiful, alternating 
chain of whirlpools called a "von Kármán vortex street."

================================================================================
*/

#include "Poisson.h"  // Contains the math to solve for the pressure field
#include "Mesh.h"     // Contains the grid definitions (the physical space)
#include "Stokes.h"   // Contains the math to move the fluid (velocity)
#include <iostream>   // For printing text to the console
#include <limits>     // For checking mathematical limits
#include <string>     // For handling file names

using namespace std;

int main()
{
    // --- 1. SIMULATION TIME PARAMETERS ---
    float dt = 0.001f;     // "Delta Time": How many seconds pass per simulation step
    float t = 0.0f;        // Current simulation time
    float tEnd = 25.0f;    // Total seconds to simulate (Long enough for vortices to form)
    int frameCount = 0;    // Keeps track of how many 3D video frames we've saved
    int step = 0;          // Keeps track of the total number of math loops calculated

    // --- 2. SOLVER PARAMETERS ---
    int max_Iterations = 500; // Max attempts the pressure solver gets to find an answer per step
    float tol = 1E-3f;        // "Tolerance": How accurate the pressure solver needs to be (Error margin)

    // --- 3. FLUID PROPERTIES ---
    float rho = 1.0f;      // Density of the fluid (1.0 = water standard)
    float nu = 0.001f;     // Kinematic Viscosity: How "sticky" or thick the fluid is
    float U_in = 1.3f;     // The speed of the fluid entering the left side of the tunnel

    // --- 4. MESH / GRID PARAMETERS ---
    int nx = 1024;         // Number of cells in the X direction (Horizontal resolution)
    int ny = 256;          // Number of cells in the Y direction (Vertical resolution)
    float Lx = 4.0f;       // Physical length of the tunnel (e.g., 4 meters)
    float Ly = 1.0f;       // Physical height of the tunnel (e.g., 1 meter)

    // --- 5. STABILITY CHECKS (CFL CONDITION) ---
    // In CFD, if fluid moves across more than one grid cell in a single time step, 
    // the math explodes. These checks ensure our 'dt' is small enough to be safe.
    float dx = Lx / nx;                     // Physical width of one grid cell
    float cfl_conv = U_in * dt / dx;        // Convective CFL number (Must be < 1.0)
    float cfl_diff = nu * dt / (dx * dx);   // Diffusive CFL number (Must be < 0.5)

    // --- 6. INITIALIZE PHYSICS OBJECTS ---
    Mesh mesh2D(nx, ny, Lx, Ly);            // Create the physical space grid
    Poisson solver(mesh2D);                 // Create the pressure math engine
    Stokes stokes(mesh2D, nu);              // Create the velocity math engine

    // --- 7. MAIN SIMULATION LOOP ---
    // This loop runs continuously until the simulated time reaches 'tEnd'
    while (t < tEnd) {

        // Step A: Set the walls and inlet speeds for the current velocities
        stokes.applyBoundary(U_in);
        
        // Step B: PREDICTOR - Calculate how momentum and viscosity move the fluid
        stokes.predict(dt, mesh2D.getDxi(), mesh2D.getDyi()); 
        
        // Step C: Set the walls and inlet speeds for the *predicted* velocities
        stokes.applyBoundaryToStar(U_in);
        
        // Step D: Calculate where the fluid is violating conservation of mass (bunching up)
        solver.buildRHS(stokes.getUStar(), stokes.getVStar(), dt, rho, mesh2D.getDxi(), mesh2D.getDyi());

        // Step E: PRESSURE SOLVER - Calculate the pressure needed to fix the mass violation
        // (We only solve every 1 step, but you could optimize to solve less frequently)
        if (step % 1 == 0){
            solver.solve(max_Iterations, tol, mesh2D.getDxi(), mesh2D.getDyi());
        }
        
        // Step F: CORRECTOR - Use the new pressure to fix the predicted velocities
        stokes.correct(solver.getP(), dt, rho, mesh2D.getDxi(), mesh2D.getDyi());
        
        // Step G: Re-apply boundaries to ensure the corrected fluid isn't flowing through walls
        stokes.applyBoundary(U_in);
        
        // Step H: EXPORT DATA TO HARD DRIVE
        // We don't save every step, otherwise our hard drive would fill up. 
        // We save a snapshot every 50 steps for the Python script to animate later.
        if (step % 50 == 0){
            string filename = "frame_" + to_string(frameCount) + ".bin";
            stokes.exportFrame(filename, frameCount, solver.getP());
            frameCount++; // Advance the frame counter for the next video frame
        }
        
        t += dt;  // Move time forward by one delta-time step
        step++;   // Increment the loop counter
    }
    
    return 0;     // End the program cleanly
}