/*
main.cpp

Final Project of 385 Numerical Modeling
Attempting a Navier Stokes solver for modeling flow around an object
in 2D first.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/23/2026

*/

#include "Poisson.h"
#include "Mesh.h"
#include "Stokes.h"
#include <iostream>
#include <limits>

using namespace std;

int main()
{
    // Simulation Parameters
    float dt = 0.001;
    float t = 0;
    float tEnd = 35;
    float U_in = 1.3;
    int frameCount = 0;
    int step = 0;

    int max_Iterations = 500;
    float tol = 1E-3;

    // Fluid Properties
    float rho = 1.0;
    float nu = 0.001;

    //MESH GRID PARAMETERS
    int nx = 1024;
    int ny = 256;
    float Lx = 4.0;
    float Ly = 1.0;

    //Chec to make sure values are within spec
    float dx = Lx / nx;
    float cfl_conv = U_in * dt / dx;
    float cfl_diff = nu * dt / (dx * dx);

    Mesh mesh2D(nx, ny, Lx, Ly); 
    Poisson solver(mesh2D);
    Stokes stokes (mesh2D, nu);

    while (t < tEnd){

        stokes.applyBoundary(U_in);
        stokes.predict(dt, mesh2D.getDxi(), mesh2D.getDyi()); 
        stokes.applyBoundaryToStar(U_in);
        solver.buildRHS(stokes.getUStar(), stokes.getVStar() , dt, rho,
                        mesh2D.getDxi(), mesh2D.getDyi());

        if (step % 1 == 0){
            solver.solve(max_Iterations, tol, mesh2D.getDxi(), mesh2D.getDyi());
        }
        stokes.correct(solver.getP(), dt, rho, mesh2D.getDxi(), mesh2D.getDyi());
        stokes.applyBoundary(U_in);
        if (step % 50 == 0){
            string filename = "frame_" + to_string(frameCount) + ".bin";
            stokes.exportFrame(filename, frameCount, solver.getP());
            
            frameCount++;
        }
        t += dt;
        step++;
    }
}
