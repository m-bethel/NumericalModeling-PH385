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
    double dt = 0.001;
    double t = 0;
    double tEnd = 10;
    double U_lid = 1.0;
    int frameCount = 0;
    int step = 0;

    int max_Iterations = 2000;
    double tol = 1E-6;

    // Fluid Properties
    double rho = 1.0;
    double nu = 0.1;

    //MESH GRID PARAMETERS
    int nx = 32;
    int ny = 16;
    double Lx = 1.0;
    double Ly = 1.0;

    //Chec to make sure values are within spec
    double dx = Lx / nx;
    double cfl_conv = U_lid * dt / dx;
    double cfl_diff = nu * dt / (dx * dx);
    cout << "CFL convective: " << cfl_conv << endl;
    cout << "CFL diffusive:  " << cfl_diff << endl;
    if (cfl_conv > 0.5 || cfl_diff > 0.5)
        cout << "WARNING: CFL condition violated!" << endl;


    Mesh mesh2D(nx, ny, Lx, Ly); 
    mesh2D.setUpMesh();
    
    Poisson solver(mesh2D);

    Stokes stokes (mesh2D, nu);

    while (t < tEnd){

        stokes.applyBoundary(U_lid);
    
        stokes.predict(dt, mesh2D.getDxi(), mesh2D.getDyi());
        
        stokes.applyBoundaryToStar(U_lid);

        solver.buildRHS(stokes.getUStar(), stokes.getVStar() , dt, rho,
                        mesh2D.getDxi(), mesh2D.getDyi());

        solver.solve(max_Iterations, tol);

        stokes.correct(solver.getP(), dt, rho, mesh2D.getDxi(), mesh2D.getDyi());
        
        stokes.applyBoundary(U_lid);

        if (step % 100 == 0){
            string filename = "frame_" + to_string(frameCount) + ".csv";
            stokes.exportFrame(filename, frameCount, solver.getP());
            frameCount++;
        }

        t += dt;
        step++;
    }


}
