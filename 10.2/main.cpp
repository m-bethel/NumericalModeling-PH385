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
    int steps = 1E6;
    double t = 0;
    double dt = 0.001;
    
    double U_lid = 1.0;
    int frameCount = 0;
    int step2 = 0;

    int max_Iterations = 2000;
    double tol = 1E-6;

    // Fluid Properties
    double rho = 1.0;
    double nu = 0.1;

    //MESH GRID PARAMETERS
    int nx = 10;
    int ny = 5;
    int imin = 1;
    int jmin = imin;
    int imax = imin + nx -1;
    int jmax = jmin + ny -1;
    double Lx = 0.1;
    double Ly = 0.05;

    Mesh mesh2D(imin, jmin, imax, jmax, Lx, Ly, nx, ny); 
    mesh2D.setUpMesh();
    
 
    Stokes stokes (nx, ny, imin, jmin, imax, jmax, mesh2D, nu);

   //Poisson solver(mesh2D);
//
    //while (t < tEnd){
//
    //    stokes.applyBoundary(U_lid);
    //
    //    stokes.predict(dt, mesh2D.getDxi(), mesh2D.getDyi());
    //    
    //    stokes.applyBoundaryToStar(U_lid);
//
    //    solver.buildRHS(stokes.getUStar(), stokes.getVStar() , dt, rho,
    //                    mesh2D.getDxi(), mesh2D.getDyi());
//
    //    solver.solve(max_Iterations, tol);
//
    //    stokes.correct(solver.getP(), dt, rho, mesh2D.getDxi(), mesh2D.getDyi());
    //    
    //    stokes.applyBoundary(U_lid);
//
    //    if (step % 100 == 0){
    //        string filename = "frame_" + to_string(frameCount) + ".csv";
    //        stokes.exportFrame(filename, frameCount, solver.getP());
    //        frameCount++;
    //    }
//
    //    t += dt;
    //    step++;
    //}


}
