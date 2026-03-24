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
    double tEnd = 100;

    // Fluid Properties
    double rho = 1.27;
    double nu = 0.1;

    //MESH GRID PARAMETERS
    int nx = 64;
    int ny = 32;
    int jmin = 2;
    int imin = 2;
    int imax = imin + nx -1;
    int jmax = jmin + ny -1;
    double Lx = 10.0;
    double Ly = 1.0;
    double dx, dy, dxi, dyi;
    //std::vector<double> x,y,xm,ym;


    Mesh mesh2D(nx, ny, Lx, Ly, imin, imax, jmin, jmax,
                 dx, dy, dxi, dyi);
    mesh2D.setUpMesh();
    
    Poisson solver(mesh2D);

    Stokes stokes ();

    while (t < tEnd){
        stokes.predict(dt);

        solver.buildRHS(stokes.getUstart(), stokes.getVStart(), dt, rho, dxi, dyi);

        solver.solve();

        stokes.correct(solver.getP(), dt, rho, dxi, dyi);

        t += dt;

    }


}
