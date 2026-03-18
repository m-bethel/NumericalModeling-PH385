/*
Solve.h

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/16/2026
*/

#ifndef SOLVE_H
#define SOLVE_H

#include <string>
#include <vector>


class Solve
{
    private:
        int m_nmax;
        int m_lmax;
        double m_dx;
        double m_dE;
        double m_tol;
        double m_bdiv;
    
    public:
        Solve(int nmax, int lmax, double dx, double dE, 
                   double tol, double bdiv);
        void output();

};

#endif //SOLVE_H