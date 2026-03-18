/*
Simulate.h

Declares the Sumlate class which handles the full simulation,
Owns a collections of Particle Object, runs the simulation loop,
and handles all file output including the graph and animation

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/2/2026
*/

#ifndef SIMULATE_H
#define SIMULATE_H

#include <string>       // std::string for file name storage
#include <vector>       // std::vector for particle collection and MSD storage
#include "Particle.h"   // Particle class used in the simulation

class Simulate
{
    private:
        std::vector<Particle> m_particles;
        std::string m_filename;
        std::vector<double> m_temp;
        std::vector<double> m_specificHeat;
        std::vector<double> m_potEnergy;
        int m_gridSize;
        int m_steps;
        int m_sampleRate;
        int m_coolingSteps;
        double m_targetQuench;
        double m_deltaEnergy;
        double m_xi;
        double m_targetTemp;
        bool m_heatingUp;
        double m_totalEnergyAdded;

    public:
        Simulate(int particles, double deltaEnergy, int gridSize);
        double calculateForces();
        
        void simulate (double deltaTime);
        void output();
};

#endif //SIMULATE_H