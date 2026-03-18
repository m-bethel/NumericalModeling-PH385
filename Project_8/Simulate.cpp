/*
Simulate.cpp

Implements the Simulate class. Manages the simulation loop across all
particles, exports results to a CSV file, and generates a gnuplot script 
to visualize the output.

Output is a png that represents Tempuerature vs. Total Internal Energy 
E = KE + PE

(README)

Adjusting Parameters:

Most of what you "should" adjust is defined as a member variable, so stick to that.

If you adjust the parameters to the point where it's not very stable, you can always add more particles
The more particles, the more stable and nicer it is. You kinda just trade computational time for stability
Feel free to use whatever number, the program will place them at the center in a grid just fine.

m_coolingSteps: Keep in mind that the m_coolingSteps should be ~ 0.1 * m_steps for best stability

m_sampleRate: m_sampleRate is also something that can be adjusted. It provides the average over 
that many steps in hopes to smooth out the massive spikes in the data.
400 is the very low end. Higher will smooth at the data, however, it will lessen 
the total datapoints accordingly. m_steps/m_sampleRate = data points.
This will effect the final temp you end at because the deltaEnergy will only add energy 
every m_sampleRate as well.
Example: If I want more detail, and for it be be smoother, I will increase 
m_steps(300000) -> m_steps(400000). I need to also increase m_coolingSteps(30000) -> m_coolingSteps(40000)
Because I increase the steps by 1/3, I could also do the same to m_sampleRate(400) -> m_sampleRate(533)
OR
navigate to the main.cpp file and reduce deltaEnergy = 0.0055 -> 0.0037
OR
a mixture of both

If curiosity gets the better of you, you can adjust a few others:

Q_thermostat(default = 4.0) the higher the value, the less aggresive the temp change is

m_xi clamps 
(default)
    if (m_xi > 0.01) m_xi = 0.01;
    if (m_xi < -0.01) m_xi = -0.01;
these prevent how hard the initial jump in velocity should be. Keep he values betwen 0.05 - 0.005. 
Outside of that doesn't work

p.capVelocity(2.0): For the absolute velocity cap, I would not go above 4.0 at the most. 
Keeping it to 3.0 and less is best for the gif visulaization. Above 4 makes it unstable.

FINAL PIECE OF ADVICE / HAIL MARY

If you can't restore the code for whatever reason, you can find it on my github repo @:

https://github.com/m-bethel/NumericalModeling-PH385/Project_8

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/1/2026
*/

#include "Simulate.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

Simulate::Simulate(int particles, double deltaEnergy, int gridSize)
    : m_deltaEnergy(deltaEnergy),
      m_gridSize(gridSize),
      m_steps(300000),
      m_xi(0.0),
      m_targetTemp(0.2),
      m_heatingUp(true),
      m_totalEnergyAdded(0.0),
      m_coolingSteps(30000),
      m_sampleRate(400),
      m_targetQuench(1E-3)
{
    double spacing = 1.122;
    int cols = ceil(sqrt(particles));
    double offset = (m_gridSize - (cols * spacing)) / 2.0;
    
    for (int i = 0; i < particles; i++)
    {
        double x = (i % cols) * spacing + offset;
        double y = (i / cols) * spacing + offset;
        m_particles.push_back(Particle(gridSize, x, y));
    }
}

void Simulate::simulate(double deltaTime)
{
    ofstream animFile("animation.dat");
    const double Q_thermostat = 4.0; // Higher value (5-10) makes simulation smoother
    
    // Initial Quenching setup
    double deltaQuench = (m_targetTemp - m_targetQuench) / m_coolingSteps;
    
    // Accumulators for averaging data
    double tempSum = 0.0;
    double energySum = 0.0;
    int sampleCount = 0;

    for (int j = 0; j < m_steps; j++)
    {
        // 1. Position Update (Velocity Verlet Step 1)
        for (Particle& p : m_particles) {
            p.savePrevForce();      // Store current force as 'prev' for velocity update
            p.updatePosition(deltaTime);
            p.periodic();           // Periodic boundary conditions
        }

        // 2. Force Calculation (Velocity Verlet Step 2)
        // Returns the new Potential Energy (PE) at the updated positions
        double currentPE = calculateForces();

        // 3. Velocity Update (Velocity Verlet Step 3)
        for (Particle& p : m_particles) {
            p.updateVelocity(deltaTime);
            // Cap speed slightly during quenching to prevent initial "explosions"
            if (j < m_coolingSteps) p.capVelocity(2.0);
        }

        // 4. Thermodynamics Calculations
        double totalKE = 0.0;
        for (Particle& p : m_particles) totalKE += p.kineticEnergy();
        double currentT = totalKE / m_particles.size();
        double currentE = (totalKE + currentPE) / m_particles.size();

        // 5. Thermostat Logic (Berendsen)
        m_xi += (currentT / m_targetTemp - 1.0) * (deltaTime / Q_thermostat);
        
        // Clamp xi to prevent the system from getting "kicked" too hard
        if (m_xi > 0.01) m_xi = 0.01;
        if (m_xi < -0.01) m_xi = -0.01;

        double scaleFactor = sqrt(fabs(1.0 - m_xi));
        for (Particle& p : m_particles) p.scaleVelocity(scaleFactor);

        // 6. Heating / Cooling Schedule
        if (j < m_coolingSteps) {
            // Initial Quench phase: cool down to the starting target
            if (m_targetTemp > m_targetQuench) m_targetTemp -= deltaQuench;
        } 
        else {
            // Main Simulation: Adjust target temperature every sample window
            if ((j - m_coolingSteps) % m_sampleRate == 0) {
                // First half of simulation is Heating, second half is Cooling
                bool isHeating = (j < (m_steps + m_coolingSteps) / 2);
                m_targetTemp += isHeating ? m_deltaEnergy : -m_deltaEnergy;
            }
        }

        // 7. Data Collection for Graphing
        // We only sample after the quenching phase is finished
        if (j >= m_coolingSteps) {
            tempSum += currentT;
            energySum += currentE;
            sampleCount++;

            // At the end of the sample window, record averages to vectors
            if ((j - m_coolingSteps) % m_sampleRate == m_sampleRate - 1) {
                m_temp.push_back(tempSum / sampleCount);      // Y-axis: Temperature
                m_potEnergy.push_back(energySum / sampleCount); // X-axis: Total Energy
                
                tempSum = 0.0;
                energySum = 0.0;
                sampleCount = 0;
            }
        }

        // 8. Animation Frame Output
        if (j % m_sampleRate == 0) {
            for (size_t i = 0; i < m_particles.size(); i++) {
                const double* pos = m_particles[i].getPosition();
                animFile << pos[0] << " " << pos[1] << "\n";
            }
            animFile << "\n\n"; // Required for Gnuplot 'index' command
        }
    }
    animFile.close();
}

double Simulate::calculateForces()
{
    double totalPE = 0.0; // Track potential energy
    for (Particle& p : m_particles) p.setForce(0.0, 0.0);

    const double cutoffSq = 6.25;
    double reducedGridSize = m_gridSize / 2.0;

    for (size_t i = 0; i < m_particles.size(); i++) {
        for (size_t j = i + 1; j < m_particles.size(); j++) {

            const double* posI = m_particles[i].getPosition();
            const double* posJ = m_particles[j].getPosition();

            double dx = posI[0] - posJ[0];
            double dy = posI[1] - posJ[1];

            if (dx >  reducedGridSize) dx -= m_gridSize;
            if (dx < -reducedGridSize) dx += m_gridSize;
            if (dy >  reducedGridSize) dy -= m_gridSize;
            if (dy < -reducedGridSize) dy += m_gridSize;

            double r2 = dx*dx + dy*dy;

            if (r2 < cutoffSq) {
                //if (r2 < 0.1) r2 = 0.1; 
                double invR2 = 1.0 / r2;
                double invR6 = invR2 * invR2 * invR2;
                double invR12 = invR6 * invR6;

                double fOverR = 24.0 * (2.0 * invR12 - invR6) * invR2;
                double fx = fOverR * dx;
                double fy = fOverR * dy;

                m_particles[i].addForce(fx, fy);
                m_particles[j].addForce(-fx, -fy);

                // Potential Energy calculation: V = 4 * ( (1/r^12) - (1/r^6) )
                totalPE += 4.0 * (invR12 - invR6);
            }
        }
    }
    return totalPE; // Return the sum
}

void Simulate::output()
{
    int halfway = (int)m_temp.size()/ 2;

    // 1. Write Data
    ofstream heatFile("heating.csv");
    ofstream coolFile("cooling.csv");
    for (int i = 0; i < halfway; i++) {
        heatFile << m_potEnergy[i]  << "," << m_temp[i] << endl;
        }
    heatFile.close();

    for (int i = halfway; i < m_temp.size(); i++) {
        coolFile << m_potEnergy[i]  << "," << m_temp[i] << endl;
        }
    coolFile.close();

    // 3. Update Gnuplot Script
    ofstream et("et_plot.gnu");
    et << "set term qt 1 title 'Temperature vs E' noraise\n";
    et << "set datafile separator ','\n";
    et << "set terminal png size 1200,600\n";
    et << "set output 'phase_transition.png'\n";
    et << "set xlabel 'Total Internal Energy (E)'\n";
    et << "set ylabel 'Temperature (k_B T)'\n";
    et << "set grid\n";
    et << "set key left top\n";
    et << "set offsets graph 0.05, 0.05, 0.05, 0.05\n";
    et << "set title 'Phase Transition'\n";
    et << "plot 'heating.csv' using 1:2 with lines lc rgb 'red' lw 1 title 'Heating', \\\n";
    et << "     'cooling.csv' using 1:2 with lines lc rgb 'blue' lw 1 title 'Cooling'\n";
    et.close();

    ofstream animScript("animation.gnu");
    animScript << "set terminal gif animate delay 15 size 900,900\n";
    animScript << "set output 'animation.gif'\n";
    animScript << "set xrange [0:" << m_gridSize << "]\n";
    animScript << "set yrange [0:" << m_gridSize << "]\n";
    animScript << "set size square\n";
    animScript << "set style fill solid 1.0\n"; 
    animScript << "unset key\n"; // Removes the filename legend
        
    int totalFrames = m_steps / m_sampleRate;
    animScript << "do for [i=0:" << totalFrames - 1 << "] {\n";
    // 'with points pt 7' = solid circles
    // 'ps 1.5' = size of the circles
    animScript << "    plot 'animation.dat' index i with points pt 7 ps 1.5 lc rgb 'blue'\n";
    animScript << "}\n";
    animScript.close();
    // Final execution sequence
    
    if (system("gnuplot animation.gnu ") != 0) {
        std::cerr << "Error: Gnuplot failed to execute." << std::endl;
        }
    if (system("gnuplot -p et_plot.gnu") != 0) {
        std::cerr << "Error: Gnuplot failed to execute." << std::endl;
        }

}