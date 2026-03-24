# Project 9 main file: finding the energy levels of quantum anharmonic
# oscillators. Note that my versions of the code write all of the output
# to stdout. In a unix environment, stdout can be piped to a file as
#
# $ python3 main.py > output.txt
#
# A note on units: we are assuming hbar=1 and m=1. This results in length
# units that are really funky and, more importantly, energy units that are
# essentially in units of hbar*omega (convenient for comparing n=2 levels
# to the analytic solutions for the harmonic oscillator).
#
# Author: Kevin Kelley (kelleyk@byui.edu)
# Date: 3/13/25

from potential import PowerPotential

# These are the parameters I'm using in my energy level searches:

nmax=10   # The maximum value of n to consider (even values only!)
lmax=10   # The maximum bound state energy to find (ground state is l=0)
dx=0.005  # The spatial discretization parameter (length units)
dE=0.02   # The initial step size to use when searching (energy units)
tol=1.e-6 # Quit searching when the energy step size gets this small
bdiv=10.0 # A difference in psi of this much between two grid points indicates
          # that the solution is diverging.

# First find the energies using the shooting method.

for n in range(2,nmax+1,2):           # Loop over even values of n
    v=PowerPotential(n,dx)            # Create the potential
    v.find_shooting(lmax,dE,tol,bdiv) # Find the energies, printing to stdout

# Now do it again using the matching method.

for n in range(2,nmax+1,2):           # Loop over even values of n
    v=PowerPotential(n,dx)            # Create the potential
    v.find_matching(lmax,dE,tol)      # Find the energies, printing to stdout
