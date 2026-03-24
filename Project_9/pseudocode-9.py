# Class to implement a quantum mechanical potential of the form V(x)=(x^n)/n
# and solve for the bound state energies.
#
# A note on units: we are assuming hbar=1 and m=1. This results in length
# units that are really funky (basically units of sqrt(1/m)), and more
# importantly, energy units that are essentially in units of hbar*omega 
# (convenient for comparing n=2 levels to the analytic solutions for the 
# harmonic oscillator).
#
# These routines write output to stdout and stderr. In a unix environment, 
# these can be piped to files in the usual way (see main.py for more info).

class PowerPotential:

    # -------------------------------------------------------------------------
    # Constructor: copy the values of n and the spatial discretization dx into
    # the object, and initialize the arrays for positions and wavefunction
    # amplitudes (called y so I don't have to keep typing psi). The constructor
    # also initializes a variable yneg, which is amplitude at x=-dx, and is
    # needed in the shooting method.
    # -------------------------------------------------------------------------

    def __init__(self,n,dx):

    # -------------------------------------------------------------------------
    # Method v: returns the value of the potential for the given position x.
    # -------------------------------------------------------------------------

    def v(self,x):

    # -------------------------------------------------------------------------
    # Method normalize: normalizes the wavefunction stored in the y array so
    # that the integral of psi^2 equals one. This is a necessary thing to do if
    # we want to compare the wavefunctions from the two methods (shooting and
    # matching), as the two methods will produce results at potentially vastly
    # different scales.
    # 
    # The normalization process is pretty straightforward. We just do a 
    # numerical integration of psi^2 (add together all the y^2*dx), then
    # divide each amplitude y by the square root of that integral.
    # -------------------------------------------------------------------------

    def normalize(self):

        # Do the integration

        # Rescale all of the amplitudes

    # -------------------------------------------------------------------------
    # Method initialize_shooting: implement the appropriate boundary conditions
    # at x=0 for the shooting method. The boundary conditions used depend on
    # the parity of the solution.
    # -------------------------------------------------------------------------

    def initialize_shooting(self,parity):

        # Reinitialize the x and y list to length zero.

        # If the parity is even (indicated by parity=0), append (x,y)=(0,1) to
        # their respective lists. Also set yneg equal to 1 (indicating that
        # the slope at x=0 is zero).

        # If the parity is odd (indicated by parity=1), append (x,y)=(0,0) to
        # their respective lists. Also set yneg equal to -dx (which makes the 
        # slope at x=0 equal to 1).

    # -------------------------------------------------------------------------
    # Method solve_shooting: use the discretized Schrodinger equation to find
    # the amplitudes at positions other than x=0. To do this, we need to know
    # the current energy (e) and also the amplitude where we assume that 
    # the solution is diverging (b). I've also added an additional parameter
    # (trim) indicating whether we should trim the diverging tail off of the
    # solution. This is useful for cleaning up the resulting wavefunction
    # before writing it to stdout.
    # -------------------------------------------------------------------------

    def solve_shooting(self,e,b,trim=False):

        # Keep going until the absolute value of the most recent amplitude
        # exceeds b.

            # If we are at x=0, we need to use yneg as the "previous"
            # amplitude. Otherwise we use y[-2].

            # Calculate the next amplitude using the discretized Schrodinger
            # equation

            # Append the new amplitude (along with it's position) to the y
            # and x arrays.

        # Define a new variable div that holds the sign (+ or -) of the 
        # last amplitude in the y array. We're going to return this value.

        # Trim the array down to the closest value to zero from the end.

        # Return the sign of the last amplitude in the array (which tells us
        # if the solution was diverging in the positive or the negative
        # direction.

    # -------------------------------------------------------------------------
    # Method find_shooting: find the bound state energies up to level lmax,
    # noting that l=0 is the ground state. The additional parameters are the
    # initial energy step used in our search (deinit), the energy tolerance at 
    # which we terminate our search (tol), and the amplitude value at which we 
    # assume the solution is diverging
    # -------------------------------------------------------------------------

    def find_shooting(self,lmax,deinit,tol,b):

        # Set the energy equal to deinit. This is a good place to start our
        # search for the ground state energy (which must be positive for these
        # potentials).

        # Loop over the bound states, from 0 up to and including lmax

            # Set the energy step de equal to the proscribed initial value
            # deinit

            # Determine the parity of the solution. I use parity=0 to indicate
            # even parity and parity=1 to indicate odd parity.

            # Set the divergence sign (div) and previous divergence sign
            # (last_div) to zero.

            # For l>0, our energy is going to initially be at the energy of the
            # next lowest state. We will keep stepping higher in energy until
            # we know we are above the energy of the next state. This is
            # indicated by the divergence sign changing from positive to
            # negative, or vice-versa.

                # Increase the energy by 2*de (faster search)

                # set the last divergence sign equal to the current one.

                # Call the shooting initialization method to set our initial 
                # conditions by calling initialize_shooting

                # Solve the Schrodinger equation numerically by calling
                # solve_shooting. Save the returned divergence sign in div

            # Now that we are above the next energy level, start refining the
            # search. Keep going until abs(de) is greater than the tolerance.

                # Call the shooting initialization method to set our initial 
                # conditions by calling initialize_shooting

                # Solve the Schrodinger equation numerically by calling
                # solve_shooting. Save the returned divergence sign in div

                # If the divergence sign is different from the previous value,
                # then we have "crossed over" the "correct" energy. We need to
                # reduce the step size by a factor of two, and then also 
                # multiply by -1 to reverse direction.

                # Add the step de to the energy, and copy the current
                # divergence sign div into last_div.

            # We should now have the (approximate) bound state energy!
            # Run the solver one last time with the energy e-de (the last
            # iteration overshot the solution), trimming the diverging tail.

            # Normalize the solution by calling the normalize method, and 
            # then print the results to stdout. I'm formatting my printing
            # so that the wavefunctions can be easing plotted using gnuplot, 
            # and the energies extracted using the grep command.
            self.normalize()
            print(f"# Shooting: n={self.n:2d}, state={l:2d}, E={e:10.3f}")
            for i in range(0,len(self.x)):
                print(self.x[i],self.y[i])
            print("\n\n")

    # -------------------------------------------------------------------------
    # Method initialize_matching: initialize the x and y lists. Since we are
    # starting at large x and solving inward, we will need to initialize the 
    # entire array. I'm assuming that at three times the half width of the
    # potential at a given energy (passed in as the parameter e) the 
    # amplitude of the wavefunction will essentially be zero.
    # -------------------------------------------------------------------------

    def initialize_matching(self,e):

        # Determine the (half) width of the well at the given energy

        # Extend the arrays to 3x the width, saving the position values in the
        # x list, and setting all the amplitudes equal to zero in the y list.

        # Set the second to last amplitude in the y array to zero (which it 
        # should be at already) and the last amplitude equal to a small 
        # negative value (-dx*1.e-4)

    # -------------------------------------------------------------------------
    # Method solve_matching: Solve the discretized Schrodinger equation for
    # the given energy (e). Note that we start the solution at large x and
    # then move toward x=0.
    # -------------------------------------------------------------------------

    def solve_matching(self,e):

        # Loop over the positions, starting at the third to last position
        # (we don't want to modify our initial conditions) and moving down to
        # x=0.

            # Calculate the "previous" (i-1) amplitude.

            # If the value of y gets too big, we need to renormalize
            # everything thus far in the array, otherwise we risk getting 
            # overflow errors. If the amplitude we just calculated is 
            # greater than 10, reduce all of the amplitudes calculated thus
            # far by a factor of 1.e-4.

    # -------------------------------------------------------------------------
    # Method find_matching: find the bound state energies up to level lmax,
    # noting that l=0 is the ground state. The additional parameters are the
    # initial energy step used in our search (deinit), the energy tolerance at 
    # which we terminate our search (tol).
    # -------------------------------------------------------------------------

    def find_matching(self,lmax,deinit,tol):

        # Set the initial energy equal to deinit. Knowing that the ground
        # state energy will be positive, this is a good place to start our
        # search.

        # Loop over states, starting at l=0 (ground state) and proceeding up to
        # and including lmax.

            # Set the energy step de equal to the proscribed initial value
            # deinit

            # Determine the parity of the solution. I use parity=0 to indicate
            # even parity and parity=1 to indicate odd parity.

            # Depending on the parity, we either need to keep track of the 
            # slopes at x=0 (even parity) or the amplitude at x=0 (odd parity)
            # The slope will be held in the variable "slope", and it's previous
            # value will be held in slope_last. When the slope at x=0 changes,
            # we know we have "crossed over" the "correct" energy. For now,
            # set these slopes equal to zero. We'll store the previous value
            # of y(0) in the variable "yz_last" (set it equal to zero for now).

            # Now, eep going up in energy until we know we are above
            # the new (next) solution. This will be indicated for even parity
            # solutions by the slope at x=0 becoming positive when y is 
            # negative, or negative when y is positive. 
            # Given the symmetry of the potential, the odd parity slopes will 
            # always match at x=0. What we need to check for odd parity is
            # whether y(0.0)=0.0. Since our search will start with y(0) being
            # positive, we terminate when y(0) becomes negative.

                # Stopping conditions depend on parity! If the parity is even,
                # break when the slope at x=0 multiplied by y(0) becomes
                # negative.

                # If the parity is odd, break when y(0) becomes negative.

                # Increase the energy by 2*de (makes the search a little faster)

                # Store the current y(0) in yz_last

                # Call the method initialize_matching to set up the arrays
                # and initial conditions.

                # Call the method solve_matching

                # Use the amplitudes at x=0 and x=dx to calculate the slope

                # Also save the slope in slope_last. We'll need that in the
                # next step.

            # Now refine the search, continuting until abs(de) is less than
            # the specified tolerance.

                # Call the method initialize_matching to set up the arrays
                # and initial conditions.

                # Call the method solve_matching

                # Use the amplitudes at x=0 and x=dx to calculate the slope

                # If the parity is even and the sign of slope is opposite the
                # sign of slope_last (their product will be negative), reduce
                # the step size by a factor of two and multiply it by -1, which
                # reverses our search direction.

                # If the parity is odd aht the sign of y(0) is opposite the
                # sign of yz_last (their product will be negative), reduce
                # the step size by a factor of two and multiply it by -1, which
                # reverses our search direction.

                # Add the step de to the energy

                # Copy the current slope into slope_last, and copy the current
                # y(0) into yz_last.

            # We have now converged on our solution! Let's run the 
            # initialization and solver one more time with energy e-de, since
            # our last iteration overshot the solution.

            # Normalize the solution by calling the normalize method, and 
            # then print the results to stdout. I'm formatting my printing
            # so that the wavefunctions can be easing plotted using gnuplot, 
            # and the energies extracted using the grep command.
            self.normalize()
            print(f"# Matching: n={self.n:2d}, state={l:2d}, E={e:10.3f}")
            for i in range(0,len(self.x)-1):
                print(self.x[i],self.y[i])
            print("\n\n")
