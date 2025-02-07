import numpy as np
from lattice_2D import random_walk

def avg_rand_walk(nt, nd, latt_type):
    """ Average over many trajectories of random walkers
        Input:
	    nt (integer) = number of jumps
	    nd (integer) = number of trials 
            latt_type (string) = type of lattice,
                                 "square" or "triangular" (TODO) implemented
        Output:
 	    rwa (list; float) = mean square displacement of each trajectory
            ree (list; float) = end-to-end distance for each trial
            sig (list; float) = relative standard deviation of <r^2>
                                of each trajectory
    """

    # Initialize the variables for calculating <R^2> and <R^4>
    rwa = np.zeros(nt+1)
    ree = np.zeros(nd)
    sig = np.zeros(nt+1)

    # Loop over trials
    for j in range(nd):
        # Call random walker code from above to generate a trajectory
        g, x, y = random_walk(nt,latt_type)

	# Calculate end-to-end distance, 2D
	# index last element
        ree[j] = np.sqrt(x[-1]**2 + y[-1]**2 )

        # Increment R^2 and R^4 at each time step (jump)
        for k in range(nt+1):
            rwa[k] += g[k]
            sig[k] += g[k] ** 2

    # Find averages by dividing by the number of trials
    # Note: rwa[0] is 0 by definition, will throw a NaN warning 
    for k in range(nt+1):
        rwa[k] /= nd
        sig[k] = (sig[k] / float(nd) - rwa[k] ** 2) / rwa[k] ** 2

    return rwa, ree, sig

