# Random walk on a 2D lattice

import numpy as np

def random_walk(nt,latt_type):
    """
    Random walk on a 2D lattice
    Inputs:
        nt [integer] = number of desired jumps (i.e., time steps)
        latt_type [string] = lattice geometry, current options 
                             are 'square'; TODO: 'triangle
    Outputs:
        rs2 [array] = square displacement at each time step
        x   [array] = x-coordinate at each time step
        y   [array] = y-coordinate at each time step
    """
    # array for x- and y-coordinates along hopping path
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    rs2 = np.zeros(nt+1)

    # particle starts at origin
    x[0] = 0 
    y[0] = 0
    rs2[0] = 0

    ## square lattice
    if latt_type == 'square':
        # create a list of random numbers from 1 to 4 with nt entries
        fd = np.floor(4 * np.random.rand(nt))
        # next two lines define the jumps on the square lattice:
        #   right, up, left, down
        delx = np.array([1, -1, 0, 0])
        dely = np.array([0, 0, 1, -1])
    else:
        raise ValueError("Lattice type not implemented! See random_walk.py")
         
    # loop over nt jumps, add the jump vector as generated randomly in fd
    #sum over nt jumps
    for j in range(nt):
        x[j+1] = x[j] + delx[int(fd[j])] # x position at j+1 jump
        y[j+1] = y[j] + dely[int(fd[j])] # y position at j+1 jump

        # square displacement position at j+1 jump in 2D
        rs2[j+1] = x[j+1]**2 + y[j+1]**2 

    return rs2, x, y
