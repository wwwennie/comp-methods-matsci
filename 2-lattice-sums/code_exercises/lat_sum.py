import numpy as np
import time # simple timer

def lat_sum1(a, n, s):
    """ Naive implementation of lattice sum with Lennard-Jones potential
        
        Input:
          a (float): cell length
          n (integer): number of atoms in simulation cell
          s (array): fractional coordinates of atom positions
        Output:
          ucell (float): total energy of simulation cell
    """
    start = time.process_time()

    ucell = 0 

    # sum over all atoms and divide by two 
    for i in range(n):
        for j in range(n):
            xij = s[j, 0] - s[i, 0]
            yij = s[j, 1] - s[i, 1]
            zij = s[j, 2] - s[i, 2]
            dist = a * np.sqrt(xij**2 + yij**2 + zij**2)
            if dist > 0:
                phi = 4 * (1 / dist**12 - 1 / dist**6)
            else:
                phi = 0
            ucell = ucell + phi
    ucell = ucell / (2 * n)
   
    end = time.process_time() - start
    return ucell, end

def lat_sum2(a, n, s):
    """ Avoid overcounting in lattice sums

        Input:
          a (float): cell length
          n (integer): number of atoms in simulation cell
          s (array): fractional coordinates of atom positions
        Output:
          ucell (float): total energy of simulation cell
    """
    start = time.process_time()

    ucell = 0
    # to avoid double counting, we change the indices over the loops
    # i will never equal j, so the if statement can be removed
    # no more factor of two!
    for i in range(n - 1):
        for j in range(i + 1, n):
            xij = s[j, 0] - s[i, 0]
            yij = s[j, 1] - s[i, 1]
            zij = s[j, 2] - s[i, 2]
            dist = a * np.sqrt(xij**2 + yij**2 + zij**2)
            phi = 4 * (1 / dist**12 - 1 / dist**6)
            ucell = ucell + phi
    ucell = ucell / n

    end = time.process_time() - start
    return ucell, end

def lat_sum3(a, n, rc, s):
    """ Inclusion of cutoff distance

        Input:
          a (float): cell length
          n (integer): number of atoms in simulation cell
          rc (float): cutoff distance
          s (array): fractional coordinates of atom positions
        Output:
          ucell (float): total energy of simulation cell
    """
    start = time.process_time()

    ucell = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            xij = s[j, 0] - s[i, 0]
            yij = s[j, 1] - s[i, 1]
            zij = s[j, 2] - s[i, 2]
            dist = a * np.sqrt(xij**2 + yij**2 + zij**2)
            if dist <= rc:
                phi = 4 * (1 / dist**12 - 1 / dist**6)
            else:
                phi = 0
            ucell = ucell + phi
    ucell = ucell / n

    end = time.process_time() - start
    return ucell, end

def lat_sum4(a, n, rc, c, s):
    """ Naive implementation of periodic boundary conditions

        Input:
          a (float): cell length
          n (integer): number of atoms in simulation cell
          rc (float): cutoff distance
          c (integer): number of periodic neighbors to left and right
          s (array): fractional coordinates of atom positions
        Output:
          ucell (float): total energy of simulation cell
    """
    start = time.process_time()
    ucell = 0
    for i in range(n):
        for j in range(n):
            for k in range(-c, c + 1):
                for l in range(-c, c + 1):
                    for m in range(-c, c + 1):
                        xij = k + s[j, 0] - s[i, 0]
                        yij = l + s[j, 1] - s[i, 1]
                        zij = m + s[j, 2] - s[i, 2]
                        dist = a * np.sqrt(xij**2 + yij**2 + zij**2)
                        if 0 < dist <= rc:
                            phi = 4 * (1 / dist**12 - 1 / dist**6)
                        else:
                            phi = 0
                        ucell = ucell + phi
    ucell = ucell / (2 * n)

    end = time.process_time() - start
    return ucell, end


def lat_sum5(a, n, rc, s):
    """ Implementation of minimum image convention 

        Input:
          a (float): cell length
          n (integer): number of atoms in simulation cell
          rc (float): cutoff distance
          s (array): fractional coordinates of atom positions
        Output:
          ucell (float): total energy of simulation cell
    """
    start = time.process_time()
    print("starting lattice sum")

    ucell = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            xij = s[j, 0] - s[i, 0]
            yij = s[j, 1] - s[i, 1]
            zij = s[j, 2] - s[i, 2]
            xij = xij - round(xij)
            yij = yij - round(yij)
            zij = zij - round(zij)
            dist = a * np.sqrt(xij**2 + yij**2 + zij**2)
            if 0 < dist <= rc:
                phi = 4 * (1 / dist**12 - 1 / dist**6)
            else:
                phi = 0
            ucell = ucell + phi
    ucell = ucell / n

    end = time.process_time() - start
    return ucell, end

