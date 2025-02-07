import numpy as np

def gr(a, n, x, y, z, nbin):
    """
    Calculate radial distribution function.

    Input:
        a (float): simulation cell dimension
        n (integer): number of atoms
        x, y, z (array of floats): atomic positions
        nbin (integer): number of bins
    Output:
        bp (array): binning used in g(r)
        ng (array): frequencies corresponding to g(r)
    """
    # cutoff, bin size
    rc = a / 2
    xb = rc / nbin
    g = np.zeros(nbin)
    bp = np.zeros(nbin)
    ng = np.zeros(nbin)

    for i in range(n - 1):  # Note limits 
        for j in range(i + 1, n):  # Note limits
            # Minimum image convention
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            dz = z[j] - z[i]
            dx -= round(dx)
            dy -= round(dy)
            dz -= round(dz)
            dist = a * np.sqrt(dx**2 + dy**2 + dz**2)
            if dist <= rc:
                ib = int(dist / xb)
                if ib < nbin:
                    g[ib] += 1
                    
    # For computing g(r)
    # Normalize and create proper distances
    factor = 2 * a**3 / (4 * np.pi * n**2 * xb)
    for i in range(nbin):
        bp[i] = (i + 0.5) * xb
        ng[i] = factor * g[i] / ((i * xb)**2) if i != 0 else 0

    return bp, ng
