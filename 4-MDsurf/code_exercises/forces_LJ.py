import numpy as np

def forces_LJ(a, n, x, y, z):
    """
    Simple lattice sum for force with cutoffs and 
    minimum image convention 

    We calculate force (fx, fy, fz), energy (u), and
    part of the pressure (w).

    """
    fx = np.zeros(n)
    fy = np.zeros(n)
    fz = np.zeros(n)
    u = 0
    w = 0

    for i in range(n - 1):  # Note limits
        ftx = 0
        fty = 0
        ftz = 0
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
                dphi = (2 / dist**12 - 1 / dist**6)
                ffx = dphi * a * dx / dist**2
                ffy = dphi * a * dy / dist**2
                ffz = dphi * a * dz / dist**2
                ftx += ffx
                fty += ffy
                ftz += ffz
                phi = (1 / dist**12 - 1 / dist**6)
                u += phi
                w += dphi

                # Add -f to sum of force on j
                fx[j] -= ffx
                fy[j] -= ffy
                fz[j] -= ffz

        # Sum up force on i (fi)
        fx[i] += ftx
        fy[i] += fty
        fz[i] += ftz

    # Need to multiply LJ by 4 and force and pressure by 24
    # Also need to correct sign in f
    u *= 4
    w *= 24
    fx *= -24
    fy *= -24
    fz *= -24

    return u, w, fx, fy, fz

