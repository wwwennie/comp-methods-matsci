import numpy as np

# 2D surface

def phisurf(acon, x, y):
    """
    Surface potential
    """
    phi = acon * np.sin(np.pi * x) * np.sin(np.pi * y)
    return phi

def fsurf(acon, x, y):
    """
    Surface force
    """
    fx = -acon * np.pi * np.cos(np.pi * x) * np.sin(np.pi * y)
    fy = -acon * np.pi * np.sin(np.pi * x) * np.cos(np.pi * y)
    return fx, fy


def initsurf(ein, acon):
    """
    Pick x and y so they sit near a well
    """
    x = -0.75 + 1.5 * np.random.rand()
    y = -0.75 + 1.5 * np.random.rand()

    # Calculate potential at x, y
    phi = phisurf(acon, x, y)

    # Find kinetic energy such that total energy is ein
    ki = ein - phi

    # Pick random velocity components
    vxs = np.random.rand()
    vys = np.random.rand()

    # Find kinetic energy from random velocities
    ks = 0.5 * (vxs**2 + vys**2)

    # Scale velocities so total energy is equal to ein
    vx = vxs * np.sqrt(ki / ks)
    vy = vys * np.sqrt(ki / ks)

    return x, y, vx, vy

