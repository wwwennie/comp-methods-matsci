import numpy as np

def MDLJ(nc,density, tin, nsteps, dt):
    """
    Initialize positions and velocities.
    Calculate some useful quantities.
    Calculate initial energy and forces.
    Now start the time stepping with the Verlet algorithm.
    
    Input:
         nc (integer): number of unit cells
         density (float): volume density in reduced units
         tin (float): intput temperature in reduced units
         nsteps (integer): total number of time steps
         dt (float): time step in reduced units
    """
    # Initialize positions and velocities
    n, x, y, z, vx, vy, vz = initLJMD(nc, tin)

    # Calculate some useful quantities
    vol = n / density
    a = vol**(1/3)
    rc = a / 2

    # Calculate initial energy and forces
    u, w, fx, fy, fz = forces_LJ(a, n, rc, s)

    # Initialize variables
    xold = np.zeros(n)
    yold = np.zeros(n)
    zold = np.zeros(n)
    xnew = np.zeros(n)
    ynew = np.zeros(n)
    znew = np.zeros(n)
    
    # Time series arrays
    un = np.zeros(nsteps)
    kn = np.zeros(nsteps)
    en = np.zeros(nsteps)
    tn = np.zeros(nsteps)
    pn = np.zeros(nsteps)

    # First find the positions at t-dt
    for i in range(n):
        xold[i] = x[i] - vx[i]*dt/a + 0.5*fx[i]*dt**2/a
        yold[i] = y[i] - vy[i]*dt/a + 0.5*fy[i]*dt**2/a
        zold[i] = z[i] - vz[i]*dt/a + 0.5*fz[i]*dt**2/a

    # Start the time steps
    for j in range(nsteps):
        k = 0
        # Find positions for time t + dt and velocities for time t
        for i in range(n):
            xnew[i] = 2*x[i] - xold[i] + fx[i]*dt**2/a
            ynew[i] = 2*y[i] - yold[i] + fy[i]*dt**2/a
            znew[i] = 2*z[i] - zold[i] + fz[i]*dt**2/a
            vx[i] = a*(xnew[i] - xold[i]) / (2*dt)
            vy[i] = a*(ynew[i] - yold[i]) / (2*dt)
            vz[i] = a*(znew[i] - zold[i]) / (2*dt)
            k += vx[i]**2 + vy[i]**2 + vz[i]**2

        k *= 0.5
        temp = 2*k / (3*n)

        # Create time series of values
        e = k + u
        un[j] = u / n
        kn[j] = k / n
        en[j] = e / n
        tn[j] = temp
        pn[j] = density*temp + w / (3*vol)

        # Reset positions for next time step
        for i in range(n):
            xold[i], yold[i], zold[i] = x[i], y[i], z[i]
            x[i], y[i], z[i] = xnew[i], ynew[i], znew[i]

        # Calculate force and energy at new positions for next cycle
        u, w, fx, fy, fz = forces_LJ(a, n, rc, x, y, z)

    return un, kn, en, tn, pn

