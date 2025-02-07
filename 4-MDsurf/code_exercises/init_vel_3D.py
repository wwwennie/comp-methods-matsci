import numpy as np

# initialize velocities and momentums using Maxwell-Boltzmann distribution
# rescale velocities to prescribed temperature Tin

def init_vel_3D(n, Tin):
    """
    Pick velocities from Maxwell-Boltzmann distribution
    for any temperature we want.
    Then we will calculate the kinetic energy and thus
    the temperature of these atoms and then we will
    rescale the velocities to the correct temperature

    Input: 
       n (integer): number of steps in trajectory
       Tin (float): initial temperature
    Output:
       vx, vy, vz (float): initial velocities
       px, py, pz (float): initial momentums
    """
    k = 0
    px = 0
    py = 0
    pz = 0

    vx = np.zeros(n)
    vy = np.zeros(n)
    vz = np.zeros(n)

    for i in range(n):
        vx[i] = np.sqrt(-2 * np.log(np.random.rand())) * np.cos(2 * np.pi * np.random.rand())
        vy[i] = np.sqrt(-2 * np.log(np.random.rand())) * np.cos(2 * np.pi * np.random.rand())
        vz[i] = np.sqrt(-2 * np.log(np.random.rand())) * np.cos(2 * np.pi * np.random.rand())
        
        px += vx[i]
        py += vy[i]
        pz += vz[i]

    # Find average momentum per atom
    px /= n
    py /= n
    pz /= n

    # Set net momentum to zero and calculate K
    for i in range(n):
        vx[i] -= px
        vy[i] -= py
        vz[i] -= pz

        k += vx[i]**2 + vy[i]**2 + vz[i]**2

    k *= 0.5

    # Kinetic energy of desired temperature (Tin)
    kin = 1.5 * n * Tin

    # Rescale velocities
    sc = np.sqrt(kin / k)
    for i in range(n):
        vx[i] *= sc
        vy[i] *= sc
        vz[i] *= sc

    return vx, vy, vz, px, py, pz

