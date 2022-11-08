"""
<Description of this module goes here: what it does, how it's used.>
"""

import numpy as np
from eos import pressure, density # fill this in
from ode import rk4
from astro_const import G, Ke # fill this in

def stellar_derivatives(m,z,mue):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure)
        mue
            ratio, nucleons to electrons.  For a carbon-oxygen white dwarf, 
            mue = 2.
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm
    """
    
    dzdm = np.zeros_like(z)

    r = z[0]

    p = z[1]

    dzdm[0] = (4*np.pi*r**2*density(p,mue))**(-1)

    dzdm[1] = -G*m/(4*np.pi*r**4)
    
    return dzdm

def central_values(Pc,delta_m,mue):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        Pc
            central pressure (units = ?)
        delta_m
            core mass (units = ?)
        mue
            nucleon/electron ratio
    
    Returns
        z = array([ r, p ])
            central values of radius and pressure (units = ?)
    """
    z = np.zeros(2)

    m = delta_m

    z[1] = Pc

    rho = density(Pc,mue)

    z[0] = ((3*m)/(4*np.pi*rho))**(1/3)

    # compute initial values of z = [ r, p ]
    return z
    
def lengthscales(m,z,mue):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = ?)
        z (array)
           [ r, p ] (units = ?)
        mue
            mean electron weight
    
    Returns
        z/|dzdm| (units = ?)
    """

    # fill this in
    H_r = 4*np.pi*z[0]**3*density(z[1],mue)
    H_p = (4*np.pi*z[0]**4*z[1])/(G*m)

    h = min(H_r,H_p)

    return h
    
def integrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
    """
    Integrates the scaled stellar structure equations
    Arguments
        Pc
            central pressure (units = ?)
        delta_m
            initial offset from center (units = ?)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (what are the units?)
    """
        
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    z = central_values(Pc,delta_m,mue)
    m = delta_m
    
    Nsteps = 0
    for step in range(max_steps):

        radius = z[0]
        pressure = z[1]
        # are we at the surface?
        if (pressure < eta*Pc):
            break
        # store the step
        m_step[step] = m
        
        r_step[step] = radius

        p_step[step] = pressure

        # set the stepsize
        h = xi*lengthscales(delta_m, z, mue)
        
        # take a step
        znew = rk4(stellar_derivatives,m,z,h,args=(mue))
        m += h
        
        z = znew

        # increment the counter
        Nsteps += 1
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]

def pressure_guess(m,mue):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation. 
    
    Arguments
        m
            mass of white dwarf (units are ?)
        mue
            mean electron mass
    
    Returns
        P
            guess for pressure
    """
    # fill this in
    Pguess = (G**5/Ke**4)*(m*mue**2)**(10/3)
    return Pguess
