""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

<Team name, members go here>
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions import pp_rate

def central_thermal(m,r,mu):
    """ 
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
            mass in solar units
        r
            radius is solar units
        mu
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
    # fill this in

    Pc = 
    rhoc = 3/(4*np.pi) * m/r**3
    Tc = (G*m/r)*mu*m_u/(2*kB)
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
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
    #sets up an empty array for dzdm (change in z values over mass)
    dzdm = np.zeros_like(z)
    
    #sets up an array of radius values which are the first elements 
    #   of the dzdm array
    r = z[0]
    
    #sets up an array of pressure values which are the second elements 
    #   of the dzdm array
    p = z[1]

    #uses the equation for radius to determine the first elements 
    #   of the dzdm array
    dzdm[0] = (4*np.pi*r**2*density(p,mue))**(-1)

    #uses the pressure equation to determine the second elements 
    #   of the dzdm array
    dzdm[1] = -G*m/(4*np.pi*r**4)

    return dzdm #returns the dzdm array

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
    
    #sets up an array of zeroes that has a shape of 2
    z = np.zeros(2)

    #sets mass as the delta_m input of the function
    m = delta_m
    
    #sets pc as the second elements of the z array
    z[1] = Pc

    #determines density by running the density function with inputs of 
    #   pc (calculated above) and mue (given for central_values)
    rho = density(Pc,mue)

    #calculates the first elements of the z array
    z[0] = ((3*m)/(4*np.pi*rho))**(1/3)
    
    # returns the z array
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

    #calculates radius
    H_r = 4*np.pi*z[0]**3*density(z[1],mue)
    
    #calculates pressure
    H_p = (4*np.pi*z[0]**4*z[1])/(G*m)

    #determines whether radius or pressure is smaller and then assigns 
    #   it to the value h (step size)
    h = min(H_r,H_p)

    #returns h
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
        
    #creates an array of zeros for m_step, r_step and p_step
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)

    # set starting conditions using central_values
    z = central_values(Pc,delta_m,mue)
    #sets mass as delta_m (given when running the function)
    m = delta_m
    
    #Initializes step value as zero
    Nsteps = 0
    #runs a loop which runs for the value of max_steps (10,000 in this case)
    for step in range(max_steps):
        
        #sets radius as the first elements of the z array
        radius = z[0]
        #sets pressure as the second elements of the z array
        pressure = z[1]
        
        # are we at the surface?
        #breaks the loop if pressure is too large
        if (pressure < eta*Pc):
            break
            
        # store the step
        #sets m as the step value in the mass_step array 
        #       (step 1, step 2, step 3, etc. = m_step[1],m_step[2],m_step[3],etc.)
        m_step[step] = m
        #same as above but for radius and radius_step
        r_step[step] = radius
        #same as above but for pressure and pressure_step
        p_step[step] = pressure

        # set the stepsize
        h = xi*lengthscales(m_step[step], z, mue)

        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=(mue))
        #updates mass by adding step to the previous mass (m = m + h)
        m += h

        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]
