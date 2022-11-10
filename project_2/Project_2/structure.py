"""
<Description of this module goes here: what it does, how it's used.>
"""

import numpy as np
from eos import pressure, density
from ode import rk4
from astro_const import G, Ke

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
    
    dzdm = np.zeros_like(z)#sets up an empty array for dzdm (change in z values over mass)

    r = z[0]#sets up an array of radius values which are the first elements of the dzdm array

    p = z[1]#sets up an array of pressure values which are the second elements of the dzdm array

    dzdm[0] = (4*np.pi*r**2*density(p,mue))**(-1) #uses the equation for radius to determine the first elements of the dzdm array

    dzdm[1] = -G*m/(4*np.pi*r**4)#uses the pressure equation to determine the second elements of the dzdm array

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
    z = np.zeros(2)#sets up an array of zeroes that has a shape of 2

    m = delta_m#sets mass as the delta_m input of the function

    z[1] = Pc#sets pc as the second elements of the z array

    rho = density(Pc,mue)#determines density by running the density function with inputs of pc (calculated above) and mue (given for central_values)

    z[0] = ((3*m)/(4*np.pi*rho))**(1/3)#calculates the first elements of the z array
    
    return z#returns the z array
    
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
    H_r = 4*np.pi*z[0]**3*density(z[1],mue)#calculates radius
    H_p = (4*np.pi*z[0]**4*z[1])/(G*m)#calculates pressure

    h = min(H_r,H_p)#determines whether radius or pressure is smaller and then assigns it to the value h (step size)

    return h#returns h
    
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
        
    m_step = np.zeros(max_steps)#creates an array of zeros for m_step
    r_step = np.zeros(max_steps)#creates an array of zeros for r_step
    p_step = np.zeros(max_steps)#creates an array of zeros for p_step

    # set starting conditions using central values
    z = central_values(Pc,delta_m,mue)#sets z by running central_values function
    m = delta_m#sets mass as delta_m (given when running the function)
    
    Nsteps = 0#Initializes step value as zero
    for step in range(max_steps):#runs a loop which runs for the value of max_steps (10,000 in this case)

        radius = z[0]#sets radius as the first elements of the z array
        pressure = z[1]#sets pressure as the second elements of the z array
        # are we at the surface?

        if (pressure < eta*Pc):#breaks the loop if pressure is too large
            
            break
        # store the step
        m_step[step] = m#sets m as the step value in the mass_step array (step 1, step 2, step 3, etc. = m_step[1],m_step[2],m_step[3],etc.)
        
        r_step[step] = radius#same as above but for radius and radius_step

        p_step[step] = pressure#same as above but for pressure and pressure_step

        # set the stepsize

        h = xi*lengthscales(m_step[step], z, mue)#sets stepsize for the actual taking of the step

        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=(mue))#calculates the new values for each step taken
        m += h#updates mass by adding step to the previous mass (m = m + h)

        # increment the counter
        Nsteps += 1#updates step counter
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')#raises an error if the number of steps exceeds the max set (10,000)
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]#returns all mass, radius, and pressure values over the whole loop in the m_step, r_step, and p_step arrays


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
    Pguess = (G**5/Ke**4)*(m*mue**2)**(10/3)#calculates pressure guess from equation 16 in instructions
    return Pguess
