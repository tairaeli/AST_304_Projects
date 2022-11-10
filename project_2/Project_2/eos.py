########################################################################
# MSU Hollow White Dwarf Society: Erin Syerson, Elias Taira, Michael 
#   Bellaver, and Joey Epley
# AST304, Fall 2020
# Michigan State University
########################################################################

"""
defines 2 functions that find the pressure and density

utilizes equation (1) in the instructions to make the calculations
    for pressure and denstity
"""

from astro_const import h, m_e, m_u
import numpy as np

def pressure(rho, mue):
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio (=2)
    
    Returns
        electron degeneracy pressure (Pascal)
    """
    
    # equation (1) from instructions
    p = 1/5 * (3/(8*np.pi))**(2/3)*h**2/(m_e)*(rho/(mue*m_u))**(5/3)
    return p

def density(p, mue):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio (=2)
        
    Returns
        mass density (kg/m**3)
    """
    
    # a rearangement of equation (1) from the instructions
    rho = (p*(1/5 * (3/(8*np.pi))**(2/3) * h**2/m_e)**-1)**(3/5) * mue * m_u
    return rho
