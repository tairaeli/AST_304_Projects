"""
defines 2 functions that find the pressure and density
"""

from astro_const import h, m_e, m_u
import numpy as np

def pressure(rho, mue):
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio
    
    Returns
        electron degeneracy pressure (Pascal)
    """
    
    # replace following lines with body of routine
    p = 1/5 * (3/(8*np.pi))**(2/3)*h**2/(m_e)*(rho/(mue*m_u))**(5/3)
    return p

def density(p, mue):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio
        
    Returns
        mass density (kg/m**3)
    """
    
    # replace following lines with body of routine
    rho = (p*(1/5 * (3/(8*np.pi))**(2/3) * h**2/m_e)**-1)**(3/5) * mue * m_u
    return rho