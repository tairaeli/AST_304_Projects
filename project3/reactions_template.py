"""
MSU Hollow Low-Mass Star Society
Elias Taira, Michael Bellaver, Joe Epley, Erin Syerson
AST 304 Project 3

Uses the PP chain to approximate the heating rate
"""

import numpy as np

def pp_rate(T,rho,XH,pp_factor=1.0):
    """
    Specific heating rate from pp chain hydrogen burning. Approximate rate 
    taken from Hansen, Kawaler, & Trimble.
    
    Arguments
        T, rho
            temperature [K] and density [kg/m**3]
        XH
            mass fraction hydrogen
        pp_factor
            multiplicative factor for rate
    Returns
        heating rate from the pp-reaction chain [W/kg]
    """
    
    # fill this out
    rate = pp_factor
    return rate
