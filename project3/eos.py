# MSU Hollow Low-Mass Star Society
# Elias Taira, Michael Bellaver, Joe Epley, Erin Syerson
# AST 304 Fall 2022
# Project 3

"""
Routines to compute an adiabatic equation of state.
"""

import numpy as np

def mean_molecular_weight(Z,A,X):
    """Computes the mean molecular weight for a fully ionized plasma with an 
    arbitrary mixture of species
    
    Arguments
        Z, A, X (either scaler or array)
            charge numbers, atomic numbers, and mass fractions
            The mass fractions must sum to 1
            
    Returns
        mu
            mean molecular weight
    """
    # arrays of the Z, A, and X values for the makeup of the star
    Zs = np.array(Z)
    As = np.array(A)
    Xs = np.array(X)
    assert np.sum(Xs) == 1.0

    # compute value of mean molecular weight
    mu = np.sum((Xs/As*(1 + Zs)))**(-1)
    
    return mu
    
def get_rho_and_T(P,P_c,rho_c,T_c):
    """
    Compute density and temperature along an adiabat of index gamma given a 
    pressure and a reference point (central pressure, density, and temperature).
    
    Arguments
        P (either scalar or array-like)
            value of pressure
    
        P_c, rho_c, T_c
            reference values; these values should be consistent with an ideal 
            gas EOS
    
    Returns
        density, temperature
    """

    # replace with computed values
    rho = rho_c*(P/P_c)**(3/5)
    T = T_c*(P/P_c)**(2/5)

    return rho, T
