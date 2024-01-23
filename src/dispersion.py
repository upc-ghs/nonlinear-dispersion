'''
Functions related to the 
calculation of dispersion
'''
import numpy as np 


def nonlinear_dispersion_function( peclet, delta, beta ):
    '''
    Function of the Peclet number employed to calculate nonlinear dispersion. 
    The output of this function needs to be multiplied by the aqueous molecular
    diffusion to obtain the dispersion coefficient.

    @params
        - peclet: the grain Peclet number v*d/D_aq.
        - delta : the ratio between the length of the pore-channel and 
                  their hydraulic radius. 
        - beta  : the exponent of the function describing incomplete mixing.

    @references
        Chiogna et al. (2010), Evidence of Compound-Dependent Hydrodynamic and 
        Mechanical Transverse Dispersion by Multitracer Laboratory Experiments
    '''

    return ( peclet**2 /(peclet + 2 + 4*delta**2 ) )**beta
