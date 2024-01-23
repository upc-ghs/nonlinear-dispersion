'''
Module storing analytical solutions
'''

import numpy as np 
import scipy.special as spc 

def domenico(
        xposition,
        yvector,
        injectionwidth,
        dispersiony,
        velocity,
        bounded=True,
        domainwidth=0,
    ):
    '''
    Analytical solution for the problem of two-dimensional 
    transport through a homogeneous porous medium (Srinivasan et al., 2007)

    @params
        - xposition      : longitudinal distance from the injection.
        - yvector        : transverse points for evaluating the solution, centered at zero.
        - injectionwidth : the injection width.
        - dispersiony    : the transverse dispersion coefficient.
        - velocity       : flow velocity.
        - bounded        : apply boundary correction to the solution, requires domainwidth.
        - domainwidth    : the width of the domain for applying the boundary correction.

    @note: validity of the solution is considered for the following conditions: 
        - t > 5  alpha_L/v
        - x > 30 alpha_L  

    @references:
        - Srinivasan et al. (2007), Domenico Solution? Is It Valid?
    '''

    output = 0.5*( 
            spc.erf( ( yvector + 0.5*injectionwidth )/( 2*np.sqrt(dispersiony*xposition/velocity) ) ) - 
            spc.erf( ( yvector - 0.5*injectionwidth )/( 2*np.sqrt(dispersiony*xposition/velocity) ) )
        )

    if bounded:
        if domainwidth==0:
            raise Exception( 'analytical.py: domenico: requires domainwidth>0 as keyword argument when bounded=True' )

        output += 0.5*( 
                spc.erf( ( yvector + 0.5*injectionwidth + domainwidth )/( 2*np.sqrt(dispersiony*xposition/velocity) ) ) - 
                spc.erf( ( yvector - 0.5*injectionwidth + domainwidth )/( 2*np.sqrt(dispersiony*xposition/velocity) ) ) + 
                spc.erf( ( yvector + 0.5*injectionwidth - domainwidth )/( 2*np.sqrt(dispersiony*xposition/velocity) ) ) - 
                spc.erf( ( yvector - 0.5*injectionwidth - domainwidth )/( 2*np.sqrt(dispersiony*xposition/velocity) ) )  
            )

    return output




if __name__ == '__main__':
    import matplotlib.pyplot as plt

    xposition      = 25
    domainwidth    = 5
    npoints        = 100 
    deltay         = domainwidth/npoints
    yvector        = domainwidth*0.5 - np.arange( 0, domainwidth + deltay, deltay )  
    injectionwidth = 0.4
    alphat         = 0.01
    velocity       = 30
    dispersiony    = alphat*velocity

    solution = domenico(
            xposition,
            yvector,
            injectionwidth, 
            dispersiony, 
            velocity, 
            bounded=False,
        )
    plt.plot( yvector, solution )

    solution = domenico(
            xposition,
            yvector,
            injectionwidth, 
            dispersiony, 
            velocity, 
            bounded=True,
            domainwidth=domainwidth
        )
    plt.plot( yvector, solution )


    plt.show(block=False)
    import pdb
    pdb.set_trace()
