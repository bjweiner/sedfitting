
# return SED flux in a filter
# default to interpolating the filter curve onto the SED curve
# optionally specify a redshift?  but then would have to also use
# luminosity distance

# bjw 5/18/2016

import numpy as np

def sedflux(swave, sflux, fwave1, fresp1, z=0.0, interpsed=0):
    # trim out empty (zero) entries in the filter array, where either
    # the wave is 0;  or the response is small? or leave those in
    iwave = (fwave1 > 0.0)
    fwave = fwave1[iwave]
    fresp = fresp1[iwave]
    filtwmin = np.amin(fwave)
    filtwmax = np.amax(fwave)
    if (z > 1.0e-5 or z < -1.0e5):
        fwave = fwave * (1.0+z)
        fresp = fresp * (1.0+z)**3
        # do something with lum distance

    if interpsed == 0:
        frespnew = np.interp(swave, fwave, fresp, left=0.0, right=0.0)
        sedfilttot = np.sum( frespnew * sflux)
        filttot = np.sum( frespnew )
    else:
        sfluxnew = np.interp(fwave, swave, sflux)
        sedfilttot = np.sum( fresp * sfluxnew )
        filttot = np.sum (fresp)

    return sedfilttot/filttot

        
    
