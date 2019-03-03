
# plot some (wave,flux) data points in filters
# a SED (probably the best-fit SED from a fit) multiplied by a normalization
# and the SED binned into the filters

# much of this is copied out of fitonesed.py
# return the wave-observed and the normalized binned SED fluxes
# wave-obs, flux-sed-scaled

# bjw 5/29/2016

import numpy as np
import matplotlib.pyplot as plt
import sedflux

def plot_data_sed(sedstruct, nbest, normfac, z, wobs, fobs, efobs, filtstruct, logplot=0):

    sedbest = sedstruct[nbest]
    wsed = sedbest['wave']
    fsed = sedbest['flux'] * normfac
    
    zp1 = 1.0 + z
    nobs = np.size(fobs)
    fpredict = np.zeros(nobs)
    nfilt = len(filtstruct)
    if nfilt != nobs:
        print 'Warning - number of fluxes and filters dont match: ',nobs,nfilt

    # How do filter units change when redshifted?
    # fnu has units like erg/cm^2/s/Hz. or E/dt/dnu
    # should fnu_obs ~ fnu_rest / (1+z) ?  see also Peebles eq 13.59
    # at some point before I thought it should be
    #  fnu_obs ~ fnu_rest / (1+z)^3, not sure why
    # This is just affecting the weighting across the width of a
    # single filter, so in practice it should make little difference

    fpredict = np.zeros(nobs)
    for i in range(nobs):
        #fpredict[i] = sedflux.sedflux(wsed*zp1, fsed, filtstruct[i]['wave'], filtstruct[i]['response'] / zp1**3)
        fpredict[i] = sedflux.sedflux(wsed*zp1, fsed, filtstruct[i]['wave'], filtstruct[i]['response'] )

    fpredict = fpredict

    # Restrict plot to the approx range of the data. This logic isn't
    # perfect and the margins aren't equally spaced in a log plot
    dw = np.amax(wobs) - np.amin(wobs)
    wplotmin = np.amin(wobs) - 0.1 * dw
    wplotmax = np.amax(wobs) + 0.1 * dw
    df = np.amax(fobs) - np.amin(fobs)
    fplotmin = np.amin(fobs) - 0.15 * df
    fplotmax = np.amax(fobs) + 0.15 * df
    plt.clf()
    if logplot==1:
        # Prevent negative plot boundaries in log space
        wplotmin = max([wplotmin, 0.5*np.amin(wobs), 1.0e-6])
        fplotmin = max([fplotmin, 0.2*np.amin(fobs), 1.0e-6])
        wplot = np.log10(wobs)
        plt.plot(wplot,np.log10(fobs),'ko')
        plt.errorbar(wplot,np.log10(fobs),yerr=efobs/fobs/2.3026,fmt='ko')
        plt.plot(wplot,np.log10(fpredict),'rx')
        plt.plot(np.log10(wsed+zp1),np.log10(fsed),'r-')
        fig = plt.axis(np.log10([wplotmin,wplotmax,fplotmin,fplotmax]))
    else:
        plt.plot(wobs,fobs,'ko')
        plt.errorbar(wobs,fobs,yerr=efobs,fmt='ko')
        plt.plot(wobs,fpredict,'rx')
        plt.plot(wsed*zp1,fsed,'r-')
        fig = plt.axis([wplotmin,wplotmax,fplotmin,fplotmax])
    #plt.show()

    return wobs, fpredict
