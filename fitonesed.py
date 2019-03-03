
# fit observed data (wavelength, flux in f_nu, error arrays),
# given a redshift, with a single specified model SED
# return chisq of fit and normalization of SED

# fixnorm allows fixing the normalization of the SED, ie normalization
# isn't free, so fnorm = 1.0

# need to pass in the filters in order matching the wavelengths
# somehow in filterstruct.  It could be a pair of 2-d arrays, but
# that's clunky. could have each filter as a dictionary with keys
# wave and response and pass in a list of these dictionaries

# returns: fnorm, chisq

# bjw 5/18/2016

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import sedflux

def fitonesed(wsed, fsed, z, wobs, fobs, efobs, filtstruct, fixnorm = 0.0, makeplot=0, minimize=0, logplot=0):
    zp1 = 1.0 + z
    nobs = np.size(fobs)
    fratio = np.zeros(nobs)
    efratio = np.zeros(nobs)
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
    fratio = fobs / fpredict 
    efratio = efobs / fpredict

    # Quick fitting of normalization by taking the mean or wt mean of
    # ratios of SED to model in each filter
    
    if fixnorm > 0.0:
        fnorm = fixnorm
        result = fnorm
    elif minimize == 0:
        # unweighted mean of ratios
        fratio_mean = np.sum(fratio) / nobs
        # weighted mean of ratios
        weight = 1.0 / efratio**2
        fratio_wtmean = np.sum(fratio*weight) / np.sum(weight)
        fnorm = fratio_wtmean
        result = fnorm
    else:
        initguess = np.sum(fratio) / nobs
        result = optimize.minimize(fluxdiff_chisq,initguess,args=(fobs,efobs,fpredict),method='SLSQP')
        fnorm = result['x']
        if result['success'] != True:
            print "Fit failed in fitonesed?"
            print result['success'],result['status'],result['message']
        
    # compute chi-sq
    fratio = fratio / fnorm
    efratio = efratio / fnorm
    chisq = np.sum( ((fratio-1.0)/efratio)**2 )

    # print wobs
    # print fobs
    # print fnorm*fpredict

    # make a plot of observed and model points
    if makeplot==1:
        plt.clf()
        if logplot==1:
            wplot = np.log10(wobs)
            plt.plot(wplot,np.log10(fobs),'ko')
            plt.errorbar(wplot,np.log10(fobs),yerr=efobs/fobs/2.3026,fmt='ko')
            plt.plot(wplot,np.log10(fnorm*fpredict),'rx')
            plt.plot(wplot,np.log10(fnorm*fpredict),'r-')
        else:
            plt.plot(wobs,fobs,'ko')
            plt.errorbar(wobs,fobs,yerr=efobs,fmt='ko')
            plt.plot(wobs,fnorm*fpredict,'rx')
            plt.plot(wobs,fnorm*fpredict,'r-')
        plt.show()

    return fnorm, chisq, result

# function to minimize - difference between model and predicted
                                   
def fluxdiff_chisq_old(x,args):
    fobs, efobs, fpredict = args
    chisq = np.sum( ((fobs - x*fpredict)/efobs)**2 )
    return chisq

def fluxdiff_chisq(x,fobs, efobs, fpredict):
    chisq = np.sum( ((fobs - x*fpredict)/efobs)**2 )
    return chisq

        
