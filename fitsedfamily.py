
# Given a family of SED models and a single object's set of
# flux data points and redshift, find the best fit SED by calling
# fitonesed repeatedly.

# the SED family is a structure, or in this case a list where each
# SED is a dictionary with wave and flux tags

# returns array of normalization and chisq for each SED, and
# number of the best one

# returns nbest, fnorm array, chisq array. nbest is 0-based

# to replace fitsedprobs, optionally return the probability array,
# which is just computed from the chisq array and deg of freedom

# bjw 5/18/2016

import numpy as np
import matplotlib.pyplot as plt
import sedflux
import fitonesed
import scipy.special

def fitsedfamily(sedstruct, z, wobs, fobs, efobs, filtstruct, fixnorm=0.0, makeplot=0, logplot=0, probability=0):

    nsed = len(sedstruct)
    nobs = len(wobs)
    fnormarray = np.zeros(nsed)
    chisqarray = np.zeros(nsed)

    for i in range(nsed):
        wsed = sedstruct[i]['wave']
        fsed = sedstruct[i]['flux']
        # fnorm1, chisq1 = fitonesed.fitonesed(wsed, fsed, z, wobs, fobs, efobs, filtstruct, fixnorm=fixnorm, makeplot=makeplot,logplot=logplot)
        # Try using minimize=1 here?
        fnorm1, chisq1, fitresult = fitonesed.fitonesed(wsed, fsed, z, wobs, fobs, efobs, filtstruct, fixnorm=fixnorm, makeplot=0)
        # fnorm1, chisq1, fitresult = fitonesed.fitonesed(wsed, fsed, z, wobs, fobs, efobs, filtstruct, fixnorm=fixnorm, makeplot=0, minimize=1)
        fnormarray[i] = fnorm1
        chisqarray[i] = chisq1

    nbest = np.argmin(chisqarray)
    print "best fit SED is ",nbest,sedstruct[nbest]['label']," chisq ",chisqarray[nbest]
    
    # plot best fit SED
    if makeplot==1:
        zp1 = 1.0+z
        wsed = sedstruct[nbest]['wave']
        fsed = sedstruct[nbest]['flux']
        fnorm = fnormarray[nbest]
        nobs = np.size(fobs)
        fpredict = np.zeros(nobs)
        for i in range(nobs):
            fpredict[i] = sedflux.sedflux(wsed*zp1, fsed, filtstruct[i]['wave'], filtstruct[i]['response'] )
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

    if probability == 0:
        return nbest, fnormarray, chisqarray
    else:
        dof = nobs-1
        probarray = scipy.special.gammaincc((dof-1)/2.0, chisqarray/2.0)
        return nbest, fnormarray, chisqarray, probarray
    # done

    

