
# fitsedprobs is now deprecated because I changed fitsedfamily
# a lot and it's easier to just optionally return the probability
# array from fitsedfamily.  5/31/2016

# Given a (large) set of SED models and a single object's set of
# flux data points and redshift, find the best fit SED by calling
# fitonesed repeatedly. This is a small mod of fitsedfamily to
# also return the fit probability for each SED to allow marginalizing

# the SED family is a structure, or in this case a list where each
# SED is a dictionary with wave and flux tags

# returns array of normalization, probability, and chisq for each SED, and
# number of the best one

# returns nbest, fnorm array, prob array, chisq array. nbest is 0-based

# bjw 5/19/2016

import numpy as np
from scipy import special

def fitsedprobs(sedstruct, z, wobs, fobs, efobs, filtstruct, fixnorm=0.0, makeplot=0):

    nobs = len(wobs)
    nsed = len(sedstruct)
    fnormarray = np.zeros(nsed)
    chisqarray = np.zeros(nsed)
    probarray = np.zeros(nsed)

    for i in range(nsed):
        wsed = sedstruct[i]['wave']
        fsed = sedstruct[i]['flux']
        fnorm1, chisq1 = fitonesed(wsed, fsed, z, wobs, fobs, efobs, filtstruct, fixnorm=fixnorm, makeplot=makeplot)
        fnormarray[i] = fnorm1
        chisqarray[i] = chisq1

    # probabilities assuming chi-squared is meaningful
    probarray = special.gammaincc((nobs-2)/2.0, chisqarray/2.0)
    
    nbest = np.argmin(chisqarray)
    print "best fit SED is ",nbest,sedstruct[nbest]['label']," chisq ",chisqarray[nbest]


    # plot best fit SED?

    return nbest, fnormarray, probarray, chisqarray

