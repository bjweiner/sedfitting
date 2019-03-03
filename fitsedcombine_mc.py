
# Derived from fitsedfamily_mc.py

# Given an SED and some error bars, generate N_mc monte carlo
# realizations of the SED, fit models to each realization with the
# fitsedcombine approach - fitting a combination and trying to make it
# sparse -
# store the predicted best fit coefficients, and the chisq,
# which will typically be pretty low

# Then compute, over the realizations, the:
#  - mean and rms of each coefficient
#  - mean and rms of the sum of the coeffs, which should be total mass

# return a list that contains statistics of the mc fits

# bjw 5/31/2016

import numpy as np
import matplotlib.pyplot as plt
import fitsedcombine

def fitsedcombine_mc(sedstruct, z, wobs, fobs, efobs, filtstruct, nmonte, penalize=0.0, initguess=0, makeplot=0, logplot=0):

    nsed = len(sedstruct)
    nobs = len(wobs)
    nbestarray = []
    bestcoeffsarray = np.zeros([nmonte,nsed])
    meancoeffarray = np.zeros(nsed)
    rmscoeffarray = np.zeros(nsed)
    sumcoeffarray = np.zeros(nmonte)
    chisqbestarray = np.zeros(nmonte)
    nzerocoeffarray = np.zeros(nmonte)

    for i in range(nmonte):
        errs_in_sigma = np.random.normal(0.0,1.0,nobs)
        fluxtmp = fobs + efobs * errs_in_sigma
        fitcoeffs, fiterrors, chisq1 = fitsedcombine.fitsedcombine(sedstruct, z, wobs, fluxtmp, efobs, filtstruct, makeplot=0, penalize=penalize, initguess=initguess)
        bestcoeffsarray[i,0:] = fitcoeffs
        chisqbestarray[i] = chisq1
        sumcoeffarray[i] = np.sum(fitcoeffs)
        nzerocoeffarray[i] = sum(fitcoeffs < 1.0e-4)
        
    for j in range(nsed):
        meancoeffarray[j] = np.mean(bestcoeffsarray[0:,j])
        rmscoeffarray[j] = np.std(bestcoeffsarray[0:,j],ddof=1)

    meansumcoeff = np.mean(sumcoeffarray)
    rmssumcoeff = np.std(sumcoeffarray,ddof=1)
    meannzerocoeff = np.mean(nzerocoeffarray)
    rmsnzerocoeff = np.std(nzerocoeffarray,ddof=1)
                                    
    print "mean fit coeffs per SED: ", meancoeffarray
    print "rms of fit coeffs: ", rmscoeffarray
    print "mean, rms of sum of coeffs: ", meansumcoeff, rmssumcoeff
    print "mean, rms of num of zero coeffs: ", meannzerocoeff, rmsnzerocoeff

    chisqbest_mean = np.mean(chisqbestarray)
    chisqbest_rms = np.std(chisqbestarray, ddof=1)

    print "avg, rms of chisq(bestfit): ",chisqbest_mean,chisqbest_rms

    # make a plot?  eg a histogram of the values of the coeffs
    if makeplot == 1:
        plt.clf()
        # set bins?
        colorarray = ['blue','cyan','green','red','magenta','black']
        for j in range(nsed):
            # Only make plot if the SED had some non-zero coeffs
            if np.max(bestcoeffsarray[0:,j]) > 1.0e-6:
                if logplot ==1:
                    xcoeff = np.log10(bestcoeffsarray[0:,j])
                else:
                    xcoeff = bestcoeffsarray[0:,j]
                    icolor = j % len(colorarray)
                    histvals1, binvals1, patches1 = plt.hist(xcoeff, color=colorarray[icolor], histtype='step')
        if logplot==1:
            fig = plt.xlabel('log M / Msun')
        else:
            fig = plt.xlabel('M / Msun')
        fig = plt.ylabel('number of realizations')
        # plt.figlegend( (patches1[0], patches2[0]), ('point estimates', 'expect values'), 'upper right')
        plt.show()
        plt.clf()
        if logplot==1:
            xsumcoeff = np.log10(sumcoeffarray)
            fig = plt.xlabel('log M / Msun')
        else:
            xsumcoeff = sumcoeffarray
            fig = plt.xlabel('Mass M')
        fig = plt.ylabel('number of realizations')
        histvals2, binvals2, patches2 = plt.hist(xsumcoeff, color='black', histtype='step')
        plt.show()

    # first list is quantities in units of the normalization, eg Msun
    result_coeffs = {'meansum':meansumcoeff, 'rmssum':rmssumcoeff, 'meancoeffs':meancoeffarray, 'rmscoeffs':rmscoeffarray, 'meannzero':meannzerocoeff, 'rmsnzero':rmsnzerocoeff }
    # , 'nzerocoeffs':nzerocoeffarray }
    # second list has to do with fit quality
    result_probs = [chisqbest_mean, chisqbest_rms]
    # third dict is values for each realization
    mc_fits = {'bestcoeffs': bestcoeffsarray, 'sumcoeff': sumcoeffarray}
    
    return result_coeffs, result_probs, mc_fits




