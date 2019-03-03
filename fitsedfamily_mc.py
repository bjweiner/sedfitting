
# Given an SED and some error bars, generate N_mc monte carlo
# realizations of the SED, fit models to each realization with the
# fitsedfamily approach - fitting one at a time to find best SED -
# store the predicted best fit normalization (ie point estimate),
# and also the expectation value and rms over the probabilities,
# ie the normalization marginalized over the models.
# Then compute, over the realizations, the:
#  - mean and rms of point estimates
#  - mean and rms of expectation values (unweighted)
#  - weighted mean of expectation values, weighted by 1/rms^2
#  - median of the (rms'es of expectation values over the probabilities)

# This allows comparing the means and the sizes of various kinds of
# error estimates.  So eg we can see if the typical
#   rms over the probabilities of one SED
# compares well with the rms of means over the realizations.

# return a list that contains statistics of the mc fits

# bjw 5/31/2016

import numpy as np
import matplotlib.pyplot as plt
import fitsedfamily

def fitsedfamily_mc(sedstruct, z, wobs, fobs, efobs, filtstruct, nmonte, fixnorm=0.0, makeplot=0, logplot=0):

    nsed = len(sedstruct)
    nobs = len(wobs)
    nbestarray = []
    bestnormarray = np.zeros(nmonte)
    expectarray = np.zeros(nmonte)
    rmsarray = np.zeros(nmonte)
    totprobarray = np.zeros(nmonte)
    chisqbestarray = np.zeros(nmonte)

    for i in range(nmonte):
        errs_in_sigma = np.random.normal(0.0,1.0,nobs)
        fluxtmp = fobs + efobs * errs_in_sigma
        nbest, fnormarray, chisqarray, probarray = fitsedfamily.fitsedfamily(sedstruct, z, wobs, fluxtmp, efobs, filtstruct, fixnorm=fixnorm, makeplot=0, probability=1)
        nbestarray.append(nbest)
        bestnormarray[i] = fnormarray[nbest]
        totprob = np.sum(probarray)
        totprobarray[i] = totprob
        chisqbestarray[i] = chisqarray[nbest]
        
        # hopefully the root-weightedmean-square formula is right
        fnormexpect = np.sum(fnormarray*probarray) / totprob
        fnormrms    = np.sqrt(sum( probarray * (fnormarray-fnormexpect)**2 ) /totprob)
        expectarray[i] = fnormexpect
        rmsarray[i] = fnormrms

    norm_mean = np.mean(bestnormarray)
    norm_rms  = np.std(bestnormarray, ddof=1)
    expect_mean = np.mean(expectarray)
    expect_rms  = np.std(expectarray, ddof=1)
    weights = 1.0 / rmsarray**2
    expect_wtmean = np.sum(weights*expectarray) / np.sum(weights)
    rms_margin_median = np.median(rmsarray)

    print "mean, rms of best fits: ", norm_mean, norm_rms
    print "mean, rms of expectation values: ", expect_mean, expect_rms
    print "weighted mean of expec values:   ", expect_wtmean
    print "median of marginalized rms'es: ", rms_margin_median

    totprob_mean = np.mean(totprobarray)
    totprob_rms = np.std(totprobarray, ddof=1)
    chisqbest_mean = np.mean(chisqbestarray)
    chisqbest_rms = np.std(chisqbestarray, ddof=1)

    print "avg of total prob and chisq(bestfit): ",totprob_mean,chisqbest_mean

    # make a plot?  eg a histogram of the values of bestnorm and expect
    if makeplot == 1:
        plt.clf()
        # set bins?
        if logplot ==1:
            xbestnorm = np.log10(bestnormarray) 
            xexpect = np.log10(expectarray)
            fig = plt.xlabel('log M / Msun')
        else:
            xbestnorm = bestnormarray
            xexpect = expectarray
            fig = plt.xlabel('M / Msun')
        histvals1, binvals1, patches1 = plt.hist(xbestnorm, color='blue', histtype='step')
        histvals2, binvals2, patches2 = plt.hist(xexpect, color='red', histtype='step')
        fig = plt.ylabel('number of realizations')
        plt.figlegend( (patches1[0], patches2[0]), ('point estimates', 'expect values'), 'upper right')
        plt.show()

    # first list is quantities in units of the normalization, eg Msun
    result_norm = [norm_mean, norm_rms, expect_mean, expect_rms, expect_wtmean, rms_margin_median]
    # second list has to do with fit quality
    result_probs = [chisqbest_mean, chisqbest_rms, totprob_mean, totprob_rms]
    # third dict is values for each realization
    mc_fits = {'nbest': nbestarray, 'best_norm': bestnormarray, 'marginal_norm': expectarray, 'marginal_rms': rmsarray}
    
    return result_norm, result_probs, mc_fits




