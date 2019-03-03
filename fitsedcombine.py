
# Given a set of SED models and a single object's set of
# flux data points and redshift, find the best fit linear combination
# of the SED models. 
# This is derived from fitsedprobs/fitsedfamily but works differently 
# since it returns 1 result and coefficients of the models
# and doesn't use fitonesed

# the SED family is a structure, or in this case a list where each
# SED is a dictionary with wave and flux tags

# returns array of coefficients, errors?
# for the best fit

# penalize is the coefficient of the L1 norm and may require some
# experimentation: set it to of order the number of data points?
# omit penalize or set it to 0 to skip penalizing by the L1 norm

# returns coef array, err_coef array, chisq. nbest is 0-based

# bjw 5/19/2016

import numpy as np
import matplotlib.pyplot as plt
import sedflux
import make_sedflux_z
import scipy.special

from scipy import optimize

def fitsedcombine(sedstruct, z, wobs, fobs, efobs, filtstruct, fixnorm=0.0, makeplot=0, logplot=0, penalize=0.0, initguess=0):

    nsed = len(sedstruct)
    nobs = len(wobs)
    #fnormarray = np.zeros(nsed)
    #chisqarray = np.zeros(nsed)
    #probarray = np.zeros(nsed)

    # Make fluxes of each SED in the filters at the z
    # here we are assuming the filtstruct already matches the wobs
    # ie the filters input match the observations.
    # fluxstruct is actually a 2-d array (nsed,nfilt)
    fluxstruct = make_sedflux_z.make_sedflux_z(sedstruct, filtstruct, z)
    print fluxstruct[0,:]

    if initguess==0:
        initvals = np.zeros(nsed)
    else:
        normalize_fac = calc_norm_fluxes(fluxstruct, fobs)
        initvals = np.ones(nsed) * normalize_fac
        
    boundlist = [ (0.0, None) ] * nsed
    # no constraints
    # result = optimize.minimize(fluxresid,initguess,args=(fobs,efobs,fluxstruct),method='SLSQP')
    # constrain coeffs to be non-negative
    penaltyfactor = penalize
    fitargs = (fobs,efobs,fluxstruct,penaltyfactor)
    # opt = {'disp': True}
    opt = {'disp': False}
    result = optimize.minimize(fluxresid_penalized,initvals,args=fitargs,method='SLSQP',bounds=boundlist,options=opt)
    # BFGS doesn't respect bounds, only L-BFGS-B, TNC and SLSQP
    # result = optimize.minimize(fluxresid_penalized,initvals,args=fitargs,method='L-BFGS-B',bounds=boundlist,options=opt)
    # result = optimize.minimize(fluxresid_penalized,initvals,args=fitargs,method='TNC',bounds=boundlist,options=opt)

    # do something with the result
    fitcoeffs = result['x']
    fluxbest = np.dot(fitcoeffs,fluxstruct)
    residbest = fobs-fluxbest
    #residvar = np.sum(residbest**2) / (nfilt-nsed)
    chisq = np.sum(((fobs-fluxbest)/efobs)**2)
    # can we get the error estimates from the Jacobian J or Hessian?
    # The covariance matrix should be something like (J^T J)^-1
    # if we minimized ((obs-model)/err)**2 then the covar matrix
    # gives the errors.
    jacob = result['jac']
    print "jacob = ",jacob
    # jacob from SLSQP is a ncoeff+1 element 1-d array, not a 2-d array
    # jacob from TNX is a ncoeff element 1-d array.  These should be
    # the gradients of the function being minimized.
    jacobt = np.transpose(jacob)
    # print "jacobt = ",jacobt
    # covar = np.linalg.inv( np.multiply(jacobt, jacob) )
    # For a least squares minimization the variance should be basically
    # covar = (J^T J)^-1.  For 1-D J vector, covar_i = 1/ J_i^2
    # so loosely, error_i = 1 / J_i.  This isn't giving me
    # sensible numbers though.
    fiterrors = np.zeros(nsed)
    ii = range(nsed)
    # fiterrors[ii] = np.sqrt(covar[ii,ii])
    fiterrors[ii] = 1.0 / abs(jacob[ii])
    # fiterrors = np.zeros(nsed)

    # one at a time fitting
    # for i in range(nsed):
    #    wsed = sedstruct[i]['wave']
    #    fsed = sedstruct[i]['flux']
    #    fnorm1, chisq1 = fitonesed(wsed, fsed, z, wobs, fobs, efobs, filtstruct, fixnorm=fixnorm, makeplot=makeplot)
    #    fnormarray[i] = fnorm1
    #    chisqarray[i] = chisq1

    # probabilities assuming chi-squared is meaningful
    # probarray = special.gammaincc((nobs-2)/2.0, chisqarray/2.0)
    
    # nbest = np.argmin(chisqarray)
    # print "best fit single SED is ",nbest,sedstruct[nbest]['label']," chisq ",chisqarray[nbest]


    # plot best fit SED
    if makeplot==1:
        zp1 = 1.0+z
        wsed = sedstruct[0]['wave']
        fsed = np.zeros(len(wsed))
        # fpredict = np.zeros(nobs)
        fluxpointspredict = fluxbest
        for i in range(nsed):
            # fpredict[i] = sedflux.sedflux(wsed*zp1, fsed, filtstruct[i]['wave'], filtstruct[i]['response'] )
            fsed = fsed + fitcoeffs[i] * sedstruct[i]['flux']
        plt.clf()
        if logplot==1:
            wplot = np.log10(wobs)
            plt.plot(wplot,np.log10(fobs),'ko')
            plt.errorbar(wplot,np.log10(fobs),yerr=efobs/fobs/2.3026,fmt='ko')
            plt.plot(wplot,np.log10(fluxpointspredict),'rx')
            plt.plot(np.log10(wsed),np.log10(fsed),'r-')
        else:
            plt.plot(wobs,fobs,'ko')
            plt.errorbar(wobs,fobs,yerr=efobs,fmt='ko')
            plt.plot(wobs,fluxpointspredict,'rx')
            plt.plot(wsed,fsed,'r-')
        plt.show()

    print "fobs = ",fobs
    print "bestmod = ", fluxbest
    print "coeffs = ", fitcoeffs
    print "c-errors = ", fiterrors
    
    return fitcoeffs, fiterrors, chisq
    #return nbest, fnormarray, probarray, chisqarray


# Function to minimize - chi-squared of data vs linear comb of models

def fluxresid(x, fobs, efobs, fluxstruct):
    #fobs, efobs, fluxstruct = args
    # the i'th x is the coefficient of the i'th SED's fluxes
    # this should be doable as a matrix multiplication to get a
    # vector fluxtot? depends on how fluxstruct is ordered
    # if fluxstruct is a 2-d array with SED as first coordinate
    # and flux as second, this should return a vector with length
    # nfilters
    fluxtot = np.dot(x,fluxstruct)
    
    chisq = np.sum( ((fobs-fluxtot)/efobs)**2 )
    return chisq

# Function to minimize - chi-squared of data vs linear comb of models
# penalized by some function of the coeffs, like L1 norm

def fluxresid_penalized(x, fobs, efobs, fluxstruct, penaltyfactor):
    #fobs, efobs, fluxstruct, penaltyfactor = args
    # the i'th x is the coefficient of the i'th SED's fluxes
    # this should be doable as a matrix multiplication to get a
    # vector fluxtot? depends on how fluxstruct is ordered
    # if fluxstruct is a 2-d array with SED as first coordinate
    # and flux as second, this should return a vector with length
    # nfilters
    fluxtot = np.dot(x,fluxstruct)
    
    chisq = np.sum( ((fobs-fluxtot)/efobs)**2 )

    norm = 'L1'
    # Summing the (absolute value of the) coeffs is the L1 norm.
    # Summing the sqrts of the coeffs is an L0.5 norm and penalizes
    # N small coeffs more than 1 big coeff.
    # coeffsum = np.sum(np.abs(x))
    if norm == 'L1':
        coeffsum = np.sum(x)
        penaltyterm = penaltyfactor * coeffsum
    elif norm == 'L0.5':
        coeffsumsqrt = np.sum(np.sqrt(abs(x)))
        penaltyterm = penaltyfactor * coeffsumsqrt
    else:
        coeffsumsq = np.sum(x**2)
        penaltyterm = penaltyfactor * coeffsumsq

    # print 'fluxresid_penalized: ',chisq,penaltyterm
    penalized = chisq + penaltyterm
    return penalized

# Order of magnitude to normalize fluxes to observations

def calc_norm_fluxes(fluxstruct, fobs):

    nfilt = len(fobs)
    nseds = len(fluxstruct)
    tmp = np.ones(nseds)
    tmpsum = np.dot(tmp, fluxstruct)
    meanfac = np.mean( fobs/tmpsum)

    return meanfac
