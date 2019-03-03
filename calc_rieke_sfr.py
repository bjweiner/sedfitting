
# calculate SFR by prescriptions of Rieke et al 2009
# also return the L_TIR.  This can take input vectors of z and logflux

# band tells what band the fluxes are observed in
# correctband allows correcting e.g. 24 um SFR via
# the fit of sfr24 vs sfr70

# logflux is assumed in Jy

# read coefficients out of a lookup table

# returns logsfr, loglir

# bjw 5/18/2016

import numpy as np
from astropy import ascii
import os

def calc_rieke_sfr(z, logflux, band, correctband = 0.0, rujo2013 = 0.0):

    homedir = os.getenv('HOME')
    fname = 'idl/sedfitting/z.sfr.mips+pacs.fitpars2.dat'
    fullfname = homedir = '/' + fname
    f = open(fullfname,'r')
    fdata = ascii.read(f,data_start=0,delimiter='\s')
    f.close()
    if rujo2013 == 1:
        fname = '/idl/sedfitting/z.sfr.mips24.rujopakarn2013.dat'
        fullfname = homedir = '/' + fname
        f = open(fullfname,'r')
        rdata = ascii.read(f,data_start=0,delimiter='\s')
        f.close()
        zrujo = rdata[0]
        arujomips24 = rdata[1]
        brujomips24 = rdata[2]

    zfits = fdata[0]
    afitsmips24 = fdata[1]
    bfitsmips24 = fdata[2]
    afitspacs70 = fdata[3]
    bfitspacs70 = fdata[4]
    afitspacs100 = fdata[5]
    bfitspacs100 = fdata[6]
    afitspacs160 = fdata[7]
    bfitspacs160 = fdata[8]
    
    if band == 24:
        if rujo2013 == 1:
            print "SFR using MIPS 24, Rujopakarn 2013 revision"
            afits = arujomips24
            bfits = brujomips24
            zfits = zrujo
        else:
            print "SFR using MIPS 24, Rieke 2009 original"
            afits = afitsmips24
            bfits = bfitsmips24
    elif band == 70:
        print "SFR using MIPS/PACS 70"
        afits = afitspacs70
        bfits = bfitspacs70
    elif band == 100:
        print "SFR using PACS 100"
        afits = afitspacs100
        bfits = bfitspacs100
    elif band == 160:
        print "SFR using MIPS/PACS 160"
        afits = afitspacs160
        bfits = bfitspacs160
    else:
        print "calc_rieke_sfr: dont have coefficients for band ",band
        afits = afitsmips24
        bfits = bfitsmips24

    a = np.interp(z,zfits,afits)
    b = np.interp(z,zfits,bfits)
    # Replace with astropy call
    dlumpc = lumdistance_lcdm(z)
    log4pidlsq = np.log10(4.0*3.1416) + 2 * (np.log10(dlummpc) + np.log10(3.086e24))

    if rujo2013 == 0:
        logsfr = a + b * (log4pidlsq + logflux - 53)

        # Apply correction based on the fit between 24 and 70 predictions
        # for galaxies in FIDEL data - Rieke2009 for both 24 and 70
        if (band == 24 and correctband == 70 and rujo2013==0):
            logsfr = 0.3745 + 0.6636 * logsfr

        loglir = 10.052 + 0.945 * logsfr
        ibright =  (logsfr > 1.003)
        # Check if this has any non-zero elements?
        loglir[ibright] = 10.096 + 0.902 * logsfr[ibright]

    else:
        # Rujopakarn 2013 uses a different equation, for loglir, 
        # which I can then turn into L24 and SFR, see eqns in the paper
        loglir = a + b * (log4pidlsq + logflux - 45)
        logl24 = (loglir - 1.096) / 0.982
        logsfr = logl24 - 9.108
        ibright = (loglir > 12.114)
        # Check if this has any non-zero elements?
        logsfr[ibright] = logl24[ibright] - 9.108 + 0.048 * (-11.208 + logl24[ibright])

    return logsfr, loglir


