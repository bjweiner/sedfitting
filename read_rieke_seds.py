
# read Rieke 2009 SEDs = average templates from a file
# return a structure,
# which is a list of dictionaries
# primarily having wave and flux.  flux is as-observed at 10 Mpc,
# units are Jy

# if flag rujo2013 is set, then renormalize the luminosities
# according to Wiphu's stretch factor table.

# bjw 5/18/2016

# Hard code locations for the moment

import numpy as np
import os
from astropy.io import ascii
import matplotlib.pyplot as plt

def read_rieke_seds(fname = '', outfile = '', rujo2013 = 0, makeplot=0):

    if fname == '':
        fname = 'mrtable2.input.dat'
        direc = 'fidel/grieke_templ4/'
        homedir = os.getenv('HOME')
        fname = homedir + '/' + direc + fname
   
    nseds = 14
    sedlabels = 9.75 + 0.25 * np.arange(nseds)

    f1 = open(fname, 'r')
    fdata = ascii.read(f1, data_start=0, delimiter='\s')
    f1.close()

    # wavelengths in first column
    wave1 = fdata.columns[0]
    sedstruct=[]
    for i in range(nseds):
        flux1 = fdata.columns[i+1]
        name1 = str(sedlabels[i])
        label1 = sedlabels[i]
        sed1 = {'label': label1,
                'name': name1,
                'wave': wave1,
                'flux': flux1}
        sedstruct.append(sed1)

    if rujo2013 == 1:
        # If labels are strings, would need to iterate through loop
        # changing labels to floats and changing back
        logstretchfac = [-0.118, 0.013, 0.173, 0.408, 0.717, 1.101,
                         1.560, 2.095, 2.704, 3.388]
        newlabels = sedlabels[0:10] + logstretchfac
        
        tmpstruct = sedstruct[0:10]
        for i in range(len(newlabels)):
            sedstruct[i]['name'] = str(newlabels[i])
            sedstruct[i]['label'] = newlabels[i]
            sedstruct[i]['flux'] = tmpstruct[i]['flux'] * 10**logstretchfac[i]

    # write SED structure to a file like a FITS table?
    if outfile != '':
        print 'Writing SED structure not implemented yet'

    # make a plot
    if makeplot != 0 :
        plt.clf()
        for i in range(nseds):
            style = 'k-'
            plt.plot(np.log10(sedstruct[i]['wave']), np.log10(sedstruct[i]['flux']), style)
        plt.show()

    return sedstruct

        
    
