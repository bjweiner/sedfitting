
# read Magdis 2012 SEDs = templates from a file
# return a structure,
# which is a list of dictionaries

# this uses a text file that I made from some save files
# I got from a Magdis website, which are for L_IR = 1 Lsun,
# hopefully in units of nu * Lnu, then converting
# them to be in units of flux at 10 Mpc like the Rieke SEDs
# this is derived from read_rieke_seds and since I formatted
# the file the same way, only difference is filename and number of SEDs
# (and number of wavelengths)

# if flag rujo2013 is set, then renormalize the luminosities
# according to Wiphu's stretch factor table.
# though this wouldn't make any sense for the Elbaz or Magdis seds which
# have no luminosity-color relation, unlike Rieke or Chary+Elbaz,
# so it's commented out.

# need to set sedstruct_label = replicate(11.0, nseds)
# although there are 10 SEDs, 9 MS for diff z ranges and 1 SB
# and they all get the same label

# bjw 5/18/2016

# Hard code locations for the moment

import numpy as np
import os
from astropy.io import ascii
import matplotlib.pyplot as plt

def read_magdis_seds(fname = '', outfile = '', makeplot=0):

    if fname == '':
        fname = 'magdis2012_ms1_9_sb.wlmic_fluxjy_10mpc.dat'
        direc = 'fidel/magdis2012/'
        homedir = os.getenv('HOME')
        fname = homedir + '/' + direc + fname
   
    nseds = 10
    #sedlabels = 9.75 + 0.25 * np.arange(nseds)
    sedlabels = 11.0 + np.zeros(nseds)

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

#    if rujo2013 == 1:
#        logstretchfac = [-0.118, 0.013, 0.173, 0.408, 0.717, 1.101,
#                         1.560, 2.095, 2.704, 3.388]
#        newlabels = sedlabels + logstretchfac
#        # This may not make a copy
#        tmpstruct = sedstruct
#        for i in range(nseds):
#            sedstruct[i]['label'] = str(newlabels[i])
#            sedstruct[i]['flux'] = tmpstruct[i]['flux'] * 10**logstretchfac[i]

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

        
    
