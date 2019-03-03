
# read a Draine & Li dust SED model from a file
# return a structure,
# as a python dictionary

# The models are normalized as 
# jnu = emissivity per H nucleon? and in the tables
# nuPnu = nu dP/dnu = frequency * power radiated per H nucleon per freq
# jnu = power radiated per H nucleon per freq per steradian

# Need to know mass of dust per H nucleon. See Table 3 of DL 2007
# which gives Mdust / M_H for each model, MW3.1_00 ... SMC, and says
# M_dust / M_gas = 1/1.36 M_dust / M_H.  (to compensate for helium)
# M_dust / M_H is 0.0100 to 0.0104 for the MW models,
#  0.00343 to 0.00359 for LMC, and 0.00206 for SMC.
# then divide by m_H = mass of one H atom.

# bjw 5/24/2016

# Hard code default file locations for the moment

import numpy as np
import os
import matplotlib.pyplot as plt
# from astropy.io import ascii

def read_one_draineli_model(fname = '', directory = '', outfile = '', makeplot = 0):

    if directory != '':
        fname = directory + '/' + fname
    if fname == '':
        fname = 'U1.00/U1.00_1e5_MW3.1_20.txt'
        directory = 'dustmass/draine_li_2007'
        homedir = os.getenv('HOME')
        fname = homedir + '/' + directory + '/' + fname
   
    f1 = open(fname, 'r')
    # test for failing to open file and return null 
    
    # fdata = ascii.read(f, data_start=0, delimiter='\s')

    # read some headers
    l = f1.readline()
    l = f1.readline()
    l = f1.readline()
    grainmodel = float(l.split()[0])
    l = f1.readline()
    a01 = float(l.split()[0])
    sigma1 = float(l.split()[1])
    bc1 = float(l.split()[2])
    bc1 = bc1*0.92
    l = f1.readline()
    a02 = float(l.split()[0])
    sigma2 = float(l.split()[1])
    bc2 = float(l.split()[2])
    bc2 = bc2*0.92
    l = f1.readline()
    umin = float(l.split()[0])
    umax = float(l.split()[1])
    beta = float(l.split()[2])
    l = f1.readline()
    uavg = float(l.split()[0])
    l = f1.readline()
    radfield = l.split()[0]
    l = f1.readline()
    powerh = float(l.split()[0])
    l = f1.readline()
    l = f1.readline()
    # read table of fluxes in filters as lists
    bandwave = []
    bandnupnu = []
    bandjnu = []
    bandname = []
    l = f1.readline()
    while (l != '' and len(l.split()) >= 2) :
        fields = l.split()
        # print fields
        bandwave.append(float(fields[0]))
        bandnupnu.append(float(fields[1]))
        bandjnu.append(float(fields[2]))
        bandname.append(fields[3])
        l = f1.readline()
    # read past the blank lines
    while (l == '' or len(l.split()) <=1):
        l = f1.readline()
    # skipping 2 header lines, then read SED as lists
    l = f1.readline()
    l = f1.readline()
    sedwave = []
    sednupnu = []
    sedjnu = []
    while (l != '' and len(l.split()) >=2):
        fields = l.split()
        # print fields
        sedwave.append(float(fields[0]))
        sednupnu.append(float(fields[1]))
        sedjnu.append(float(fields[2]))
        l = f1.readline()
    f1.close()
    bandwave = np.array(bandwave)
    bandnupnu = np.array(bandnupnu)
    bandjnu = np.array(bandjnu)
    sedwave_array = np.array(sedwave)
    sednupnu_array = np.array(sednupnu)
    sedjnu_array = np.array(sedjnu)
    
    modelstruct = {'filename': fname, 'grainmodel': grainmodel,
                   'a01': a01, 'sigma1': sigma1, 'bc1': bc1,
                   'a02': a02, 'sigma2': sigma2, 'bc2': bc2,
                   'umin': umin, 'umax': umax, 'beta': beta,
                   'uavg': uavg, 'radfield': radfield, 'powerh': powerh,
                   'bandwave': bandwave, 'bandnupnu': bandnupnu, 'bandjnu': bandjnu, 'bandname': bandname,
                   'sedwave': sedwave_array, 'sednupnu': sednupnu_array, 'sedjnu': sedjnu_array}

    # write model structure to a file like a FITS table?
    if outfile != '':
        print 'Writing model structure not implemented yet'

    # print len(bandwave),len(bandnupnu),len(sedwave),len(sednupnu)
    # print modelstruct['sedwave']
    # print modelstruct['sednupnu']
    
    # plot the model
    if makeplot != 0: 
        plt.clf()
        style = 'k-'
        plt.plot(np.log10(modelstruct['sedwave']), np.log10(modelstruct['sednupnu']), style)
        plt.show()

    return modelstruct

        
    
