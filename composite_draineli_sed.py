
# Take two sets of Draine & Li SEDs and add them with varying
# mixes using gamma, according to
#   j_nu = (1-gamma)*j_nu[umin,umin] + gamma*j_nu[umin,umax]
# so it's intended that
# sedset1 is the power law distrib of U
# sedset2 is for the diffuse medium and has a single U=Umin
# these are SEDs produced by convert_draineli_sed

# gamma is an array; for each pair, produce several SEDs for
# diff values of gamma

# indexes1 and indexes2 allow you to specify which of set2 matches
# which of set1, thus, if gamma has 4 values,
#   output[0] = gamma[0] * sedset1[indexes1[0]] + (1-gamma[0])*sedset2[indexes2[0]]
#   output[1] = gamma[1] * sedset1[indexes1[0]] + (1-gamma[1])*sedset2[indexes2[0]]
#  and so on

# interpolate onto the wavelength array of sedset1

# return a similar sedset structure

import numpy as np
import matplotlib.pyplot as plt

def composite_draineli_sed(sedset1, sedset2, gamma, indexes1=0, indexes2=0, makeplot=0):

    if indexes1==0:
        indexes1 = range(len(sedset1))
    if indexes2==0:
        indexes2 = range(len(sedset2))
    
    if len(indexes1) != len(indexes2):
        print "Warning: sed set lengths should match in composite_draineli_sed. Unpredictable results."
        nsed = min(len(indexes1),len(indexes2))
    else:
        nsed = len(indexes1)

    ngamma = len(gamma)

    sedstruct = []
    for i in range(nsed):
        for j in range(ngamma):
            gam = gamma[j]
            sed1 = sedset1[indexes1[i]]
            sed2 = sedset2[indexes2[i]]
            # Draine's wavelengths are decreasing so need to reverse arrays
            # for interp, try using [::-1] for reversed view of array
            sed2fluxinterp = np.interp(sed1['wave'][::-1], sed2['wave'][::-1],
                                       sed2['flux'], left=0.0, right=0.0)
            # plt.clf()
            # plt.plot(np.log10(sed1['wave']),np.log10(sed1['flux']),'k-')
            # plt.plot(np.log10(sed2['wave']),np.log10(sed2['flux']),'b-')
            # plt.plot(np.log10(sed1['wave']),np.log10(sed2fluxinterp),'bx')
            # plt.show()
            sedfluxnew = gam * sed1['flux'] + (1.0-gam)*sed2fluxinterp
            labelnew = gam*sed1['label'] + (1.0-gam)*sed2['label']
            namenew = sed1['name'] + '_gamma' + str(gam) + '_' + sed2['name']
            struct1 = {'label':labelnew, 'name':namenew, 'wave':sed1['wave'], 'flux':sedfluxnew}
            sedstruct.append(struct1)

    if makeplot != 0:
        plt.clf()
        for i in range(len(sedstruct)):
            style = 'k-'
            plt.plot(np.log10(sedstruct[i]['wave']), np.log10(sedstruct[i]['flux']), style)
        plt.show()

    return sedstruct

