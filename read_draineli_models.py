
# read a list of Draine & Li (2007) model filenames and
# read each into a structure. Returns a list of structures (dictionaries)

# bjw 5/26/2016

import numpy as np
import os
import matplotlib.pyplot as plt
from read_one_draineli_model import *

def read_draineli_models(fname, dir='', makeplot=0):

    if fname == '':
        print 'read_draineli_models: need to supply a list filename'

    # Use dir to prepend to each of the names in the DL filelist,
    # not to the list name itself
    #if dir != '':
    #    fname = dir + '/' + fname
        
    flist = open(fname,'r')
    # test for failing to open file

    modelstructs = []
    line = flist.readline()
    while line != '':
        if dir != '':
            line = dir + '/' + line
        struct1 = read_one_draineli_model(line.strip())
        modelstructs.append(struct1)
        line = flist.readline()
    flist.close()
    nstruct = len(modelstructs)
    print 'Read ',nstruct,' Draine & Li models'

    # make plot
    if makeplot != 0:
        plt.clf()
        for i in range(nstruct):
            style = 'k-'
            plt.plot(np.log10(modelstructs[i]['sedwave']), np.log10(modelstructs[i]['sednupnu']), style)
        plt.show()


    return modelstructs

