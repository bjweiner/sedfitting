
# Read IR filters into a structure, which is just a list of
# dictionaries.

# Hard code the filter file locations for the moment

# bjw 5/18/2016

from astropy.io import ascii
import os
import numpy as np
import matplotlib.pyplot as plt

def read_ir_filters(outfile='', mipsonly=0, iras=0, makeplot=0):

    if mipsonly == 1:
        nfilt = 3
        labels = [ 24.0, 70.0, 160.0 ]
        fnames = ['/deep/kcorr/filters/mips24filt.resp.edit',
                  '/deep/kcorr/filters/pacs70filt.resp.edit', 
                  '/deep/kcorr/filters/pacs160filt.resp.edit']
    elif iras ==1:
        nfilt = 4
        labels = [ 12.0, 25.0, 60.0, 100.0 ]
        fnames = ['/deep/kcorr/filters/IRAS-IRAS.12mu.mic.dat',
                  '/deep/kcorr/filters/IRAS-IRAS.25mu.mic.dat',
                  '/deep/kcorr/filters/IRAS-IRAS.60mu.mic.dat',
                  '/deep/kcorr/filters/IRAS-IRAS.100mu.mic.dat']
    else:
        nfilt = 8
        labels = [ 8.0, 24.0, 70.0, 100.0, 160.0, 250.0, 350.0, 500.0 ]
        fnames = ['/deep/kcorr/filters/irac.ch4.resp.edit',
                  '/deep/kcorr/filters/mips24filt.resp.edit',
                  '/deep/kcorr/filters/pacs70filt.resp.edit',
                  '/deep/kcorr/filters/pacs100filt.resp.edit',
                  '/deep/kcorr/filters/pacs160filt.resp.edit',
                  '/deep/kcorr/filters/Herschel-SPIRE.PSW.250.dat',
                  '/deep/kcorr/filters/Herschel-SPIRE.PMW.350.dat',
                  '/deep/kcorr/filters/Herschel-SPIRE.PLW.500.dat']

    filtstruct = []
    homedir = os.getenv('HOME')
    # or
    # homedir = os.path.expanduser("~")
    for i in range(nfilt):
        f1 = open(homedir+fnames[i], 'r')
        fdata = ascii.read(f1,data_start=0,delimiter='\s')
        f1.close()
        filt1 = {'label': labels[i],
                 'name': str(labels[i]),
                 'wave': np.array(fdata.columns[0]),
                 'response': np.array(fdata.columns[1]) }
        filtstruct.append(filt1)

    # write filter structure to a file like a FITS table?
    if outfile != '':
        print 'Writing filter structure not implemented yet'

    # make a plot
    wmin = 1.0e6
    wmax = 0.0
    respmin = 1.0e6
    respmax = -1.0e6
    for i in range(nfilt):
        wmin = min(wmin, np.amin(filtstruct[i]['wave']))
        wmax = max(wmax, np.amax(filtstruct[i]['wave']))
        respmin = min(respmin, np.amin(filtstruct[i]['response']))
        respmax = max(respmax, np.amax(filtstruct[i]['response']))
    # could I just do this?
    # wmin = np.amin(filtstruct['wave'])
    # wmax = np.amax(filtstruct['wave'])
    # print wmin,wmax, respmin,respmax

    # make a plot
    if makeplot != 0:
        plt.clf()
        fig = plt.axis([wmin,wmax,respmin,respmax])
        for i in range(nfilt):
            style = 'k-'
            plt.plot(filtstruct[i]['wave'], filtstruct[i]['response'], style)
        plt.show()

    return filtstruct


    
        
