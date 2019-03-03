
# return luminosity distance in Mpc for an array of z values

# this is the simple way where I read it off a pre existing look up table

# should be able to replace this with an astropy.cosmology call

import numpy as np
# from astropy.read import ascii
import os

def lumdistance_lcdm(zarray):

    homedir = os.getenv('HOME')
    fname = 'idl/sedfitting/distcalc.out.lcdm.dz001'
    fullfname = homedir + '/' + fname

    f = open(fullfname,'r')
    # fdata = ascii.read(f,data_start=0,delimiter='\s')
    fdata = np.genfromtxt(fullfname)
    f.close()

    #zgrid = fdata.columns[0]
    #dcgrid = fdata.columns[1]
    #dmgrid = fdata.columns[2]
    #dagrid = fdata.columns[3]
    #dlgrid = fdata.columns[4]
    #dmodgrid = fdata.columns[5]
    #dvcgrid = fdata.columns[6]
    #vcgrid = fdata.columns[7]
    #tlookgrid = fdata.columns[8]

    zgrid = fdata[0:,0]
    dcgrid = fdata[0:,1]
    dmgrid = fdata[0:,2]
    dagrid = fdata[0:,3]
    dlgrid = fdata[0:,4]
    dmodgrid = fdata[0:,5]
    dvcgrid = fdata[0:,6]
    vcgrid = fdata[0:,7]
    tlookgrid = fdata[0:,8]

    dlarray = np.interp(zarray, zgrid, dlgrid)
    return dlarray


