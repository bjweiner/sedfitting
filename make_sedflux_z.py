
# Given an SED structure, a filter structure, and a z, make
# a structure of the log fluxes of the SEDs in the filters at that z
# derived from make_sed_tracks

# refdist is the fiducial distance that the SED fluxes are normalized
# to, in Mpc; 10 Mpc for Rieke et al 2009

# returns fluxstruct (list of dictionaries?
# or should I return a list of lists here, or an array)
# returning a 2-d array

# bjw 5/20/2016

import numpy as np
import matplotlib.pyplot as plt
import lumdistance_lcdm
import sedflux

def make_sedflux_z(sedstruct, filtstruct, z, refdist=10.0, return_log=0):
    
  nfilt = len(filtstruct)
  nsed = len(sedstruct)
  dummy = ''
  log4pi = np.log10(4.0 * 3.141593)

  fluxstruct = np.zeros((nsed,nfilt))

  fluxes = np.zeros(nfilt)
  # logcolors = np.zeros(nfilt-1)

  lumdist = lumdistance_lcdm.lumdistance_lcdm(z)
  # convert to Mpc if needed, refdist is 10 Mpc usually

  logdistfac = 2.0 * np.log10(lumdist/refdist)
  for k in range(nsed):
    for j in range(nfilt):
      # should be multiplying filter (?) by 1+z here for the frequency units?
      fluxes[j] = sedflux.sedflux(sedstruct[k]['wave']*(1.0+z),
                          sedstruct[k]['flux'],
                          filtstruct[j]['wave'],
                          filtstruct[j]['response'] * (1.0+z)) 
    logfluxes = np.log10(fluxes) - logdistfac
    fluxes = 10**logfluxes
    #logcolors = logfluxes[0:nfilt-1] - logfluxes[1:nfilt]
    if return_log != 0:
      #struct1 = {'logfluxes': logfluxes}
      struct1 = logfluxes
    else:
      #struct1 = {'fluxes': fluxes}
      struct1 = fluxes
    #fluxstruct.append(struct1)
    fluxstruct[k] = struct1

  # make a plot?
  
  return fluxstruct

