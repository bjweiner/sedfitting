
# make color and luminosity tracks of SEDs

# zlims is an array of zmin, zmax, dz like [0, 2, 0.05]

# refdist is the fiducial distance that the SED fluxes are normalized
# to, in Mpc; 10 Mpc for Rieke et al 2009

# this returns a structure that has tracks with redshift:
# sedlogflux is the observed log flux of each sed in each filter
# sedcolor is the observed log flux ratio of the Nth to the
# N+1'th filter
# ie if you pass it filters 24, 70, 160 um, you should get
# 3 logfluxes in 24, 70, 160, and 2 log ratios of 24/70 and 70/160

# See ~/data/bootes/brownmatch/Readme.fitsed
# for some use of the IDL version

# returns trackstruct (list of dictionaries?)
# This is not very convenient and I should think about returning
# 2-d arrays or something - decided to have make_sedflux_z
# return a 2-d array, so this could be a 3-d array, or a list
# of 2-d arrays?

# bjw 5/20/2016

import numpy as np

def make_sed_tracks(sedstruct, filtstruct, zlims, refdist=10.0, return_log=1):
    
  zmin = zlims[0]
  zmax = zlims[1]
  dz = zlims[2]
  nz = int((zmax-zmin)/dz) + 1
  zarray = zmin + dz * np.range(nz)
  nfilt = len(filtstruct)
  nsed = len(sedstruct)
  dummy = ''
  log4pi = np.log10(4.0 * 3.141593)

  trackstruct = []
  #fluxstruct = np.zeros(nsed,nfilt)

  fluxes = np.zeros(nfilt)
  logcolors = np.zeros(nfilt-1)

  lumdist = lumdistance_lcdm(zarray)
  # convert to Mpc if needed, refdist is 10 Mpc usually

  for i in range(nz):
      z = zarray[i]
      logdistfac = 2.0 * np.log10(lumdist[i]/refdist)
      fluxarray = np.zeros(nsed,nfilt)
      colorarray = np.zeros(nsed,nfilt-1)
      for k in range(nsed):
          for j in range(nfilt):
              # should be multiplying filter (?) by 1+z here for the frequency units?
              fluxes[j] = sedflux(sedstruct[k]['wave']*(1.0+z),
                                  sedstruct[k]['flux'],
                                  filtstruct[j]['wave'],
                                  filtstruct[j]['response'] * (1.0+z))
          logfluxes = np.log10(fluxes) - logdistfac
          fluxes = 10**logfluxes
          logcolors = logfluxes[0:nfilt-1] - logfluxes[1:nfilt]
          if return_log != 0:
            struct1 = {'z': z, 'logfluxes': logfluxes, 'logcolors': logcolors}
            fluxarray[k] = logfluxes
          else:
            struct1 = {'z': z, 'fluxes': fluxes, 'logcolors': logcolors}
            fluxarray[k] = fluxes
          colorarray[k] = logcolors
      #trackstruct.append(struct1)
      struct2 = {'z': z, 'fluxes': fluxarray, 'logcolors': colorarray}
      trackstruct.append(struct2)
      
  # this makes trackstruct a list of dictionaries, one dictionary for
  # each redshift

  # make plot?
  
  return trackstruct

