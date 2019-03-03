#
# use ureka to get newer scipy
# set import to read from python/sedfitting
# eg
#   ur_setup
#   PYTHONPATH=/Users/bjw/software/ureka/Ureka/python/lib/python2.7/site-packages/
#   export PYTHONPATH=$PYTHONPATH:$HOME/python/sedfitting
#   python makeplots1.py
# or
#   python ~/text/conf/cmu-stat-jun16/makeplots1.py

# see python/sedfitting/ Readme.testing, Readme.montecarlo

import numpy as np
import matplotlib.pyplot as plt

import scipy.special
from scipy import optimize

# maybe:   from sedfitting import ...

# Change or add chary in place of rieke - Feb 2019

import read_rieke_seds, read_ir_filters
import read_chary_seds
import sedflux, make_sedflux_z
import lumdistance_lcdm
import read_one_draineli_model, read_draineli_models
import convert_draineli_sed
import composite_draineli_sed
import fitonesed
import fitsedfamily, fitsedfamily_mc
import fitsedcombine, fitsedcombine_mc

import plot_data_sed

# upper LIR limit for LIR plots
# 13.2 to show everything, 12.7 to suppress the top 2 SEDs that George
# says are total extrapolations
# loglir_uplim = 13.2
loglir_uplim = 12.7

plotdir='chary_plots_v2/'

# read filters and seds

filt1 = read_ir_filters.read_ir_filters(makeplot=0)
filtwaves = np.zeros(len(filt1))
for i in range(len(filt1)):
    filtwaves[i] = filt1[i]['label']

# sedrieke = read_rieke_seds.read_rieke_seds(makeplot=0)
# rieke_loglir = 9.75 + 0.25 * np.arange(len(sedrieke))

# I made the Chary file to be like the Rieke file and loglir are the same
sedchary = read_chary_seds.read_chary_seds(makeplot=0)
chary_loglir = 9.75 + 0.25 * np.arange(len(sedchary))
nsed_max = np.size(np.where(chary_loglir < loglir_uplim))
# nsed_touse = len(sedchary)
nsed_touse = nsed_max

ztest = 0.003
# flux1 = make_sedflux_z.make_sedflux_z(sedrieke,filt1,ztest)
flux1 = make_sedflux_z.make_sedflux_z(sedchary,filt1,ztest)

# read DL models

fname = '/Users/bjw/dustmass/draine_li_2007/list.U1.00.model_subset1'
direc = '/Users/bjw/dustmass/draine_li_2007'
dlmodels_set1 = read_draineli_models.read_draineli_models(fname,dir=direc,makeplot=0)

dlseds_set1 = convert_draineli_sed.convert_draineli_sed(dlmodels_set1)

fname = '/Users/bjw/dustmass/draine_li_2007/list.models.largesubset1'
direc = '/Users/bjw/dustmass/draine_li_2007'
dlmodels_set2 = read_draineli_models.read_draineli_models(fname,dir=direc,makeplot=0)

dlseds_set2 = convert_draineli_sed.convert_draineli_sed(dlmodels_set2)

# make composite DL models

fname = '/Users/bjw/dustmass/draine_li_2007/list.composite_models1.umax'
direc = '/Users/bjw/dustmass/draine_li_2007'
dlmodels_part1 = read_draineli_models.read_draineli_models(fname,dir=direc,makeplot=0)
dlseds_part1 = convert_draineli_sed.convert_draineli_sed(dlmodels_part1)

fname = '/Users/bjw/dustmass/draine_li_2007/list.composite_models1.umin'
direc = '/Users/bjw/dustmass/draine_li_2007'
dlmodels_part2 = read_draineli_models.read_draineli_models(fname,dir=direc,makeplot=0)
dlseds_part2 = convert_draineli_sed.convert_draineli_sed(dlmodels_part2)

gammavals = [0.0, 0.1, 0.2, 0.3]
dlseds_composite = composite_draineli_sed.composite_draineli_sed(dlseds_part1, dlseds_part2, gammavals, makeplot=0)

import copy

# renormalize some DL models to be in units of 1e6 msun
# This makes fitting much more stable since coeffs are near 1.
# are there syntax problems here?
dlseds_composite_renorm = copy.deepcopy(dlseds_composite)
for i in range(len(dlseds_composite_renorm)):
  tmp1 = 1e6 * dlseds_composite_renorm[i]['flux']
  dlseds_composite_renorm[i]['flux'] = tmp1



####
# plot the Chary template spectra
# change to obey loglir_uplim

plt.clf()
for i in range(nsed_touse):
    linestyle = 'k-'
    plt.plot(np.log10(sedchary[i]['wave']), np.log10(sedchary[i]['flux']), linestyle)
    dotstyle = 'ko'
    plt.plot(np.log10(filtwaves),np.log10(flux1[i,0:]),dotstyle)

# fig = plt.xlim(0.3,3.0)
fig = plt.axis([0.5,3.0,-1.0,5.5])
fig = plt.xlabel('log wavelength, microns')
fig = plt.ylabel('Chary template flux, Jy')
plt.savefig(plotdir + 'chary_templ_logflux.pdf')

####
# plot some Draine & Li models

# dlseds_plot = dlseds_composite_best9
# dlseds_plot = dlseds_composite
dlseds_plot = dlseds_composite_renorm

plt.clf()
for i in range(len(dlseds_plot)):
    plt.plot(np.log10(dlseds_plot[i]['wave']),np.log10(dlseds_plot[i]['flux']),'k-')
#plt.xlim(0.3,3.0)
#plt.axis([0.3,3.0,-11.0,-4.0])
plt.axis([0.3,3.0,-5.0,2.0])
fig = plt.xlabel('log wavelength, microns')
fig = plt.ylabel('DL07 model flux for 10^6 Msun')
plt.savefig(plotdir + 'dlseds_composite_renorm_flux.pdf')

# moved the plotting of the best-9 models to later

####
# fit best single SED to a series of Chary models and plot each

# dlseds_fituse = dlseds_set2
dlseds_fituse = dlseds_composite_renorm

iobs_array = [1,3,5,7,9]
lir_name = ['10','10.5','11','11.5','12']

for ii in range(len(iobs_array)):
    # iobs = 5
    iobs = iobs_array[ii]
    plotname = 'chary_lir' + lir_name[ii] + '_onesed_fit.pdf'
    # iobs = 5
    testwave = filtwaves
    zobs = 0.003
    testflux = flux1[iobs,0:]
    testferr = 0.1 * testflux
    nbest, fnormarray, chisqarray = fitsedfamily.fitsedfamily(dlseds_fituse, zobs, testwave, testflux, testferr, filt1, makeplot=0, logplot=1)
    nobs=len(testwave)
    probarray = scipy.special.gammaincc((nobs-2)/2.0, chisqarray/2.0)
    totprob = sum(probarray)

    fitwave_model = dlseds_fituse[nbest]['wave'] * (1+zobs)
    fitflux_model = fnormarray[nbest] * dlseds_fituse[nbest]['flux']
    fpredict = np.zeros(len(testwave))
    for i in range(len(fpredict)):
        fpredict[i] = fnormarray[nbest] * sedflux.sedflux(fitwave_model, dlseds_fituse[nbest]['flux'], filt1[i]['wave'], filt1[i]['response'])
        fitwave = testwave

    plt.clf()
    plt.plot(np.log10(sedchary[iobs]['wave']), np.log10(sedchary[iobs]['flux']), 'k-')
    plt.plot(np.log10(testwave), np.log10(testflux), 'ko')
    plt.errorbar(np.log10(testwave), np.log10(testflux), yerr=testferr/testflux/2.3026, fmt='ko')
    plt.plot(np.log10(fitwave), np.log10(fpredict), 'rx')
    plt.plot(np.log10(fitwave_model), np.log10(fitflux_model), 'r-')
    # plt.xlim(0.3,3.0)
    ax = plt.axis([0.3,3.0,0.0,4.0])
    plt.xlabel('log wavelength')
    plt.ylabel('log flux, Chary template + 1-model fit')
    plt.savefig(plotdir + plotname)

# plt.savefig(plotdir + 'chary_lir11_onesed_fit.pdf')

#wobs2, fpredict2 = plot_data_sed.plot_data_sed(dlseds_composite, nbest, fnormarray[nbest], zobs, testwave, testflux, testferr, filt1, logplot=1)
#plt.savefig(plotdir + 'chary_lir11_onesed_fitv2.pdf')

##########
# monte carlo one model SED at a time to Chary templates

# dlseds_fituse = dlseds_set2
dlseds_fituse = dlseds_composite_renorm

ztest = 0.003

# nsed = len(sedchary)
nsed = nsed_touse
testwave = filtwaves
zobs = 0.003
# nmonte = 100
nmonte = 20
result_norm_list = []
result_prob_list = []
result_mc_list = []
for iobs in range(nsed):
   testflux = flux1[iobs,0:]
   testferr = 0.1 * testflux
   # testferr = 0.2 * testflux
   result_norm, result_prob, result_mcfits = fitsedfamily_mc.fitsedfamily_mc(dlseds_fituse, zobs, testwave, testflux, testferr, filt1, nmonte, makeplot=0)
   result_norm_list.append(result_norm)
   result_prob_list.append(result_prob)
   result_mc_list.append(result_mcfits)

print " best     rms      expect    rms     wtmean    median-err"
for i in range(nsed):
   tmp1 = np.array(result_norm_list[i])
   print '%7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f' % tuple(tmp1)

result_mc_list_10 = copy.deepcopy(result_mc_list)

# Look at result_mc_list to see what SEDs are used
nbest10all = []
for i in range(len(result_mc_list_10)):
   # nbest10all.append(result_mc_list_10[i]['nbest'])
   nbest10all = nbest10all + result_mc_list_10[i]['nbest']
#nbest20all = []
#for i in range(len(result_mc_list_20)):
#   # nbest20all.append(result_mc_list_20[i]['nbest'])
#   nbest20all = nbest20all + result_mc_list_20[i]['nbest']

# plot histogram of which spectra get used as best fits
bestmodels10, counts10 = np.unique(np.array(nbest10all),return_counts=True)
#bestmodels20, counts20 = np.unique(np.array(nbest20all),return_counts=True)
counts10sort = sorted(counts10, reverse=True)
#counts20sort = sorted(counts20, reverse=True)
ntotmodels = sum(counts10)
counts10sortfrac = np.array(counts10sort)/float(ntotmodels)
#counts20sortfrac = np.array(counts20sort)/float(ntotmodels)

# print the indexes of the N most often used models
model_index_counts = zip(bestmodels10, counts10)
tmp1 = sorted(model_index_counts, key=lambda elem: elem[1])
# this undoes the zip, making two tuples sorted by the counts
tmp2 = zip(*tmp1)
model_index_sorted = tmp2[0]
models_10mostfrequent = model_index_sorted[0:10]
print "indexes of most used models: ", models_10mostfrequent
print "count frac of most used models:  ", counts10sortfrac[0:10]

plt.clf()
plt.xlabel('DL07 SED models ordered by fit popularity')
plt.ylabel('Fraction of best fits that are model N')
ax1 = plt.step(range(len(counts10)),counts10sortfrac,'b-')
#ax2 = plt.step(range(len(counts20)),counts20sortfrac,'r-')
plt.text(7,0.25,'10% flux errors',color='blue')
#plt.text(7,0.2,'20% flux errors',color='red')
#plt.figlegend( (ax1[0], ax2[0]), ('10% flux errors', '20% flux errors'), 'upper right')
plt.savefig(plotdir + 'hist_modelcounts_fitonesed.pdf')

mass_result_best = np.zeros(nsed)
mass_result_bestrms = np.zeros(nsed)
mass_result_expect = np.zeros(nsed)
mass_result_expectrms = np.zeros(nsed)
mass_result_marginalrms = np.zeros(nsed)
for i in range(nsed):
    mass_result_best[i] = result_norm_list[i][0]
    mass_result_bestrms[i] = result_norm_list[i][1]
    mass_result_expect[i] = result_norm_list[i][2]
    mass_result_expectrms[i] = result_norm_list[i][3]
    mass_result_marginalrms[i] = result_norm_list[i][5]

#mass_result_best = result_norm_list[0:,0]
#mass_result_bestrms = result_norm_list[0:,1]
#mass_result_expect = result_norm_list[0:,2]
#mass_result_expectrms = result_norm_list[0:,3]
#mass_result_marginalrms = result_norm_list[0:,5]

# plot log Mdust estimate as fn of log LIR
plt.clf()
plt.plot(chary_loglir[0:nsed], np.log10(mass_result_expect)+6, 'ko')
plt.errorbar(chary_loglir[0:nsed], np.log10(mass_result_expect)+6, yerr=mass_result_expectrms/mass_result_expect/2.3026,fmt='ko')
fig = plt.xlabel('Chary log IR luminosity, Lsun')
fig = plt.ylabel('log dust mass, Msun')
fig = plt.axis([9.55,loglir_uplim,7.0,9.5])
plt.savefig(plotdir + 'loglir_logmdust_onesed_expect.pdf')

# plot ratio of Mdust error from MC to median error est from marginalizing
# over probablilities in single sim
# this may not be meaningful if the MC realizations mostly get stuck
# on the same SED

# suppress the log lir = 12.75 and 13.0 points because the SEDs are
# extrapolations and the fit failed

# itoplot = range(len(chary_loglir))
# ntoplot = len(chary_loglir) - 2
ntoplot = nsed_touse
itoplot = range(ntoplot)

plt.clf()
plt.plot(chary_loglir[itoplot], mass_result_expectrms[itoplot]/mass_result_marginalrms[itoplot], 'ko')
plt.plot(chary_loglir[itoplot], mass_result_expectrms[itoplot]/mass_result_marginalrms[itoplot], 'k-')
fig = plt.xlabel('Chary log IR luminosity, Lsun')
fig = plt.ylabel('Mdust error: MC RMS / marginal RMS')
# plt.xlim(9.5,13.25)
fig = plt.axis([9.55,loglir_uplim,0.0,10.0])
plt.savefig(plotdir + 'loglir_mdust_onesed_error_ratio.pdf')

##########
# most frequently used SEDs in some fits I did earlier, for Rieke
#
# These are old.
# indexbest3 = [ 56, 76, 149]
# indexbest9 = [ 56, 76, 113, 116, 133, 136, 149, 181, 201]
# Use the indexes from the sorted list of most popular
indexbest3 = model_index_sorted[0:3]
indexbest9 = model_index_sorted[0:9]

# There's probably a better way but this works
dlseds_composite_best9 = []
for i in indexbest9:
   dlseds_composite_best9.append(dlseds_composite[i])
dlseds_composite_best3 = []
for i in indexbest3:
   dlseds_composite_best3.append(dlseds_composite[i])

dlseds_composite_renorm9 = copy.deepcopy(dlseds_composite_best9)
for i in range(len(dlseds_composite_renorm9)):
  tmp1 = 1e6 * dlseds_composite_renorm9[i]['flux']
  dlseds_composite_renorm9[i]['flux'] = tmp1

dlseds_composite_renorm3 = copy.deepcopy(dlseds_composite_best3)
for i in range(len(dlseds_composite_renorm3)):
  tmp1 = 1e6 * dlseds_composite_renorm3[i]['flux']
  dlseds_composite_renorm3[i]['flux'] = tmp1

dlseds_plot = dlseds_composite_renorm9

plt.clf()
for i in range(len(dlseds_plot)):
    plt.plot(np.log10(dlseds_plot[i]['wave']),np.log10(dlseds_plot[i]['flux']),'k-')
#plt.xlim(0.3,3.0)
#plt.axis([0.3,3.0,-11.0,-4.0])
plt.axis([0.3,3.0,-5.0,2.0])
fig = plt.xlabel('log wavelength, microns')
fig = plt.ylabel('DL07 model flux for 10^6 Msun')
plt.savefig(plotdir + 'dlseds_composite_renorm9_flux.pdf')

##########
# fit a combination to a single Chary SED and plot

dlseds_fituse = dlseds_composite_renorm9

iobs_array = [1,3,5,7,9]
lir_name = ['10','10.5','11','11.5','12']

for ii in range(len(iobs_array)):
    # iobs = 5
    iobs = iobs_array[ii]
    plotname = 'chary_lir' + lir_name[ii] + '_combine_fit.pdf'
    testwave = filtwaves
    zobs = 0.003
    testflux = flux1[iobs,0:]
    testferr = 0.1 * testflux
    fitcoeffs, fiterrors, chisq = fitsedcombine.fitsedcombine(dlseds_fituse, zobs, testwave, testflux, testferr, filt1, penalize=1.0, initguess=0.0, makeplot=0, logplot=1)
    nobs=len(testwave)

    plt.clf()
    plt.plot(np.log10(sedchary[iobs]['wave']), np.log10(sedchary[iobs]['flux']), 'k-')
    plt.plot(np.log10(testwave), np.log10(testflux), 'ko')
    plt.errorbar(np.log10(testwave), np.log10(testflux), yerr=testferr/testflux/2.3026, fmt='ko')
    fpredict = np.zeros(len(testwave))
    fsum = np.zeros(len(testwave))
    fsum_model = np.zeros(len(dlseds_fituse[0]['wave']))
    for j in range(len(dlseds_fituse)):
        if fitcoeffs[j] > 1.0e-6:
            fitwave_model = dlseds_fituse[j]['wave']
            fitflux_model = fitcoeffs[j] * dlseds_fituse[j]['flux']
            for k in range(nobs):
                fpredict[k] = fitcoeffs[j] * sedflux.sedflux(fitwave_model, dlseds_fituse[j]['flux'], filt1[k]['wave'], filt1[k]['response'])
                fitwave=testwave
                plt.plot(np.log10(fitwave), np.log10(fpredict), 'bx')
                plt.plot(np.log10(fitwave_model), np.log10(fitflux_model), 'b-')
            fsum = fsum + fpredict
            fsum_model = fsum_model + fitflux_model

    plt.plot(np.log10(fitwave), np.log10(fsum), 'rx')
    plt.plot(np.log10(fitwave_model), np.log10(fsum_model), 'r-')
    #plt.xlim(0.3,3.0)
    ax = plt.axis([0.3,3.0,0.0,4.0])
    plt.xlabel('log wavelength')
    plt.ylabel('log flux, Chary template + combined fit')
    plt.savefig(plotdir + plotname)

# plt.savefig(plotdir + 'chary_lir11.5_combine_fit.pdf')

#stop

####
# monte carlo of fitting combination over renorm best9 modes to all
# Chary templates

# try penalize=1.0, initguess=0, and SLSQP

dlseds_fituse = dlseds_composite_renorm9

#nsed = len(sedchary)
nsed = nsed_touse
testwave = filtwaves
zobs = 0.003
# nmonte = 100
# nmonte = 20
nmonte = 40
result_coeffs_list = []
result_prob_list = []
result_mc_list = []
for iobs in range(nsed):
   testflux = flux1[iobs,0:]
   testferr = 0.1 * testflux
   # testferr = 0.2 * testflux
   result_coeffs, result_prob, result_mcfits = fitsedcombine_mc.fitsedcombine_mc(dlseds_fituse, zobs, testwave, testflux, testferr, filt1, nmonte, penalize=1.0, initguess=0, makeplot=0, logplot=0)
   result_coeffs_list.append(result_coeffs)
   result_prob_list.append(result_prob)
   result_mc_list.append(result_mcfits)

result_mc_list_10 = result_mc_list[:]

for i in range(nsed):
   tmp1 = (result_coeffs_list[i]['meansum'], result_coeffs_list[i]['rmssum'],
           result_coeffs_list[i]['meannzero'], result_coeffs_list[i]['rmsnzero'])
   print '%6.2f  %5.2f  %5.2f  %5.2f' % tmp1

mass_combine_mean = np.zeros(nsed)
mass_combine_rms = np.zeros(nsed)
for i in range(nsed):
    mass_combine_mean[i] = result_coeffs_list[i]['meansum']
    mass_combine_rms[i] = result_coeffs_list[i]['rmssum']
    
# plot log Mdust estimate as fn of log LIR
plt.clf()
plt.plot(chary_loglir[0:nsed], np.log10(mass_combine_mean)+6, 'ko')
plt.errorbar(chary_loglir[0:nsed], np.log10(mass_combine_mean)+6, yerr=mass_combine_rms/mass_combine_mean/2.3026,fmt='ko')
fig = plt.xlabel('log IR luminosity, Lsun')
fig = plt.ylabel('log combined dust mass, Msun')
fig = plt.axis([9.55,loglir_uplim,7.0,9.5])
plt.savefig(plotdir + 'loglir_logmdust_combine_mean.pdf')

# plot fit to one template showing all fitted components - see above


####

#stop


####
# plot comparing the log LIR estimates from combine and onesed

plt.clf()
plt.subplot(1,1,1)
plt.plot(chary_loglir[itoplot], np.log10(mass_result_expect[itoplot])+6, 'ro')
plt.plot(chary_loglir[itoplot], np.log10(mass_result_expect[itoplot])+6, 'r-')
plt.errorbar(chary_loglir[itoplot], np.log10(mass_result_expect[itoplot])+6, yerr=mass_result_expectrms[itoplot]/mass_result_expect[itoplot]/2.3026,fmt='ro')
fig = plt.xlabel('Chary log IR luminosity, Lsun')
fig = plt.ylabel('log dust mass, Msun')
fig = plt.text(9.8,9.1,'one SED fit',color='red')
#fig = plt.axis([9.55,loglir_uplim,7.0,9.5])
#plt.subplot(2,1,2)
plt.plot(chary_loglir[itoplot], np.log10(mass_combine_mean[itoplot])+6, 'bo')
plt.plot(chary_loglir[itoplot], np.log10(mass_combine_mean[itoplot])+6, 'b-')
plt.errorbar(chary_loglir[itoplot], np.log10(mass_combine_mean[itoplot])+6, yerr=mass_combine_rms[itoplot]/mass_combine_mean[itoplot]/2.3026,fmt='bo')
#fig = plt.xlabel('Chary log IR luminosity, Lsun')
#fig = plt.ylabel('log combined dust mass, Msun')
fig = plt.text(9.8,8.7,'combined SED fit',color='blue')
fig = plt.axis([9.55,loglir_uplim,7.0,9.5])
plt.savefig(plotdir + 'loglir_logmdust_combine_and_onesed.pdf')

logmassdiff = np.log10(mass_result_expect[itoplot]) - np.log10(mass_combine_mean[itoplot])

plt.subplot(1,1,1)

####
# plot comparing the error estimates

plt.clf()
plt.subplot(1,1,1)
#logerror_ratio = (mass_combine_rms/mass_combine_mean/2.3026) / (mass_result_expectrms/mass_result_expect/2.3026)
error_ratio = mass_result_expectrms / mass_combine_rms
error_ratio_dex = np.log10(error_ratio)
plt.plot(chary_loglir[itoplot], error_ratio, 'ko')
plt.plot(chary_loglir[itoplot], error_ratio, 'k-')
fig = plt.xlabel('Chary log IR luminosity, Lsun')
fig = plt.ylabel('error estimate ratio, 1 SED / combination')
fig = plt.axis([9.55,loglir_uplim,0.0,8.0])
plt.savefig(plotdir + 'loglir_errmdust_ratio.pdf')

#####

print "Summary:"
print "number of DL models in U1 set and in large subset: ",len(dlmodels_set1),len(dlmodels_set2)
print "number of DL models, composite: ",len(dlseds_composite)
print "number of DL models plotted: ",len(dlseds_composite_renorm)," and",len(dlseds_plot)
print "number of DL models used in 1-model fits: ",len(dlseds_fituse)
print "indexes of most used models: ", models_10mostfrequent
print "count frac of most used models:  ", counts10sortfrac[0:10]
print "indexes of models I was using in fit: ", indexbest9
print "number of galaxy templates fitted: ", nsed
print "log mass offset, onesed - combine: ", logmassdiff
print "mean log mass offset over the SEDs: ", np.mean(logmassdiff)
