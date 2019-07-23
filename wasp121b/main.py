import os, sys

import numpy as np
import scipy
import scipy.interpolate
from scipy.interpolate import UnivariateSpline
from scipy.signal import lombscargle

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy import constants as c

# allesfitter
import allesfitter
from allesfitter.exoworlds_rdx.lightcurves.index_transits import index_eclipses
from allesfitter.exoworlds_rdx.lightcurves.lightcurve_tools import rebin_err, phase_fold

# own modules
import tdpy.mcmc
from tdpy.util import summgene

# Max's modules
#from exoworlds.lightcurves.lightcurve_tools import phase_fold
from exoworlds.tess import extract_SPOC_data, extract_QLP_data

import emcee
import dynesty

def retr_lpos(para, meanphas, data, datastdv, modltype, indxoccl):
    
    offs = para[0]
    amfo = para[1:8]
    tempdayy = para[8]
    tempnigh = para[9]
    if modltype == 'shft':
        shft = para[10]
    else:
        shft = 0.
    
    if ((para < limtpara[0, :]) | (para > limtpara[1, :])).any():
        lpos = -np.inf
    else:
        llik, deptatmo = retr_llikfull(meanphas, data, datastdv, offs, amfo, tempdayy, tempnigh, shft, modltype, indxoccl)
        lpri = - 0.5 * ((amfo[3] / 0.1)**2)
        lpos = llik + lpri
        
    return lpos


def retr_modl(meanphas, offs, amfo, tempdayy, tempnigh, shft, modltype, indxoccl):
    
    phasshft = (meanphas - shft * np.pi / 180.) * 2. * np.pi
    occl = retr_occl(tempdayy, tempstar)
    deptnigh = retr_occl(tempnigh, tempstar)
    deptatmo = occl - deptnigh
    modl = offs + 0.5 * 1e-6 * deptatmo * np.cos(phasshft)
    for k in range(3):
        modl += 0.5e-6 * amfo[k] * np.sin((k + 1) * meanphas * 2. * np.pi) 
    for k in range(2):
        modl += 0.5e-6 * amfo[k+5] * np.cos((k + 2) * meanphas * 2. * np.pi)
    
    modl[indxoccl] -= occl * 1e-6
    
    return modl, deptatmo


def retr_occl(tempdayy, tempstar):

    specplan = thpt / meanlamb**5 / (np.exp(cons / meanlamb / tempdayy) - 1.)
    specstar = thpt / meanlamb**5 / (np.exp(cons / meanlamb / tempstar) - 1.)
    fluxplan = np.sum(difflamb * specplan)
    fluxstar = np.sum(difflamb * specstar)
    occl = fluxplan / fluxstar
    
    return occl


def retr_llik(meanphas, data, datastdv, offs, amfo, tempdayy, tempstar, shft, modltype, indxoccl):
    
    llik, deptatmo = retr_llikfull(meanphas, data, datastdv, offs, amfo, tempdayy, tempstar, shft, modltype, indxoccl)

    return llik


def retr_llikfull(meanphas, data, datastdv, offs, amfo, tempdayy, tempstar, shft, modltype, indxoccl):
    
    modl, deptatmo = retr_modl(meanphas, offs, amfo, tempdayy, tempstar, shft, modltype, indxoccl)
    llik = -0.5 * np.sum((modl - data)**2 / datastdv**2)
    
    return llik, deptatmo


def icdf(para):
    
    icdf = limtpara[0, :] + para * (limtpara[1, :] - limtpara[0, :])

    return icdf


def retr_mask(time, flag, strgmask):
    
    frst = time[0]
    four = time[-1]
    mean = np.mean(time)
    diff = time - mean
    seco = diff[np.where(diff < 0)][-1]
    thrd = diff[np.where(diff > 0)][0]
    
    indx = np.where(~((time < frst + 0.5) | ((time < thrd + 0.5) & (time > mean)) | (time > four - 0.25) | \
                                                                    ((time < 2458511.8) & (time > 2458511.3))))[0]
    
    if strgmask == 'ful1':
        indx = np.setdiff1d(indx, np.where(time > 2458504)[0])
    if strgmask == 'ful2':
        indx = np.setdiff1d(indx, np.where(time < 2458504)[0])
    
    return indx

# paths
pathbase = '/Users/tdaylan/DropboxMITTTT/wasp121b/tansu/'
pathstor = pathbase + 'shft/stor/'
pathdata = '/Users/tdaylan/Desktop/temp/WASP-121/data/'
pathimag = pathstor + 'imag/'

# parameters
# WASP-121
## period
peri = 1.2749255 # [day]
## epoch
epoc = 2456635.70832 # [BJD]
## Half of the transit duration
dura = 0.1203 / 2. # [day] Delrez2016

# extract SPOC data
extract_SPOC_data([pathdata + 'data_preparation/archival_data/tess2019006130736-s0007-0000000022529346-0131-s_lc.fits'], 
                  outdir=pathdata+'data_preparation/TESS_SAP', PDC=False, auto_correct_dil=True, extract_centd=True, extract_dil=True)

extract_SPOC_data([pathdata + 'data_preparation/archival_data/tess2019006130736-s0007-0000000022529346-0131-s_lc.fits'], 
                  outdir=pathdata+'data_preparation/TESS_PDCSAP', PDC=True, auto_correct_dil=True, extract_centd=True, extract_dil=True)

# get TESS throughput
data = np.loadtxt(pathdata + 'band.csv', delimiter=',', skiprows=9)
meanlambband = data[:, 0] * 1e-3
thptband = data[:, 1]

# plot PDCSAP and SAP light curves
time, lcurpdcc, stdvlcurpdcc = np.genfromtxt(pathdata + 'data_preparation/TESS_PDCSAP/TESS.csv', delimiter=',', unpack=True)
time, lcursapp, stdvlcursapp = np.genfromtxt(pathdata + 'data_preparation/TESS_SAP/TESS.csv', delimiter=',', unpack=True)

figr, axis = plt.subplots(2, 1, figsize=(12, 5))
axis[0].plot(time - 2457000, lcurpdcc, 'k.')
axis[1].plot(time - 2457000, lcursapp, '.', color='black')
axis[0].set_xticklabels([])
axis[1].set_xlabel('Time [BJD - 2457000]')
axis[0].set_ylabel('Relative Flux')
axis[1].set_ylabel('Relative Flux')
path = pathimag + 'lcur.pdf'
plt.subplots_adjust(hspace=0.)
plt.savefig(path)
plt.close()

# fit a spline to the SAP light curve
timebrek = [0., 2458504., np.inf]
numbbrek = len(timebrek)
lcurdetr = []
lcurspln = []
print 'Fitting the spline...'
for i in np.arange(numbbrek - 1):
    ind = np.where((time>=timebrek[i]) & (time<=timebrek[i+1]))[0]
    time1 = time[ind]
    flux1 = lcursapp[ind]
    flux_err1 = stdvlcursapp[ind]
    ind_ecl1, ind_ecl2, ind_out = index_eclipses(time1, epoc, peri, width_1=4./24., width_2=4./24.)
    spl = UnivariateSpline(time1[ind_out], flux1[ind_out]-1.)
    print 'i'
    print i
    print '$\beta$:', spl.get_coeffs()
    print '$t_k$:', spl.get_knots()
    xx = np.linspace(time1[0], time1[-1], 1000)
    axis[0].plot(time1 - 2457000, flux1, 'k.')
    axis[0].plot(time1[ind_out] - 2457000, flux1[ind_out], ls='', color='black')
    axis[0].plot(xx - 2457000, spl(xx)+1., 'm-' )
    axis[1].plot(time1 - 2457000, flux1-spl(time1), '.', color='black')
    lcurdetr += list(flux1-spl(time1))
    lcurspln += list(spl(time1))
lcurdetr = np.array(lcurdetr)
lcurspln = np.array(lcurspln)

figr, axis = plt.subplots(2, 1, figsize=(12, 8))
axis[0].set_ylabel('Relative Flux')
axis[1].set_ylabel('Relative Flux')
axis[0].set_xticklabels([])
axis[1].set(xlabel='Time [BJD - 2457000]')
plt.subplots_adjust(hspace=0.)
path = pathimag + 'lcur_spln.pdf'
plt.savefig(path)
plt.close()
    
ind_ecl1, ind_ecl2, ind_out = index_eclipses(time, epoc, peri, width_1=4./24., width_2=2./24.)
offset = np.mean(lcurdetr[ind_ecl2])-1.
lcurdetr -= offset

# plot the phase-folded splined light curve
## phase-fold the splined light curve
phasfold, pcur, stdvpcur, N, phas = phase_fold(time, lcurdetr, peri, epoc, dt = 0.005, ferr_type='meansig', ferr_style='sem', sigmaclip=True)
#phasfold, pcur, stdvpcur, N, phas = phase_fold(time, lcurdetr, peri, epoc, dt=0.001, ferr_type='meansig', ferr_style='sem', sigmaclip=True)
#phasfold, pcur, stdvpcur, N, phas = phase_fold(time, lcurdetr, peri, epoc, dt=0.001/5., ferr_type='meansig', ferr_style='sem', sigmaclip=True)
figr, axis = plt.subplots(2, 1, figsize=(12, 8))
axis[0].plot(phas, lcurdetr, '.', color='grey', alpha=0.3)
axis[0].plot(phasfold, pcur, 'k.')
axis[1].plot(phas, lcurdetr, '.', color='grey', alpha=0.3)
axis[1].plot(phasfold, (pcur - 1) * 1e6, 'k.', ms=10)
axis[0].set(ylabel='Relative Flux')
axis[1].set(ylabel='Relative Flux - 1 [ppm]')
axis[1].set(ylim=[-800, 1000])
axis[0].set_xticklabels([])
axis[1].set_xlabel('Phase')
plt.subplots_adjust(hspace=0.)
path = pathimag + 'pcur.pdf'
plt.savefig(path)
plt.close()

data = np.column_stack((time, lcurdetr, stdvlcursapp))
np.savetxt(pathdata + 'data_preparation/.csv', data, delimiter=',', header='time,flux,flux_err')

figr, axis = plt.subplots(figsize=(12, 5))
axis.plot(phas, lcurdetr, '.', color='grey', alpha=0.3)
axis.errorbar(phasfold, pcur, yerr=stdvpcur, ls='', fmt='black', capsize=0, ms=8)
axis.set(xlabel='Phase', ylabel='Flux')
path = pathimag + 'pcur.pdf'
plt.savefig(path)
plt.close()

# read the allesfitter posterior
companion = 'b'
inst = 'TESS'
key = 'flux'
pathalle = pathdata + 'allesfits/allesfit_full_gp/'
alles = allesfitter.allesclass(pathalle)

radirati = alles.posterior_params['b_rr']
smajaxis = 0.02544 * 2092.51 # [RJ]
radiplan = 1.865 # [RJ]
radistar = 14.51 # [RJ]
rsma = alles.posterior_params['b_rsuma']
albe = alles.posterior_params['b_geom_albedo_TESS']

fluxrati = albe * (radirati * rsma / (1. + radirati))**2
sbrgrati = alles.posterior_params['b_sbratio_TESS'] * radirati**2
numbsamp = radirati.size
listsamp = np.empty((numbsamp, 3))
listsamp[:, 0] = sbrgrati * 1e6
listsamp[:, 1] = (sbrgrati + fluxrati) * 1e6
listsamp[:, 2] = fluxrati * 1e6

listlabl = ['Nightside [ppm]', 'Secondary [ppm]', 'Modulation [ppm]']
tdpy.mcmc.plot_grid(pathimag, 'post', listsamp, listlabl)

# plot a lightcurve from the posteriors
time = np.concatenate((alles.data['TESS']['time'], alles.data['TESS']['time']+alles.posterior_params_median['b_period']))
flux = np.concatenate((alles.data['TESS']['flux'], alles.data['TESS']['flux']))
model_time = 1.*time
model_flux = alles.get_posterior_median_model(inst, key, xx=model_time)
baseline = alles.get_posterior_median_baseline(inst, key, xx=time)
bintime, binflux, binflux_err, N = rebin_err(time, flux-baseline, ferr=None, dt = 0.02, phasefolded=False, ferr_type='meansig', ferr_style='sem', sigmaclip=True)

figr, axis = plt.subplots(3, 1, figsize=(12, 14))
for k in range(len(axis)):
    if k == 0:
        varbtemp = flux-baseline
    else:
        varbtemp = (flux-baseline - 1.) * 1e6
    axis[k].plot(time, varbtemp, '.', color='grey', alpha=0.3, label='Raw data')
    if k == 0:
        varbtemp = binflux
        varbtemperrr = binflux_err
    else:
        varbtemp = (binflux - 1.) * 1e6
        varbtemperrr = binflux_err * 1e6
    axis[k].errorbar(bintime, varbtemp, marker='o', yerr=varbtemperrr, linestyle='none', capsize=0, ls='', color='black', label='Binned data')
    if k == 0:
        varbtemp = model_flux
    else:
        varbtemp = (model_flux - 1.) * 1e6
    axis[k].plot(model_time, varbtemp, color='b', lw=2, label='Total Model')
    axis[k].set(xlabel='Time (BJD)')
axis[0].set(ylabel='Relative Flux')
axis[1].set(ylabel='Relative Flux - 1 [ppm]')
axis[2].set(ylabel='Relative Flux - 1 [ppm]')
axis[1].set(ylim=[-800,1000])
axis[2].set(xlim=[0.2*peri,0.8*peri], ylim=[-800, 1000])

alles = allesfitter.allesclass(pathalle)
alles.posterior_params_median['b_sbratio_TESS'] = 0
alles.settings['host_shape_TESS'] = 'sphere'
alles.settings['b_shape_TESS'] = 'sphere'
alles.posterior_params_median['b_geom_albedo_TESS'] = 0
alles.posterior_params_median['host_gdc_TESS'] = 0
alles.posterior_params_median['host_bfac_TESS'] = 0
model_flux2 = alles.get_posterior_median_model(inst, key, xx=model_time)
for k in range(len(axis)):
    if k == 0:
        varbtemp = model_flux2
    else:
        varbtemp = (model_flux2 - 1.) * 1e6
    axis[k].plot(model_time, varbtemp, lw=2, color='m', label='PT', ls='--')

alles = allesfitter.allesclass(pathalle)
#alles.posterior_params_median['b_sbratio_TESS'] = 0
alles.settings['host_shape_TESS'] = 'sphere'
alles.settings['b_shape_TESS'] = 'sphere'
alles.posterior_params_median['b_geom_albedo_TESS'] = 0
alles.posterior_params_median['host_gdc_TESS'] = 0
alles.posterior_params_median['host_bfac_TESS'] = 0
model_flux2 = alles.get_posterior_median_model(inst, key, xx=model_time)
for k in range(len(axis)):
    if k == 0:
        varbtemp = model_flux2
    else:
        varbtemp = (model_flux2 - 1.) * 1e6
    axis[k].plot(model_time, varbtemp, lw=2, color='g', label='PIL', ls='--', zorder=11)

alles = allesfitter.allesclass(pathalle)
alles.posterior_params_median['b_sbratio_TESS'] = 0
alles.settings['host_shape_TESS'] = 'sphere'
alles.settings['b_shape_TESS'] = 'sphere'
alles.posterior_params_median['b_geom_albedo_TESS'] = 0
#alles.posterior_params_median['host_gdc_TESS'] = 0
alles.posterior_params_median['host_bfac_TESS'] = 0
model_flux2 = alles.get_posterior_median_model(inst, key, xx=model_time)
for k in range(len(axis)):
    if k == 0:
        varbtemp = model_flux2
    else:
        varbtemp = (model_flux2 - 1.) * 1e6
    axis[k].plot(model_time, varbtemp, lw=2, color='orange', ls='--', label='SGD', zorder=11)

alles = allesfitter.allesclass(pathalle)
alles.posterior_params_median['b_sbratio_TESS'] = 0
#alles.settings['host_shape_TESS'] = 'sphere'
#alles.settings['b_shape_TESS'] = 'sphere'
alles.posterior_params_median['b_geom_albedo_TESS'] = 0
alles.posterior_params_median['host_gdc_TESS'] = 0
alles.posterior_params_median['host_bfac_TESS'] = 0
model_flux2 = alles.get_posterior_median_model(inst, key, xx=model_time)
for k in range(len(axis)):
    if k == 0:
        varbtemp = model_flux2
    else:
        varbtemp = (model_flux2 - 1.) * 1e6
    axis[k].plot(model_time, varbtemp, lw=2, color='r', ls='--', label='SEV', zorder=11)

alles = allesfitter.allesclass(pathalle)
alles.posterior_params_median['b_sbratio_TESS'] = 0
alles.settings['host_shape_TESS'] = 'sphere'
alles.settings['b_shape_TESS'] = 'sphere'
#alles.posterior_params_median['b_geom_albedo_TESS'] = 0
alles.posterior_params_median['host_gdc_TESS'] = 0
alles.posterior_params_median['host_bfac_TESS'] = 0
model_flux2 = alles.get_posterior_median_model(inst, key, xx=model_time)
for k in range(len(axis)):
    if k == 0:
        varbtemp = model_flux2
    else:
        varbtemp = (model_flux2 - 1.) * 1e6
    axis[k].plot(model_time, varbtemp, lw=2, color='m', ls='--', label='PRL', zorder=11)

axis[0].legend()
plt.subplots_adjust(hspace=0.)

path = pathimag + 'datamodl.pdf'
plt.savefig(path)
plt.close()


# plot the spherical limits
figr, axis = plt.subplots(2, 1, figsize=(12, 10))

alles.settings['host_shape_TESS'] = 'roche'
alles.settings['b_shape_TESS'] = 'roche'
model_fluxfull = alles.get_posterior_median_model(inst, key, xx=model_time)
axis[0].plot(model_time, model_fluxfull, 'r-', lw=2)
axis[0].set_xticklabels([])
axis[0].set_ylabel('Relative Flux')

alles.settings['host_shape_TESS'] = 'sphere'
alles.settings['b_shape_TESS'] = 'sphere'
model_flux = alles.get_posterior_median_model(inst, key, xx=model_time)
axis[1].plot(model_time, (model_flux - model_fluxfull) * 1e6, lw=2, label='Spherical star and planet residual')

alles.settings['host_shape_TESS'] = 'roche'
alles.settings['b_shape_TESS'] = 'sphere'
model_flux = alles.get_posterior_median_model(inst, key, xx=model_time)
axis[1].plot(model_time, (model_flux - model_fluxfull) * 1e6, lw=2, label='Spherical planet residual')
axis[1].legend()
axis[1].set(ylabel='Relative flux [ppm]')

plt.subplots_adjust(hspace=0.)

path = pathimag + 'pcurmodldiff.pdf'
plt.savefig(path)
plt.close()

# calculate prior on the mass ratio (Stassun+2017)
Mp = np.random.normal(loc=(375.99289 *c.M_earth/c.M_sun).value, scale=(20.34112*c.M_earth/c.M_sun).value, size=10000)
Ms = np.random.normal(loc=1.52644, scale=0.361148, size=10000)
q = Mp / Ms
print( 'q', np.mean(q), np.std(q), np.percentile(q, [16,50,84] ) )
print( 'q', np.percentile(q,50), np.percentile(q,50)-np.percentile(q,16), np.percentile(q,84)-np.percentile(q,50) )
fig = plt.figure()
plt.hist(Mp)
fig = plt.figure()
plt.hist(Ms)
fig = plt.figure()
plt.hist(q)





binslamb = np.linspace(0.6, 1., 1001)
meanlamb = (binslamb[1:] + binslamb[:-1]) / 2.
difflamb = (binslamb[1:] - binslamb[:-1]) / 2.
cons = 0.0143877735e6 # [um K]

## read eclipse data
liststrgfile = ['ContribFuncArr.txt', \
                'EmissionDataArray.txt', \
                #'RetrievalParamSamples.txt', \
                'ContribFuncWav.txt', \
                'EmissionModelArray.txt', \
                'RetrievalPTSamples.txt', \
                'pdependent_abundances/', \
                ]


# grid plot, along with C/O
#path = pathdata + 'ascii_output/RetrievalParamSamples.txt'
#listsamp = np.loadtxt(path)
#print 'listsamp'
#summgene(listsamp)
#tdpy.mcmc.plot_grid(pathimag, 'post', listsamp, ['$\kappa_IR$', '$\gamma$', '$\psi$', '[M/H]', '[C/H]', '[O/H]'])

# plot spectrum, depth, brightness temp
path = pathdata + 'ascii_output/ContribFuncWav.txt'
wlen = np.loadtxt(path)
# emission
# Columns are:
#   1. Wavelength (micron)
#   2. Best-fit retrieval eclipse depths (ppm)
#   3. Best-fit blackbody eclipse depths (ppm)
#   4. Blackbody envelope upper eclipse depths (ppm)
#   5. Blackbody envelope lower eclipse depths (ppm)
#   6. Best-fit retrieval planet surface flux (erg/s/cm^2)
#   7. Best-fit blackbody planet surface flux (erg/s/cm^2)
#   8. Blackbody envelope upper planet surface flux (erg/s/cm^2)
#   9. Blackbody envelope lower planet surface flux (erg/s/cm^2)
#  10. Stellar surface flux (erg/s/cm^2)
# data
# Columns are:
#  1. Central wavelength (micron)
#  2. Channel half-width (micron)
#  3. Eclipse depth value (ppm)
#  4. Eclipse depth uncertainty (ppm)
#  5. Brightness temperature value (Kelvin)
#  6. Brightness temperature uncertainty (Kelvin)
#  7. Planet surface flux value (erg/s/cm^2)
#  8. Planet surface flux uncertainty (erg/s/cm^2)
# Rows are:
#  1. TESS
#  2. Zprime (ground photometry)
#  3. Ks (ground photometry)
#  4. IRAC 3.6 micron
#  5. IRAC 4.5 micron
#  6-22. WFC3 G102
#  23-50. WFC3 G141
figr, axis = plt.subplots(4, 1, figsize=(12, 14))
path = pathdata + 'ascii_output/EmissionModelArray.txt'
arrymodl = np.loadtxt(path)
print 'arrymodl'
summgene(arrymodl)
path = pathdata + 'ascii_output/EmissionDataArray.txt'
arrydata = np.loadtxt(path)

data = np.loadtxt(pathdata + 'band.csv', delimiter=',', skiprows=9)
meanlambband = data[:, 0] * 1e-3 
thptband = data[:, 1] 

## stellar spectrum and TESS throughput
axis[0].plot(arrymodl[:, 0], 1e-9 * arrymodl[:, 9], label='Star', color='k')
axis[0].set_ylabel(r'$F_{\nu}$ [10$^9$ erg/s/cm$^2$]')
axistwin = axis[0].twinx()
axistwin.plot(meanlambband, thptband, color='grey')
axistwin.set_ylabel(r'TESS Throughput')

## eclipse depths
### model
axis[1].plot(arrymodl[:, 0], arrymodl[:,1], label='Retrieval')
axis[1].plot(arrymodl[:, 0], arrymodl[:,2], label='Blackbody', color='orange')
axis[1].fill_between(arrymodl[:, 0], arrymodl[:, 3], arrymodl[:, 4], alpha=0.1, color='orange')
### data
for k in range(8):
    labl = None
    axis[1].errorbar(arrydata[k, 0], arrydata[k, 2], xerr=arrydata[k, 1], yerr=arrydata[k, 3], label=labl, color=None, ls='', marker='o')
axis[1].set_ylabel(r'Depth [ppm]')
axis[1].set_xticklabels([])

## spectra
### model
axis[2].plot(arrymodl[:, 0], 1e-9 * arrymodl[:, 5], label='Retrieval', color='b')
axis[2].plot(arrymodl[:, 0], 1e-9 * arrymodl[:, 6], label='Blackbody', color='skyblue', alpha=0.3)
axis[2].fill_between(arrymodl[:, 0], 1e-9 * arrymodl[:, 7], 1e-9 * arrymodl[:, 8], color='skyblue', alpha=0.3)
### data
for k in range(8):
    axis[2].errorbar(arrydata[k, 0],  1e-9 * arrydata[k, 2], xerr=arrydata[k, 1], yerr=1e-9*arrydata[k, 3], \
                                                                                        label=None, color=None, ls='', marker='o')
axis[2].set_ylabel(r'$F_{\nu}$ [10$^9$ erg/s/cm$^2$]')
axis[2].set_xticklabels([])

## brightness temperature
### data
for k in range(8):
    if k == 0:
        labl = 'TESS (This work)'
    if k == 1:
        labl = 'Z$^\prime$'
    if k == 2:
        labl = '$K_s$'
    if k == 3:
        labl = 'IRAC 3.6 $\mu$m'
    if k == 4:
        labl = 'IRAC 4.5 $\mu$m'
    if k == 5:
        labl = 'HST/WFC3 G102'
    if k == 6:
        labl = 'HST/WFC3 G141'
    else:
        labl = None
#  1. TESS
#  2. Zprime (ground photometry)
#  3. Ks (ground photometry)
#  4. IRAC 3.6 micron
#  5. IRAC 4.5 micron
#  6-22. WFC3 G102
#  23-50. WFC3 G141
    axis[3].errorbar(arrydata[k, 0], arrydata[k, 4], xerr=arrydata[k, 1], yerr=arrydata[k, 5], label=labl, color=None, ls='', marker='o')
axis[3].set_ylabel(r'$T_B$ [K]')
axis[3].set_xlabel(r'$\lambda$ [$\mu$m]')
axis[3].legend(fancybox=True, bbox_to_anchor=[0.8, 0.8, 0.2, 0.2])
plt.subplots_adjust(hspace=0., wspace=0.)
path = pathimag + 'spec.pdf'
print 'Writing to %s...' % path
plt.savefig(path)
plt.close()


# plot pressure-temperature, contribution
path = pathdata + 'ascii_output/RetrievalPTSamples.txt'
dataptem = np.loadtxt(path)

path = pathdata + 'ascii_output/ContribFuncArr.txt'
dataconf = np.loadtxt(path)

path = pathdata + 'ascii_output/ContribFuncArr.txt'
data = np.loadtxt(path)
liststrgcomp = ['CH4.txt', 'CO.txt', 'FeH.txt', 'H+.txt', 'H.txt', 'H2.txt', 'H2O.txt', 'H_.txt', 'He.txt', 'K+.txt', \
                                                    'K.txt', 'NH3.txt', 'Na+.txt', 'Na.txt', 'TiO.txt', 'VO.txt', 'e_.txt']
listdatacomp = []
for strg in liststrgcomp:
    path = pathdata + 'ascii_output/pdependent_abundances/' + strg
    listdatacomp.append(np.loadtxt(path))

## contibution/PT/abun
figr = plt.figure(figsize=(12, 4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
axis = [[] for a in range(2)]
for a in range(2):
    axis[a] = plt.subplot(gs[0, a])
### contribution function
### pressure temperature
#axis[1].plot(dataptem[1, :], dataptem[0, :])
#axis[1].plot(dataptem[2, :], dataptem[0, :])
#axis[1].plot(dataptem[3, :], dataptem[0, :])
numbsamp = dataconf.shape[0] - 1
indxsamp = np.arange(numbsamp)
for i in indxsamp[::100]:
    axis[0].plot(dataptem[i, :], dataptem[0, :], color='b', alpha=0.1)
axis[0].set_xlabel('$T$ [K]')
axis[0].set_yscale('log')
axis[0].set_ylabel('$P$ [mbar]')
axis[0].invert_yaxis()
axistwin = axis[0].twinx()
axistwin.plot(dataconf[100, :], dataconf[0, :])
axistwin.set_yticklabels([])
## abundance
numbcomp = len(listdatacomp)
indxcomp = np.arange(numbcomp)
for k in indxcomp:
    axis[1].plot(listdatacomp[k][:, 1], listdatacomp[k][:, 0])
axis[1].set_yticklabels([])
axis[1].set_xscale('log')
axis[1].set_xlabel('VMR')
axis[1].set_yscale('log')
axis[1].invert_yaxis()
plt.subplots_adjust(hspace=0., wspace=0.)
path = pathimag + 'ptem.pdf'
print 'Writing to %s...' % path
plt.savefig(path)
plt.close()


# cross-section
figr, axis = plt.subplots(figsize=(6, 4))
for k in indxcomp:
    print 'k'
    print k
    print 'listdatacomp[k][:, 1]'
    summgene(listdatacomp[k][:, 1])
    print 'listdatacomp[k][:, 0]'
    summgene(listdatacomp[k][:, 0])
    #axis.plot(listdatacomp[k][:, 1], listdatacomp[k][:, 0])
axis.set_xlabel(r'$\lambda$ [$\mu$m]$')
axis.set_ylabel(r'$\sigma$ [cm$^2$]')
axis.set_yscale('log')
path = pathimag + 'csec.pdf'
print 'Writing to %s...' % path
plt.savefig(path)
plt.close()


#thpt = scipy.interpolate.interp1d(meanlambband, thptband)(meanlamb)
#specplan = thpt / meanlamb**5 / (np.exp(cons / meanlamb / tempplan) - 1.)
#specstar = thpt / meanlamb**5 / (np.exp(cons / meanlamb / tempstar) - 1.)
#fluxplan = np.sum(difflamb * specplan)
#fluxstar = np.sum(difflamb * specstar)
#fluxrati = fluxplan / fluxstar
#
## plotting
### Half of the duration plotting window
## temp
#timetole = 3. * dura

# type of sampling
samptype = 'emce'

alph = 0.2
        
listaminchi2runs = []
listsamppararuns = []
listlistlablpara = []
listlistlablparafull = []
liststrgruns = []
listmodltyperuns = []
    
# get data
path = pathdata + 'TESS.csv'
arryfrst = np.loadtxt(path, delimiter=',', skiprows=1)
    
liststrgmask = ['umsk', 'ful1', 'ful2', 'full']

for strgmask in liststrgmask:
    print 'strgmask'
    print strgmask
    
    # mask the data, if necessary
    if strgmask == 'umsk':
        indxgood = np.arange(arryfrst.shape[0])
    else:
        indxgood = retr_mask(arryfrst[:, 0], indxgood, strgmask)
    arryseco = arryfrst[indxgood, :]
    if strgmask != 'umsk':
        path = pathdata + 'TESS_%s.csv' % strgmask
        print 'Writing to %s...' % path
        np.savetxt(path, arryseco, delimiter=',')
            
    if strgmask == 'full':
        liststrgbins = ['0100', 'ubnd']
    else:
        liststrgbins = ['ubnd']

    ## plot the light curve
    figr, axis = plt.subplots(figsize=(12, 4))
    axis.set_ylabel('Relative flux')
    axis.set_xlabel('Time [day]')
    axis.plot(arryseco[:, 0], arryseco[:, 1], ls='', marker='o', markersize=5, alpha=0.3)
    temp, listcaps, temp = axis.errorbar(arryseco[:, 0],arryseco[:, 1], yerr=arryseco[:, 2], \
                                                                    color='black', ls='', markersize=4, marker='o', lw=3, alpha=alph)
    for caps in listcaps:
        caps.set_markeredgewidth(3)
    plt.tight_layout()
    path = pathimag + 'lcur_%s.pdf' % (strgmask)
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()

    ## plot Lomb Scargle periodogram
    figr, axis = plt.subplots(figsize=(12, 4))
    axis.set_ylabel('Power')
    axis.set_xlabel('Frequency [1/day]')
    arryfreq = np.linspace(0.1, 10., 2000)
    for a in range(2):
        indxtemp = np.arange(arryseco.shape[0])
        if a == 0:
            colr = 'g'
        if a == 1:
            colr = 'r'
            for k in range(1400, 1500):
                indxtemp = np.setdiff1d(indxtemp, np.where(abs(arryseco[:, 0] - k * peri - epoc) < dura * 2)[0])
        ydat = scipy.signal.lombscargle(arryseco[indxtemp, 0], arryseco[indxtemp, 1], arryfreq)
        axis.plot(arryfreq * 2. * np.pi, ydat, ls='', marker='o', markersize=5, alpha=0.3, color=colr)
    for a in range(4):
        axis.axvline(a / peri, ls='--', color='black')
    plt.tight_layout()
    path = pathimag + 'lspd_%s.pdf' % (strgmask)
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()

    for strgbins in liststrgbins:

        print 'strgbins'
        print strgbins
        
        if strgbins != 'ubnd':
            numbbins = int(strgbins)
        else:
            numbbins = None

        # bin the data, if necessary
        if strgbins != 'ubnd':
            arrythrd = tesstarg.util.flbn(arryseco, epoc, peri, numbbins)
        else:
            arrythrd = tesstarg.util.fold(arryseco, epoc, peri)
            
        # parse the data
        meanphas = (arrythrd[:, 0] - epoc) / peri
        data = arrythrd[:, 1]
        datastdv = arrythrd[:, 2]
        
        meanphas = (meanphas) % 1.
        indxsort = np.argsort(meanphas)
        meanphas = meanphas[indxsort]
        data = data[indxsort]
        datastdv = datastdv[indxsort]
        
        # mask out the primary transit
        indx = np.where(abs(meanphas - 0.5) < (1. - timetole / peri) / 2.)[0]
        data = data[indx]
        datastdv = datastdv[indx]
        meanphas = meanphas[indx]
        
        numbphas = data.size
        indxphas = np.arange(numbphas)
        
        # make plot
        arrythrd

        # list of models
        if strgmask == 'full' and strgbins == 'ubnd':
            listmodltype = ['simp', 'shft']
        
        for modltype in listmodltype:
            if modltype == 'simp':
                numbpara = 10
            if modltype == 'shft':
                numbpara = 11
            
            listmodltyperuns.append('%s' % modltype)
            if strgbins != 'ubnd':
                strgruns = 'binned'
            if modltype == 'umsk':
                strgruns = 'no mask'
            if modltype == 'ful1':
                strgruns = 'first orbit'
            if modltype == 'ful1':
                strgruns = 'second orbit'
            if modltype == 'shft':
                strgruns = 'phase offset'
            liststrgruns.append(strgruns)
            if samptype == 'emce':
                numbsamp = 10000
                numbsampburn = 1000
            
            numbdoff = numbphas - numbpara
            
            listlablpara = ['Constant', 'A1 Amplitude [ppm]', 'A2 Amplitude [ppm]', 'A3 Amplitude [ppm]', \
                            'Excess emission amplitude', 'Reflected B1 amplitude [ppm]', 'B2 Amplitude [ppm]', \
                            'B3 Amplitude [ppm]', 'Dayside temperature [K]', 'Nightside temperature [K]']
            liststrgpara = ['cons', 'ama1', 'ama2', 'ama3', 'exce', 'amb1', 'amb2', 'amb3', 'dayt', 'nigt']
            if modltype == 'shft':
                listlablpara += ['Phase offset [degree]']
                liststrgpara += ['offs']
            listlablparafull = listlablpara[:]
            liststrgparafull = liststrgpara[:]
            listlablparafull += ['Secondary transit depth [ppm]', 'Night side [ppm]', 'Atmospheric [ppm]', 'Geometric albedo']
            liststrgparafull += ['seco', 'nigh', 'atmo', 'galb']
            
            numbparafull = len(liststrgparafull)
            indxpara = np.arange(numbpara)
            indxparafull = np.arange(numbparafull)
            limtpara = np.empty((2, numbpara))
            # offs
            limtpara[0, 0] = 0.9
            limtpara[1, 0] = 1.1
            # amfo
            limtpara[0, 1:8] = -1e3
            limtpara[1, 1:8] = 1e3
            # tempdayy
            limtpara[0, 8] = 0.
            limtpara[1, 8] = 5000.
            # tempstar
            limtpara[0, 9] = 0.
            limtpara[1, 9] = 10000.
            if modltype == 'shft':
                # phas
                limtpara[0, 10] = -10.
                limtpara[1, 10] = 10.
            
            indxoccl = np.where(abs(meanphas - 0.5) < dura / peri)[0]
            dictllik = [meanphas, data, datastdv, modltype, indxoccl]
            dicticdf = []
            
            if samptype == 'emce':
                numbwalk = 50
                indxwalk = np.arange(numbwalk)
                parainit = []
                for k in indxwalk:
                    parainit.append(np.empty(numbpara))
                    meannorm = (limtpara[0, :] + limtpara[1, :]) / 2.
                    stdvnorm = (limtpara[0, :] - limtpara[1, :]) / 10.
                    parainit[k]  = (scipy.stats.truncnorm.rvs((limtpara[0, :] - meannorm) / stdvnorm, \
                                                                (limtpara[1, :] - meannorm) / stdvnorm)) * stdvnorm + meannorm
                numbsampwalk = numbsamp / numbwalk
                numbsampwalkburn = numbsampburn / numbwalk
                objtsamp = emcee.EnsembleSampler(numbwalk, numbpara, retr_lpos, args=dictllik)
                parainitburn, prob, state = objtsamp.run_mcmc(parainit, numbsampwalkburn)
                objtsamp.reset()
                objtsamp.run_mcmc(parainitburn, numbsampwalk)
                objtsave = objtsamp
            else:
            
                sampler = dynesty.NestedSampler(retr_llik, icdf, numbpara, logl_args=dictllik, ptform_args=dicticdf, bound='single', dlogz=1000.)
                sampler.run_nested()
                results = sampler.results
                results.summary()
                objtsave = results
                
            if samptype == 'emce':
                numbsamp = objtsave.flatchain.shape[0]
                indxsampwalk = np.arange(numbsampwalk)
            else:
                numbsamp = objtsave['samples'].shape[0]
            
            indxsamp = np.arange(numbsamp)
            
            # resample the nested posterior
            if samptype == 'nest':
                weights = np.exp(results['logwt'] - results['logz'][-1])
                samppara = dynesty.utils.resample_equal(results.samples, weights)
                assert samppara.size == results.samples.size
            
            if samptype == 'emce':
                arrysamp = objtsave.flatchain
            else:
                arrysamp = samppara
            
            sampllik = objtsave.lnprobability.flatten()
            aminchi2 = (-2. * np.amax(sampllik) / numbdoff)
            listaminchi2runs.append(aminchi2)
            
            arrysamptemp = np.copy(arrysamp)
            arrysamp = np.empty((numbsamp, numbparafull))
            arrysamp[:, :-(numbparafull - numbpara)] = arrysamptemp
            for k in indxsamp:
                arrysamp[k, -4] = retr_occl(arrysamptemp[k, 8], tempstar)
                arrysamp[k, -3] = retr_occl(arrysamptemp[k, 9], tempstar)
                arrysamp[k, -2] = arrysamp[k, -4] - arrysamp[k, -3]
                arrysamp[:, -1] = (smaj**2 / radiplan**2) * (1e-6 * arrysamp[:, 5])

            listsamppararuns.append(arrysamp)
            listlistlablpara.append(listlablpara)
            listlistlablparafull.append(listlablparafull)

            # plot the posterior
            
            ### histogram
            for k in indxparafull:
                figr, axis = plt.subplots()
                if samptype == 'emce':
                    axis.hist(arrysamp[:, k]) 
                else:
                    axis.hist(samppara[:, k]) 
                axis.set_xlabel(listlablparafull[k])
                path = pathimag + 'diag/hist_%s_%s_%s_%s.pdf' % (liststrgparafull[k], modltype, strgmask, strgbins)
                plt.tight_layout()
                print 'Writing to %s...' % path
                plt.savefig(path)
                plt.close()
            
            if samptype == 'nest':
                for keys in objtsave:
                    if isinstance(objtsave[keys], np.ndarray) and objtsave[keys].size == numbsamp:
                        figr, axis = plt.subplots()
                        axis.plot(indxsamp, objtsave[keys])
                        path = pathimag + '%s_%s.pdf' % (keys, modltype)
                        print 'Writing to %s...' % path
                        plt.savefig(path)
            else:
                ## log-likelihood
                figr, axis = plt.subplots()
                if samptype == 'emce':
                    for i in indxwalk:
                        axis.plot(indxsampwalk[::10], objtsave.lnprobability[::10, i])
                else:
                    axis.plot(indxsamp, objtsave['logl'])
                path = pathimag + 'diag/llik_%s_%s_%s.pdf' % (modltype, strgmask, strgbins)
                plt.tight_layout()
                print 'Writing to %s...' % path
                plt.savefig(path)
                plt.close()
            
                chi2 = -2. * objtsave.lnprobability
                print 'Posterior-mean chi2: '
                print np.mean(chi2)
                print 'Posterior-mean chi2 per dof: '
                print np.mean(chi2) / numbdoff
            
            ### sample model phas
            numbsampplot = 100
            indxsampplot = np.random.choice(indxsamp, numbsampplot, replace=False)
            yerr = datastdv
            
            numbphasfine = 1000
            meanphasfine = np.linspace(0.1, 0.9, numbphasfine)
            phasmodlfine = np.empty((numbsamp, numbphasfine))
            indxphasfineoccl = np.where(abs(meanphasfine - 0.5) < dura / peri)[0]
            for k, indxsampplottemp in enumerate(indxsampplot):
                if samptype == 'emce':
                    objttemp = objtsave.flatchain
                else:
                    objttemp = samppara
                offs = objttemp[indxsampplottemp, 0]
                amfo = objttemp[indxsampplottemp, 1:8]
                tempdayy = objttemp[indxsampplottemp, 8]
                tempstar = objttemp[indxsampplottemp, 9]
                if modltype == 'shft':
                    shft = objttemp[indxsampplottemp, 10]
                else:
                    shft = np.zeros_like(tempdayy)
                phasmodlfine[k, :], deptatmofine = retr_modl(meanphasfine, offs, amfo, tempdayy, tempstar, shft, modltype, indxphasfineoccl)
            
            figr, axis = plt.subplots(figsize=(12, 4))
            axis.errorbar(meanphas, data, yerr=yerr, color='black', marker='o', ls='', alpha=0.2, markersize=1)
            for k, indxsampplottemp in enumerate(indxsampplot):
                axis.plot(meanphasfine, phasmodlfine[k, :], alpha=0.05, color='b')
            axis.set_xlabel('Phase')
            axis.set_ylabel('Relative flux')
            plt.tight_layout()
            path = pathimag + 'phas_modl_%s_%s_%s.pdf' % (modltype, strgmask, strgbins)
            print 'Writing to %s...' % path
            plt.savefig(path)
            plt.close()
            
            ### nested sampling specific
            if samptype == 'nest':
                figr, axes = dynesty.plotting.runplot(results)
                path = pathimag + 'dyne_runs_%s.pdf' % modltype
                print 'Writing to %s...' % path
                plt.savefig(path)
                plt.close()
                
                figr, axes = dynesty.plotting.traceplot(results)
                path = pathimag + 'dyne_trac_%s.pdf' % modltype
                print 'Writing to %s...' % path
                plt.savefig(path)
                plt.close()
                
                figr, axes = dynesty.plotting.cornerplot(results)
                path = pathimag + 'dyne_corn_%s.pdf' % modltype
                print 'Writing to %s...' % path
                plt.savefig(path)
                plt.close()

numbruns = len(liststrgruns)
indxruns = np.arange(numbruns)
for b in indxruns:
    tdpy.mcmc.plot_grid(pathimag, '%d' % b, listsamppararuns[b], listlistlablparafull[b])

fileoutp = open(pathdata + 'post.csv', 'w')

numbrrunsdivi = 3
numbdivi = numbruns / numbrrunsdivi + 1
indxdivi = np.arange(numbdivi)
for a in indxdivi:
    numbleft = numbruns - numbrrunsdivi * a
    if numbleft >= numbrrunsdivi:
        numbrunsdivithis = numbrrunsdivi
    else:
        numbrunsdivithis = numbleft
    indxrunsdivithis = np.arange(numbrunsdivithis)
    fileoutp.write(' & ')
    for k in indxrunsdivithis:
        b = a * numbrrunsdivi + k
        fileoutp.write('%s' % liststrgruns[b])
        if k != numbrunsdivithis - 1:
            fileoutp.write(' & ')
    fileoutp.write('\\\\\n')
    fileoutp.write('\\hline\n')
    fileoutp.write('\\hline\n')
    fileoutp.write('\\hline\n')
    fileoutp.write('$\chi^2_{\\nu}$ & ')
    for k in indxrunsdivithis:  
        b = a * numbrrunsdivi + k
        fileoutp.write('$%.3g$' % listaminchi2runs[b])
        if k != numbrunsdivithis - 1:
            fileoutp.write(' & ')
    fileoutp.write('\\hline\n')
    fileoutp.write('\\hline\n')
    fileoutp.write('\\\\\n')
    for n in indxparafull:
        fileoutp.write('%s & ' % listlablparafull[n])
        for k in indxrunsdivithis:
            b = a * numbrrunsdivi + k
            if len(listlistlablparafull[b]) < numbparafull:
                continue

            ydat = np.median(listsamppararuns[b][:, n])
            uerr = np.percentile(listsamppararuns[b][:, n], 84.) - ydat
            lerr = ydat - np.percentile(listsamppararuns[b][:, n], 16.)
            fileoutp.write('$%.3g \substack{+%.3g \\\\ -%.3g}$' % (ydat, uerr, lerr))
            if k != numbrunsdivithis - 1:
                fileoutp.write(' & ')
        fileoutp.write('\\\\\n')
    fileoutp.write('\\hline\n')
fileoutp.close()


# make a plot of WASP-121b properties

# irradiation vs. 




