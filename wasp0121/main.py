import os, sys

import numpy as np
import scipy
import scipy.interpolate
from scipy.interpolate import UnivariateSpline
from scipy.signal import lombscargle

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import seaborn as sns

from astropy import constants as c

# allesfitter
import allesfitter
from allesfitter.exoworlds_rdx.lightcurves.index_transits import index_eclipses
from allesfitter.exoworlds_rdx.lightcurves.lightcurve_tools import rebin_err, phase_fold

# own modules
import tdpy.mcmc
from tdpy.util import summgene
import ephesus.util
from tdpy.util import prnt_list

# Max's modules
from exoworlds.tess import extract_SPOC_data, extract_QLP_data

import emcee
import dynesty

matplotlib.rc('text', usetex = True)
matplotlib.rc('font', **{'family' : "sans-serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)


def retr_lpos_albb(para, gdat):
    
    albb = para[0]
    epsi = para[1]
   
    if albb < 0 or albb > 1 or epsi > 1 or epsi < 0:
        lpos = -np.inf
    else:
        psiimodl = (1 - albb)**.25 * (1. - 5. / 8. * epsi)**.25 / 1.5**.25
        lpos = gdat.likeintp(psiimodl)

    return lpos


def retr_lpos_four(para, gdat, modltype):
    
    cons = para[0]
    deptnigh =para[1]
    amfo = para[2:8]
    if modltype == 'shft':
        shft = para[8]
    else:
        shft = 0.
    
    if ((para < gdat.limtpara[0, :]) | (para > gdat.limtpara[1, :])).any():
        lpos = -np.inf
    else:
        llik, deptpmod = retr_llik_four(gdat, cons, deptnigh, amfo, shft, modltype)
        #lpri = - 0.5 * ((amfo[3] / 0.1)**2)
        #lpos = llik + lpri
        lpos = llik

    return lpos


def retr_specbbod(gdat, tmpt, wlen):
    
    spec = 3.742e11 / wlen**5 / (np.exp(gdat.conswlentmpt / (wlen * tmpt)) - 1.)
    
    return spec


def retr_modl_spec(gdat, tmpt, booltess=False, strgtype='intg'):
    
    if booltess:
        thpt = scipy.interpolate.interp1d(gdat.meanwlenband, gdat.thptband)(wlen)
    else:
        thpt = 1.
    
    if strgtype == 'intg':
        spec = retr_specbbod(gdat, tmpt, gdat.meanwlen)
        spec = np.sum(gdat.diffwlen * spec)
    if strgtype == 'diff' or strgtype == 'logt':
        spec = retr_specbbod(gdat, tmpt, gdat.cntrwlen)
        if strgtype == 'logt':
            spec *= gdat.cntrwlen
    
    return spec


def retr_lpos_spec(para, gdat):
    
    if ((para < gdat.limtpara[0, :]) | (para > gdat.limtpara[1, :])).any():
        lpos = -np.inf
    else:
        tmpt = para[0]
        specboloplan = retr_modl_spec(gdat, tmpt, booltess=False, strgtype='intg')
        deptplan = gdat.radirati**2 * specboloplan / gdat.specstarintg * 1e6
        llik = -0.5 * np.sum((deptplan - gdat.deptobsd)**2 / gdat.varideptobsd)
        #print('tmpt')
        #print(tmpt)
        #print('deptplan')
        #print(deptplan)
        #print('gdat.deptobsd')
        #print(gdat.deptobsd)
        #print('gdat.varideptobsd')
        #print(gdat.varideptobsd)
        #print
        lpos = llik
    
    return lpos


def retr_llik_four(gdat, cons, deptnigh, amfo, shft, modltype):
    
    modl, deptpmod = retr_modl_four(gdat, gdat.meanphas, gdat.indxphasseco, cons, deptnigh, amfo, shft, modltype)
    
    llik = -0.5 * np.sum((modl - gdat.data)**2 / gdat.datastdv**2)
    
    return llik, deptpmod


def retr_modl_four(gdat, phas, indxphasseco, cons, deptnigh, amfo, shft, modltype):
    
    phasshft = phas + shft * np.pi / 180.
    deptpmod = amfo[0]
    modl = cons + deptnigh * 1e-6 + 0.5 * 1e-6 * deptpmod * np.cos(phasshft * 2. * np.pi)
    for k in range(2):
        modl += 0.5e-6 * amfo[k+1] * np.cos((k + 2) * phas * 2. * np.pi)
    for k in range(3):
        modl += 0.5e-6 * amfo[k+3] * np.sin((k + 1) * phas * 2. * np.pi) 
    deptseco = -deptpmod + deptnigh
    modl[indxphasseco] -= deptseco * 1e-6
    
    return modl, deptpmod


def icdf(para):
    
    icdf = gdat.limtpara[0, :] + para * (gdat.limtpara[1, :] - gdat.limtpara[0, :])

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


def init():

    gdat = tdpy.util.gdatstrt()
    
    # paths
    pathbase = os.environ['TDGU_DATA_PATH'] + '/wasp0121/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    
    # plotting
    gdat.figrsize = np.empty((5, 2))
    gdat.figrsize[0, :] = np.array([12., 4.])
    gdat.figrsize[1, :] = np.array([12., 6.])
    gdat.figrsize[2, :] = np.array([12., 10.])
    gdat.figrsize[3, :] = np.array([12., 14.])
    gdat.figrsize[4, :] = np.array([6., 6.])
    boolpost = False
    if boolpost:
        gdat.figrsize /= 1.5
    
    ## TESS throughput 
    gdat.data = np.loadtxt(pathdata + 'band.csv', delimiter=',', skiprows=9)
    gdat.meanwlenband = gdat.data[:, 0] * 1e-3
    gdat.thptband = gdat.data[:, 1]
    
    ## other data
    massplan = 1.184 # [MJ]
    massstar = 1417. # [MJ]
    radistar = 14.51 # [RJ]
    incl = 87.6 * np.pi / 180.
    smax = 54.437 # [RJ]
   
    smajaxis = 0.02544 * 2092.51 # [RJ]
    radiplan = 1.865 # [RJ]
    radistar = 14.51 # [RJ]
    gdat.radirati = radiplan / radistar

    gdat.timetess = 2457000.

    ## parameters
    # WASP-121
    ## period
    gdat.peri = 1.2749255 # [day]
    ## epoch
    gdat.epoc = 2456635.70832 # [BJD]
    ## Half of the transit duration
    dura = 0.1203 / 2. # [day] Delrez2016
    
    ## expected ellipsoidal variation (EV) and Doppler beaming (DB)
    print('Predicting the ellipsoidal variation and Doppler beaming amplitudes...')
    ### EV
    #### limb and gravity darkening coefficients from Claret2017
    u = 0.4
    g = 0.2
    alphelli = 0.15 * (15 + u) * (1 + g) / (3 - u)
    deptelli = alphelli * massplan * np.sin(incl)**2 / massstar * (radistar / smax)**3
    ### DB
    kkrv = 181. # [m/s]
    binswlenbeam = np.linspace(0.6, 1., 101)
    meanwlenbeam = (binswlenbeam[1:] + binswlenbeam[:-1]) / 2.
    diffwlenbeam = (binswlenbeam[1:] - binswlenbeam[:-1]) / 2.
    x = 2.248 / meanwlenbeam
    f = .25 * x * np.exp(x) / (np.exp(x) - 1.)
    deptbeam = 4. * kkrv / 3e8 * np.sum(diffwlenbeam * f)

    # extract SPOC data
    extract_SPOC_data([pathdata + 'data_preparation/archival_data/tess2019006130736-s0007-0000000022529346-0131-s_lc.fits'], 
                      outdir=pathdata+'data_preparation/TESS_SAP', PDC=False, auto_correct_dil=True, extract_centd=True, extract_dil=True)
    
    extract_SPOC_data([pathdata + 'data_preparation/archival_data/tess2019006130736-s0007-0000000022529346-0131-s_lc.fits'], 
                      outdir=pathdata+'data_preparation/TESS_PDCSAP', PDC=True, auto_correct_dil=True, extract_centd=True, extract_dil=True)
    
    
    # plot PDCSAP and SAP light curves
    gdat.timepdcc, lcurpdcc, stdvlcurpdcc = np.genfromtxt(pathdata + 'data_preparation/TESS_PDCSAP/TESS.csv', delimiter=',', unpack=True)
    gdat.timesapp, lcursapp, stdvlcursapp = np.genfromtxt(pathdata + 'data_preparation/TESS_SAP/TESS.csv', delimiter=',', unpack=True)
    
    figr, axis = plt.subplots(2, 1, figsize=gdat.figrsize[1, :])
    axis[0].plot(gdat.timepdcc - gdat.timetess, lcurpdcc, 'k.')
    axis[1].plot(gdat.timesapp - gdat.timetess, lcursapp, '.', color='k')
    for a in range(2):
        axis[a].set_ylabel('Relative Flux')
    #axis[0].set_xticklabels([])
    axis[1].set_xlabel('Time [BJD - 2457000]')
    path = pathimag + 'lcur.pdf'
    plt.subplots_adjust(hspace=0.)
    plt.savefig(path)
    plt.close()
    
    # fit a spline to the SAP light curve
    print('Fitting a spline to the SAP light curve...')
    timesplnedge = [0., 2458504., np.inf]
    numbsplnedge = len(timesplnedge)
    numbsplnknot = numbsplnedge - 1
    indxsplnknot = np.arange(numbsplnknot)

    # produce a table for the spline coefficients
    fileoutp = open(pathdata + 'coef.csv', 'w')
    fileoutp.write(' & ')
    
    figr, axis = plt.subplots(2, 1, figsize=gdat.figrsize[1, :])
    gdat.lcurspln = []
    lcursplnbase = []
    listobjtspln = []
    wdth = 4. / 24. # [days]
    for i in indxsplnknot:
        indxtimesplnregi = np.where((gdat.timesapp >= timesplnedge[i]) & (gdat.timesapp <= timesplnedge[i+1]))[0]
        timesplnregi = gdat.timesapp[indxtimesplnregi]
        lcursplnbaseregi = lcursapp[indxtimesplnregi]
        flux_err1 = stdvlcursapp[indxtimesplnregi]
        temp, indxtimesplntran, indxtimesplnoutt = index_eclipses(timesplnregi, gdat.epoc, gdat.peri, width_1=wdth, width_2=wdth)
        objtspln = UnivariateSpline(timesplnregi[indxtimesplnoutt], lcursplnbaseregi[indxtimesplnoutt]-1.)
        timesplnregifine = np.linspace(timesplnregi[0], timesplnregi[-1], 1000)
        gdat.lcurspln += list(lcursplnbaseregi-objtspln(timesplnregi))
        lcursplnbase += list(objtspln(timesplnregi))
        listobjtspln.append(objtspln)
    
        print('i')
        print(i)
        print('$\beta$:', listobjtspln[i].get_coeffs())
        print('$t_k$:', listobjtspln[i].get_knots())
        
        # plot the masked and detrended light curves
        axis[0].plot(timesplnregi[indxtimesplnoutt] - gdat.timetess, lcursplnbaseregi[indxtimesplnoutt], 'k.', color='k')
        axis[0].plot(timesplnregifine - gdat.timetess, listobjtspln[i](timesplnregifine)+1., 'b-', lw=3)
        axis[1].plot(timesplnregi - gdat.timetess, lcursplnbaseregi-listobjtspln[i](timesplnregi), '.', color='k')
    for a in range(2):
        axis[a].set_ylabel('Relative Flux')
    axis[0].set_xticklabels([])
    axis[1].set_xlabel('Time [BJD - 2457000]')
    path = pathimag + 'lcur_spln.pdf'
    plt.subplots_adjust(hspace=0.)
    plt.savefig(path)
    plt.close()
   


    figr, axis = plt.subplots(4, 1, figsize=gdat.figrsize[2, :])
    axis[0].plot(gdat.timesapp - gdat.timetess, lcursapp, '.', color='k')
    axis[1].plot(gdat.timepdcc - gdat.timetess, lcurpdcc, 'k.')
    for a in range(4):
        axis[a].set_ylabel('Relative Flux')
    #axis[0].set_xticklabels([])
    
    # fit a spline to the SAP light curve
    print('Fitting a spline to the SAP light curve...')
    timesplnedge = [0., 2458504., np.inf]
    numbsplnedge = len(timesplnedge)
    numbsplnknot = numbsplnedge - 1
    indxsplnknot = np.arange(numbsplnknot)

    # produce a table for the spline coefficients
    fileoutp = open(pathdata + 'coef.csv', 'w')
    fileoutp.write(' & ')
    
    gdat.lcurspln = []
    lcursplnbase = []
    listobjtspln = []
    wdth = 4. / 24. # [days]
    for i in indxsplnknot:
        indxtimesplnregi = np.where((gdat.timesapp >= timesplnedge[i]) & (gdat.timesapp <= timesplnedge[i+1]))[0]
        timesplnregi = gdat.timesapp[indxtimesplnregi]
        lcursplnbaseregi = lcursapp[indxtimesplnregi]
        flux_err1 = stdvlcursapp[indxtimesplnregi]
        temp, indxtimesplntran, indxtimesplnoutt = index_eclipses(timesplnregi, gdat.epoc, gdat.peri, width_1=wdth, width_2=wdth)
        objtspln = UnivariateSpline(timesplnregi[indxtimesplnoutt], lcursplnbaseregi[indxtimesplnoutt]-1.)
        timesplnregifine = np.linspace(timesplnregi[0], timesplnregi[-1], 1000)
        gdat.lcurspln += list(lcursplnbaseregi-objtspln(timesplnregi))
        lcursplnbase += list(objtspln(timesplnregi))
        listobjtspln.append(objtspln)
    
        print('i')
        print(i)
        print('$\beta$:', listobjtspln[i].get_coeffs())
        print('$t_k$:', listobjtspln[i].get_knots())
        
        # plot the masked and detrended light curves
        axis[2].plot(timesplnregi[indxtimesplnoutt] - gdat.timetess, lcursplnbaseregi[indxtimesplnoutt], 'k.', color='k')
        axis[2].plot(timesplnregifine - gdat.timetess, listobjtspln[i](timesplnregifine)+1., 'b-', lw=3)
        axis[3].plot(timesplnregi - gdat.timetess, lcursplnbaseregi-listobjtspln[i](timesplnregi), '.', color='k')
    for a in range(2):
        axis[a].set_ylabel('Relative Flux')
    axis[2].set_xticklabels([])
    axis[3].set_xlabel('Time [BJD - 2457000]')
    path = pathimag + 'lcur_totl.pdf'
    plt.subplots_adjust(hspace=0.)
    plt.savefig(path)
    plt.close()
   










    gdat.lcurspln = np.array(gdat.lcurspln)
    lcursplnbase = np.array(lcursplnbase)
    temp, indxtimesplntran, indxtimesplnoutt = index_eclipses(gdat.timesapp, gdat.epoc, gdat.peri, width_1=4./24., width_2=2./24.)
    offset = np.mean(gdat.lcurspln[indxtimesplntran]) - 1.
    gdat.lcurspln -= offset
    
    fileoutp.write('\\hline\n')
    fileoutp.close()
    
    gdat.time = gdat.timesapp

    ## phase-fold and save the detrended light curve
    dt = 0.01
    ferr_type = 'medsig'
    ferr_style = 'sem'
    sigmaclip = False
        
    gdat.timebind, gdat.lcursplnbind, gdat.stdvlcursplnbind, N = rebin_err(gdat.time, gdat.lcurspln, \
                                                                                ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    gdat.phasbind, gdat.pcursplnbind, gdat.stdvpcursplnbind, N, gdat.phas = phase_fold(gdat.time, gdat.lcurspln, gdat.peri, gdat.epoc, \
                                                                                ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    
    data = np.column_stack((gdat.time, gdat.lcurspln, stdvlcursapp))
    path = pathdata + 'data_preparation/TESS.csv'
    print('Writing to %s' % path)
    np.savetxt(path, data, delimiter=',', header='time,flux,flux_err')
    
    strgalletype = 'offs'
    
    # read the allesfitter posterior
    companion = 'b'
    strginst = 'TESS'
    pathalle = pathdata + 'allesfits/allesfit_%s/' % strgalletype
    alles = allesfitter.allesclass(pathalle)
    

    listradirati = alles.posterior_params['b_rr']
    listrsma = alles.posterior_params['b_rsuma']
    if strgalletype == 'ther':
        listphasshft = 360. * alles.posterior_params['b_thermal_emission_timeshift_TESS'] / alles.posterior_params['b_period']
        
        prnt_list(listphasshft, 'Phase shift')
        return

    else:
        listalbgalle = alles.posterior_params['b_geom_albedo_TESS']
        listdeptnigh = alles.posterior_params['b_sbratio_TESS'] * listradirati**2
    
    numbsamp = listrsma.size

    listradirati = listradirati[int(numbsamp/4):]
    listalbgalle = listalbgalle[int(numbsamp/4):]
    listdeptnigh = listdeptnigh[int(numbsamp/4):]
    listrsma = listrsma[int(numbsamp/4):]
    numbsamp = listrsma.size
    
    # calculate nightside, secondary and planetary modulation from allesfitter output
    ## what allesfitter calls 'geometric albedo' is not the actual geometric albedo
    listdeptpmod = listalbgalle * (listradirati * listrsma / (1. + listradirati))**2
    listsamp = np.empty((numbsamp, 6))
    listsamp[:, 0] = listdeptnigh * 1e6
    listsamp[:, 1] = (listdeptnigh + listdeptpmod) * 1e6
    listsamp[:, 2] = listdeptpmod * 1e6
    listsamp[:, 3] = listdeptpmod * np.random.rand(listdeptpmod.size) * 1e6
    listsamp[:, 4] = listdeptpmod * 1e6 - listsamp[:, 3]
    listsamp[:, 5] = listsamp[:, 4] * 1e-6 * (1. + listradirati)**2 / listrsma**2 / listradirati**2
    listalbgtess = listsamp[:, 5]
    listlabl = ['Nightside [ppm]', 'Secondary [ppm]', 'Modulation [ppm]', 'Thermal [ppm]', 'Reflected [ppm]', 'Geometric Albedo']
    tdpy.mcmc.plot_grid(pathimag, 'post_alle', listsamp, listlabl, plotsize=2.5)
    
    listdeptseco = listsamp[:, 1]
    medideptseco = np.median(listsamp[:, 1])
    stdvdeptseco = np.std(listsamp[:, 1])
    medideptnigh = np.median(listsamp[:, 0])
    stdvdeptnigh = np.std(listsamp[:, 0])

    path = pathdata + 'PC-Solar-NEW-OPA-TiO-LR.dat'
    arryvivi = np.loadtxt(path, delimiter=',')
    phasvivi = (arryvivi[:, 0] / 360. + 0.75) % 1. - 0.25
    deptvivi = arryvivi[:, 4] * 1e6
    indxphasvivisort = np.argsort(phasvivi)
    phasvivi = phasvivi[indxphasvivisort]
    deptvivi = deptvivi[indxphasvivisort]
    print('phasvivi')
    summgene(phasvivi)
    path = pathdata + 'PC-Solar-NEW-OPA-TiO-LR-AllK.dat'
    arryvivi = np.loadtxt(path, delimiter=',')
    wlenvivi = arryvivi[:, 1]
    specvivi = arryvivi[:, 2]

    # plot a lightcurve from the posteriors
    gdat.lcurmodl = alles.get_posterior_median_model(strginst, 'flux', xx=gdat.time)
    gdat.lcurgpro = alles.get_posterior_median_baseline(strginst, 'flux', xx=gdat.time)
    
    gdat.lcurdetr = gdat.lcurspln - gdat.lcurgpro

    gdat.timebind, gdat.lcurmodlbind, gdat.stdvlcurmodlbind, N = rebin_err(gdat.time, gdat.lcurmodl, \
                                                                            ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    gdat.phasbind, gdat.pcurmodlbind, gdat.stdvpcurmodlbind, N, gdat.phas = phase_fold(gdat.time, gdat.lcurmodl, gdat.peri, gdat.epoc, \
                                                                            ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    
    gdat.indxphassort = np.argsort(gdat.phas)

    gdat.timebind, gdat.lcurdetrbind, gdat.stdvlcurdetrbind, N = rebin_err(gdat.time, gdat.lcurdetr, \
                                                                            ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    gdat.phasbind, gdat.pcurdetrbind, gdat.stdvpcurdetrbind, N, gdat.phas = phase_fold(gdat.time, gdat.lcurdetr, gdat.peri, gdat.epoc, \
                                                                            ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    
    gdat.indxtimegapp = np.argmax(gdat.time[1:] - gdat.time[:-1]) + 1
    figr = plt.figure(figsize=gdat.figrsize[3, :])
    axis = [[] for k in range(3)]
    axis[0] = figr.add_subplot(3, 1, 1)
    axis[1] = figr.add_subplot(3, 1, 2)
    axis[2] = figr.add_subplot(3, 1, 3, sharex=axis[1])
    
    for k in range(len(axis)):
        
        if k > 0:
            xdat = gdat.phasbind
        else:
            xdat = gdat.time - gdat.timetess
        
        if k > 0:
            ydat = gdat.pcurdetrbind
        else:
            ydat = gdat.lcurdetr
        if k == 2:
            ydat = (ydat - 1.) * 1e6 + medideptnigh
        axis[k].plot(xdat, ydat, '.', color='grey', alpha=0.3, label='Raw data')
        
        if k > 0:
            xdat = gdat.phasbind
            ydat = gdat.pcurdetrbind
            yerr = np.copy(gdat.stdvpcurdetrbind)
        else:
            xdat = gdat.timebind - gdat.timetess
            ydat = gdat.lcurdetrbind
            yerr = gdat.stdvlcurdetrbind
        if k == 2:
            ydat = (ydat - 1.) * 1e6 + medideptnigh
            yerr *= 1e6
        axis[k].errorbar(xdat, ydat, marker='o', yerr=yerr, capsize=0, ls='', color='k', label='Binned data')
        
        if k > 0:
            xdat = gdat.phas[gdat.indxphassort]
            ydat = gdat.lcurmodl[gdat.indxphassort]
        else:
            xdat = gdat.time - gdat.timetess
            ydat = gdat.lcurmodl
        if k == 2:
            ydat = (ydat - 1.) * 1e6 + medideptnigh
        
        if k == 0:
            axis[k].plot(xdat[:gdat.indxtimegapp], ydat[:gdat.indxtimegapp], color='b', lw=2, label='Total Model')
            axis[k].plot(xdat[gdat.indxtimegapp:], ydat[gdat.indxtimegapp:], color='b', lw=2)
        else:
            axis[k].plot(xdat, ydat, color='b', lw=2, label='Model (This work)')
        
        if k == 2:
            axis[k].plot(phasvivi, deptvivi, color='orange', lw=2, label='GCM (Parmentier+2018)')

            axis[k].axhline(0., ls='-.', alpha=0.3, color='grey')

        if k == 0:
            axis[k].set(xlabel='Time (BJD)')
        if k > 0:
            axis[k].set(xlabel='Phase')
    axis[0].set(ylabel='Relative Flux')
    axis[1].set(ylabel='Relative Flux')
    axis[2].set(ylabel='Relative Flux - 1 [ppm]')
    #axis[1].set(ylim=[-800,1000])
    axis[2].set(ylim=[-400, 1000])
    
    ## plot components in the zoomed panel
    ### EV
    alles = allesfitter.allesclass(pathalle)
    alles.posterior_params_median['b_sbratio_TESS'] = 0
    #alles.settings['host_shape_TESS'] = 'sphere'
    #alles.settings['b_shape_TESS'] = 'sphere'
    alles.posterior_params_median['b_geom_albedo_TESS'] = 0
    alles.posterior_params_median['host_gdc_TESS'] = 0
    alles.posterior_params_median['host_bfac_TESS'] = 0
    gdat.lcurmodlcomp = alles.get_posterior_median_model(strginst, 'flux', xx=gdat.time)
    gdat.lcurmodlevvv = np.copy(gdat.lcurmodlcomp)
    #gdat.phasbind, gdat.pcurmodlcompbind, gdat.stdvpcurmodlcompbind, N, gdat.phas = phase_fold(gdat.time, gdat.lcurmodlcomp, gdat.peri, gdat.epoc, \
    #                                                                        ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    

    xdat = gdat.phas[gdat.indxphassort]
    ydat = (gdat.lcurmodlcomp[gdat.indxphassort] - 1.) * 1e6
    indxfrst = np.where(xdat < -0.07)[0]
    indxseco = np.where(xdat > 0.07)[0]
    axis[2].plot(xdat[indxfrst], ydat[indxfrst], lw=2, color='r', label='Ellipsoidal', ls='--')
    axis[2].plot(xdat[indxseco], ydat[indxseco], lw=2, color='r', ls='--')
    
    objtalle = allesfitter.allesclass(pathalle)
    alles.posterior_params_median['b_sbratio_TESS'] = 0
    alles.posterior_params_median['b_geom_albedo_TESS'] = 0
    alles.posterior_params_median['host_gdc_TESS'] = 0
    alles.posterior_params_median['host_bfac_TESS'] = 0

    ### planetary modulation
    alles = allesfitter.allesclass(pathalle)
    alles.posterior_params_median['b_sbratio_TESS'] = 0
    alles.settings['host_shape_TESS'] = 'sphere'
    alles.settings['b_shape_TESS'] = 'sphere'
    #alles.posterior_params_median['b_geom_albedo_TESS'] = 0
    alles.posterior_params_median['host_gdc_TESS'] = 0
    alles.posterior_params_median['host_bfac_TESS'] = 0
    gdat.lcurmodlcomp = alles.get_posterior_median_model(strginst, 'flux', xx=gdat.time)
    #gdat.phasbind, gdat.pcurmodlcompbind, gdat.stdvpcurmodlcompbind, N, gdat.phas = phase_fold(gdat.time, gdat.lcurmodlcomp, gdat.peri, gdat.epoc, \
    #                                                                        ferr_type=ferr_type, ferr_style=ferr_style, dt=dt, sigmaclip=sigmaclip)
    #axis[2].plot(gdat.phas[gdat.indxphassort], medideptnigh + (gdat.lcurmodlcomp[gdat.indxphassort] - 1.) * 1e6, \
    #                                                            lw=2, color='g', label='Planetary Modulation', ls='--', zorder=11)
    
    axis[2].legend(ncol=2)
    
    path = pathimag + 'pcur_alle.pdf'
    plt.savefig(path)
    plt.close()
    
    
    # plot the spherical limits
    figr, axis = plt.subplots(figsize=gdat.figrsize[0, :])
    
    alles = allesfitter.allesclass(pathalle)
    alles.posterior_params_median['b_sbratio_TESS'] = 0
    alles.settings['host_shape_TESS'] = 'sphere'
    alles.settings['b_shape_TESS'] = 'roche'
    alles.posterior_params_median['b_geom_albedo_TESS'] = 0
    alles.posterior_params_median['host_gdc_TESS'] = 0
    alles.posterior_params_median['host_bfac_TESS'] = 0
    lcurmodltemp = alles.get_posterior_median_model(strginst, 'flux', xx=gdat.time)
    print('(gdat.lcurmodlevvv - lcurmodltemp)[gdat.indxphassort] * 1e6')
    summgene((gdat.lcurmodlevvv - lcurmodltemp)[gdat.indxphassort] * 1e6)
    axis.plot(gdat.phas[gdat.indxphassort], (gdat.lcurmodlevvv - lcurmodltemp)[gdat.indxphassort] * 1e6, lw=2, label='Spherical star')
    
    alles = allesfitter.allesclass(pathalle)
    alles.posterior_params_median['b_sbratio_TESS'] = 0
    alles.settings['host_shape_TESS'] = 'roche'
    alles.settings['b_shape_TESS'] = 'sphere'
    alles.posterior_params_median['b_geom_albedo_TESS'] = 0
    alles.posterior_params_median['host_gdc_TESS'] = 0
    alles.posterior_params_median['host_bfac_TESS'] = 0
    lcurmodltemp = alles.get_posterior_median_model(strginst, 'flux', xx=gdat.time)
    print('(gdat.lcurmodlevvv - lcurmodltemp)[gdat.indxphassort] * 1e6')
    summgene((gdat.lcurmodlevvv - lcurmodltemp)[gdat.indxphassort] * 1e6)
    axis.plot(gdat.phas[gdat.indxphassort], (gdat.lcurmodlevvv - lcurmodltemp)[gdat.indxphassort] * 1e6, lw=2, label='Spherical planet')
    axis.legend()
    axis.set_ylim([-100, 100])
    axis.set(xlabel='Phase')
    axis.set(ylabel='Relative flux [ppm]')
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
  
    # use psi posterior to infer Bond albedo and heat circulation efficiency

    # get data
    ## from Tom
    path = pathdata + 'ascii_output/EmissionModelArray.txt'
    arrymodl = np.loadtxt(path)
    path = pathdata + 'ascii_output/EmissionDataArray.txt'
    arrydata = np.loadtxt(path)
    # update Tom's array with the new secondary depth
    arrydata[0, 2] = medideptseco
    arrydata[0, 3] = stdvdeptseco
    ### add the nightsiide emission
    arrydata = np.concatenate((arrydata, np.array([[arrydata[0, 0], arrydata[0, 1], medideptnigh, stdvdeptnigh, 0, 0, 0, 0]])), axis=0)
    ### spectrum of the host star, WASP-121
    gdat.meanwlenthomraww = arrymodl[:, 0]
    gdat.specstarthomraww = arrymodl[:, 9]
    
    # make a contour plot of geometric albedo without and with thermal component prior
    ## calculate the geometric albedo with the ATMO prior
    wlenmodl = arrymodl[:, 0]
    deptmodl = arrymodl[:, 1]
    indxwlenmodltess = np.where((wlenmodl > 0.6) & (wlenmodl < 0.95))[0]
    deptmodlther = np.mean(deptmodl[indxwlenmodltess])
    deptmodl = medideptseco + np.random.randn(numbsamp) * stdvdeptseco
    deptmodlrefl = np.random.randn(numbsamp) * 50 + 490 - deptmodlther
    
    print('HACKING')
    listalbg = deptmodlrefl * 1e-6 * (smajaxis / radiplan)**2
    listalbg = listalbg[listalbg > 0]
    
    prnt_list(listalbg, 'Albedo')
    prnt_list(listalbgtess, 'Albedo TESS only')
    prnt_list(listdeptseco, 'Secondary depth [ppm]')
    
    figr, axis = plt.subplots(figsize=gdat.figrsize[4, :])

    binsalbg = np.linspace(min(np.amin(listalbgtess), np.amin(listalbg)), max(np.amax(listalbgtess), np.amax(listalbg)), 100)
    meanalbg = (binsalbg[1:] + binsalbg[:-1]) / 2.
    pdfnalbgtess = scipy.stats.gaussian_kde(listalbgtess, bw_method=.2)(meanalbg)
    pdfnalbg = scipy.stats.gaussian_kde(listalbg, bw_method=.2)(meanalbg)
    #pdfnalbgtess = np.histogram(listalbgtess, bins=binsalbg)[0] / float(listalbgtess.size)
    #pdfnalbg = np.histogram(listalbg, bins=binsalbg)[0] / float(listalbg.size)
    axis.plot(meanalbg, pdfnalbgtess, label='TESS only', lw=2)
    axis.plot(meanalbg, pdfnalbg, label='TESS + ATMO prior', lw=2)
    axis.set_xlabel('$A_g$')
    axis.set_ylabel('$P(A_g)$')
    axis.legend()
    plt.subplots_adjust()
    path = pathimag + 'pdfn_albg.pdf'
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()
    
    
    # wavelength axis
    gdat.conswlentmpt = 0.0143877735e6 # [um K]
    
    ## read eclipse data
    liststrgfile = ['ContribFuncArr.txt', \
                    'EmissionDataArray.txt', \
                    #'RetrievalParamSamples.txt', \
                    'ContribFuncWav.txt', \
                    'EmissionModelArray.txt', \
                    'RetrievalPTSamples.txt', \
                    'pdependent_abundances/', \
                    ]
    
    
    ## get posterior on irradiation efficiency
    path = pathdata + 'ascii_output/RetrievalParamSamples.txt'
    listsampatmo = np.loadtxt(path)
    gdat.listpsii = listsampatmo[:, 2]
    histpsii, binspsii = np.histogram(gdat.listpsii)
    meanpsii = (binspsii[1:] + binspsii[:-1]) / 2.
    gdat.likeintp = scipy.interpolate.interp1d(meanpsii, np.log(histpsii), fill_value=-np.inf, bounds_error=False)
    
    ## sample from Bond albedo and circulation efficiency
    numbwalk = 10
    dictalbb = [gdat]
    indxwalk = np.arange(numbwalk)
    parainit = []
    for i in indxwalk:
        parainit.append(np.random.randn(2) * 0.001 + 0.02)
    objtsamp = emcee.EnsembleSampler(numbwalk, 2, retr_lpos_albb, args=dictalbb)
    objtsamp.run_mcmc(parainit, 10000)
    listalbb = objtsamp.flatchain[:, 0]
    listepsi = objtsamp.flatchain[:, 1]

    listsamp = np.empty((listalbb.size, 2))
    listsamp[:, 0] = listalbb
    listsamp[:, 1] = listepsi
    tdpy.mcmc.plot_grid(pathimag, 'post_albbespi', listsamp, ['$A_b$', r'$\varepsilon$'], plotsize=2.5)

    # plot ATMO posterior
    tdpy.mcmc.plot_grid(pathimag, 'post_atmo', listsampatmo, ['$\kappa_{IR}$', '$\gamma$', '$\psi$', \
                                                                        '[M/M$_{\odot}$]', '[C/C$_{\odot}$]', '[O/O$_{\odot}$]'], plotsize=2.5)
   

    # plot spectrum, depth, brightness temp
    path = pathdata + 'ascii_output/ContribFuncWav.txt'
    wlen = np.loadtxt(path)
    listcolr = ['k', 'm', 'purple', 'olive', 'olive', 'r', 'g']
    for i in range(15):
        listcolr.append('r')
    for i in range(28):
        listcolr.append('g')
   
    # infer brightness temperatures
    numbpara = 1
    numbsampwalk = 1000
    numbsampburnwalk = 100
    numbwalk = 5 * numbpara
    indxwalk = np.arange(numbwalk)
    numbsamp = numbsampwalk * numbwalk
    numbsampburn = numbsampburnwalk * numbwalk
    indxsampwalk = np.arange(numbsampwalk)
    gdat.numbdatatmpt = arrydata.shape[0]
    gdat.indxdatatmpt = np.arange(gdat.numbdatatmpt)
    tmpt = np.empty((numbsamp, gdat.numbdatatmpt))
    specarry = np.empty((2, 3, gdat.numbdatatmpt))
    for k in gdat.indxdatatmpt:
        
        gdat.minmwlen = arrydata[k, 0] - arrydata[k, 1]
        gdat.maxmwlen = arrydata[k, 0] + arrydata[k, 1]
        gdat.binswlen = np.linspace(gdat.minmwlen, gdat.maxmwlen, 100)
        gdat.meanwlen = (gdat.binswlen[1:] + gdat.binswlen[:-1]) / 2.
        gdat.diffwlen = (gdat.binswlen[1:] - gdat.binswlen[:-1]) / 2.
        gdat.cntrwlen = np.mean(gdat.meanwlen)
        if not (k == 0 or k == gdat.numbdatatmpt - 1):
            continue
        gdat.indxenerdata = k

        gdat.specstardiff = retr_modl_spec(gdat, 6396., strgtype='diff')
        gdat.specplandiff = retr_modl_spec(gdat, 2700., strgtype='diff')
        gdat.specstarlogt = retr_modl_spec(gdat, 6396., strgtype='logt')
        gdat.specplanlogt = retr_modl_spec(gdat, 2700., strgtype='logt')
        gdat.specstarintg = retr_modl_spec(gdat, 6396., strgtype='intg')
        gdat.specplanintg = retr_modl_spec(gdat, 2700., strgtype='intg')
        
        gdat.specstarthomlogt = scipy.interpolate.interp1d(gdat.meanwlenthomraww, gdat.specstarthomraww)(gdat.cntrwlen)
        gdat.specstarthomdiff = gdat.specstarthomlogt / gdat.cntrwlen
        gdat.specstarthomintg = np.sum(gdat.diffwlen * \
                                    scipy.interpolate.interp1d(gdat.meanwlenthomraww, gdat.specstarthomraww)(gdat.meanwlen) / gdat.meanwlen)

        gdat.deptobsd = arrydata[k, 2]
        gdat.stdvdeptobsd = arrydata[k, 3]
        gdat.varideptobsd = gdat.stdvdeptobsd**2
    
        numbdoff = gdat.numbdatatmpt - numbpara
        listlablpara = ['Temperature']
        liststrgpara = ['tmpt']
        indxpara = np.arange(numbpara)
        gdat.limtpara = np.empty((2, numbpara))
        # cons
        gdat.limtpara[0, 0] = 1000.
        gdat.limtpara[1, 0] = 4000.
        dictspec = [gdat]
        indxwalk = np.arange(numbwalk)
        parainit = []
        for i in indxwalk:
            parainit.append(np.empty(numbpara))
            meannorm = (gdat.limtpara[0, :] + gdat.limtpara[1, :]) / 2.
            stdvnorm = (gdat.limtpara[0, :] - gdat.limtpara[1, :]) / 10.
            parainit[i]  = (scipy.stats.truncnorm.rvs((gdat.limtpara[0, :] - meannorm) / stdvnorm, \
                                                        (gdat.limtpara[1, :] - meannorm) / stdvnorm)) * stdvnorm + meannorm
        if numbsampburnwalk == 0:
            raise Exception('')
        objtsamp = emcee.EnsembleSampler(numbwalk, numbpara, retr_lpos_spec, args=dictspec)
        parainitburn, prob, state = objtsamp.run_mcmc(parainit, numbsampburnwalk)
        objtsamp.reset()
        objtsamp.run_mcmc(parainitburn, numbsampwalk)
        tmpt[:, k] = objtsamp.flatchain[:, 0]

    posttmpttess = tmpt[:, np.array([0, gdat.numbdatatmpt - 1])]
    
    prnt_list(posttmpttess[:, 0], 'Dayside temperature')
    prnt_list(posttmpttess[:, 1], 'Nightside temperature')
    
    listtmptcont = (posttmpttess[:, 0] - posttmpttess[:, 1]) / posttmpttess[:, 0]
    prnt_list(listtmptcont, 'Temperature contrast')
    
    figr, axis = plt.subplots(4, 1, figsize=gdat.figrsize[3, :], sharex=True)
    ## stellar spectrum and TESS throughput
    axis[0].plot(arrymodl[:, 0], 1e-9 * arrymodl[:, 9], label='Host star, WASP-121', color='grey')
    axis[0].plot(0., 0., ls='--', label='TESS Throughput', color='grey')
    axis[0].set_ylabel(r'$\nu F_{\nu}$ [10$^9$ erg/s/cm$^2$]')
    axis[0].legend(fancybox=True, bbox_to_anchor=[0.8, 0.17, 0.2, 0.2])
    axistwin = axis[0].twinx()
    axistwin.plot(gdat.meanwlenband, gdat.thptband, color='grey', ls='--', label='TESS')
    axistwin.set_ylabel(r'Throughput')
    
    ## eclipse depths
    ### model
    axis[1].plot(arrymodl[:, 0], arrymodl[:,1], label='1D Retrieval (This work)')
    axis[1].plot(arrymodl[:, 0], arrymodl[:,2], label='Blackbody (This work)', alpha=0.3, color='skyblue')
    axis[1].fill_between(arrymodl[:, 0], arrymodl[:, 3], arrymodl[:, 4], alpha=0.3, color='skyblue')
    objtplotvivi, = axis[1].plot(wlenvivi, specvivi * 1e6, color='orange', alpha=0.6, lw=2)
    ### data
    for k in range(5):
        axis[1].errorbar(arrydata[k, 0], arrydata[k, 2], xerr=arrydata[k, 1], yerr=arrydata[k, 3], ls='', marker='o', color=listcolr[k])
    axis[1].errorbar(arrydata[5:22, 0], arrydata[5:22, 2], xerr=arrydata[5:22, 1], yerr=arrydata[5:22, 3], ls='', marker='o', color='r')
    axis[1].errorbar(arrydata[22:-1, 0], arrydata[22:-1, 2], xerr=arrydata[22:-1, 1], yerr=arrydata[22:-1, 3], ls='', marker='o', color='g')
    axis[1].set_ylabel(r'Depth [ppm]')
    axis[1].set_xticklabels([])
    
    ## spectra
    ### model
    objtplotretr, = axis[2].plot(arrymodl[:, 0], 1e-9 * arrymodl[:, 5], label='1D Retrieval (This work)', color='b')
    objtplotmblc, = axis[2].plot(arrymodl[:, 0], 1e-9 * arrymodl[:, 6], label='Blackbody (This work)', color='skyblue', alpha=0.3)
    objtploteblc = axis[2].fill_between(arrymodl[:, 0], 1e-9 * arrymodl[:, 7], 1e-9 * arrymodl[:, 8], color='skyblue', alpha=0.3)
    axis[2].legend([objtplotretr, (objtplotmblc, objtploteblc), objtplotvivi], \
                                                                ['1D Retrieval (This work)', 'Blackbody (This work)', 'GCM (Parmentier+2018)'], \
                                                                            bbox_to_anchor=[0.8, 1.3, 0.2, 0.2])
    ### data
    for k in range(5):
        axis[2].errorbar(arrydata[k, 0],  1e-9 * arrydata[k, 6], xerr=arrydata[k, 1], yerr=1e-9*arrydata[k, 7], ls='', marker='o', color=listcolr[k])
    axis[2].errorbar(arrydata[5:22, 0], 1e-9 * arrydata[5:22, 6], xerr=arrydata[5:22, 1], yerr=1e-9*arrydata[5:22, 7], ls='', marker='o', color='r')
    axis[2].errorbar(arrydata[22:-1, 0], 1e-9 * arrydata[22:-1, 6], xerr=arrydata[22:-1, 1], yerr=1e-9*arrydata[22:-1, 7], ls='', marker='o', color='g')
    

    axis[2].set_ylabel(r'$\nu F_{\nu}$ [10$^9$ erg/s/cm$^2$]')
    axis[2].set_xticklabels([])
    
    ## brightness temperature
    ### data
    for k in range(5):
        if k == 0:
            labl = 'TESS (This work)'
        if k == 1:
            labl = 'Z$^\prime$ (Delrez+2016)'
        if k == 2:
            labl = '$K_s$ (Kovacs\&Kovacs2019)'
        if k == 3:
            labl = 'IRAC $\mu$m (Garhart+2019)'
        #if k == 4:
        #    labl = 'IRAC 4.5 $\mu$m (Garhart+2019)'
        axis[3].errorbar(arrydata[k, 0], arrydata[k, 4], xerr=arrydata[k, 1], yerr=arrydata[k, 5], label=labl, ls='', marker='o', color=listcolr[k])
    axis[3].errorbar(arrydata[5:22, 0], arrydata[5:22, 4], xerr=arrydata[5:22, 1], yerr=arrydata[5:22, 5], label='HST G102 (Evans+2019)', ls='', marker='o', color='r')
    axis[3].errorbar(arrydata[22:-1, 0], arrydata[22:-1, 4], xerr=arrydata[22:-1, 1], yerr=arrydata[22:-1, 5], label='HST G141 (Evans+2017)', ls='', marker='o', color='g')
    #axis[3].errorbar(arrydata[:, 0], np.median(tmpt, 0), xerr=arrydata[:, 1], yerr=np.std(tmpt, 0), label='My calc', ls='', marker='o', color='c')
    axis[3].set_ylabel(r'$T_B$ [K]')
    axis[3].set_xlabel(r'$\lambda$ [$\mu$m]')
    axis[3].legend(fancybox=True, bbox_to_anchor=[0.8, 3.8, 0.2, 0.2], ncol=2)
    
    axis[1].set_ylim([20, None])
    axis[1].set_yscale('log')
    for i in range(4):
        axis[i].set_xscale('log')
    axis[3].set_xlim([0.5, 5])
    axis[3].xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    axis[3].xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.subplots_adjust(hspace=0., wspace=0.)
    path = pathimag + 'spec.pdf'
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()
    
    
    # get contribution function
    # interpolate the throughput
    path = pathdata + 'ascii_output/ContribFuncArr.txt'
    ctrb = np.loadtxt(path)
    presctrb = ctrb[0, :]
    path = pathdata + 'ascii_output/ContribFuncWav.txt'
    wlenctrb = np.loadtxt(path, skiprows=1)
    gdat.thptbandctrb = scipy.interpolate.interp1d(gdat.meanwlenband, gdat.thptband, fill_value=0, bounds_error=False)(wlenctrb)
    numbwlenctrb = wlenctrb.size
    indxwlenctrb = np.arange(numbwlenctrb)
    numbpresctrb = presctrb.size
    indxpresctrb = np.arange(numbpresctrb)
    ctrbtess = np.empty(numbpresctrb)
    for k in indxpresctrb:
        ctrbtess[k] = np.sum(ctrb[1:, k] * gdat.thptbandctrb)
    ctrbtess *= 1e-12 / np.amax(ctrbtess)

    # plot pressure-temperature, contribution
    path = pathdata + 'ascii_output/RetrievalPTSamples.txt'
    dataptem = np.loadtxt(path)
    
    liststrgcomp = ['CH4.txt', 'CO.txt', 'FeH.txt', 'H+.txt', 'H.txt', 'H2.txt', 'H2O.txt', 'H_.txt', 'He.txt', 'K+.txt', \
                                                        'K.txt', 'NH3.txt', 'Na+.txt', 'Na.txt', 'TiO.txt', 'VO.txt', 'e_.txt']
    listlablcomp = ['CH$_4$', 'CO', 'FeH', 'H$^+$', 'H', 'H$_2$', 'H$_2$O', 'H$^-$', 'He', 'K$^+$', \
                                                        'K', 'NH$_3$', 'Na$^+$', 'Na', 'TiO', 'VO', 'e$^-$']
    listdatacomp = []
    for strg in liststrgcomp:
        path = pathdata + 'ascii_output/pdependent_abundances/' + strg
        listdatacomp.append(np.loadtxt(path))
    
    ## contibution/PT/abun
    figr, axis = plt.subplots(nrows=1, ncols=2, sharey=True, gridspec_kw={'width_ratios': [1, 2]}, figsize=gdat.figrsize[0, :])
    ### pressure temperature
    numbsamp = dataptem.shape[0] - 1
    indxsamp = np.arange(numbsamp)
    for i in indxsamp[::100]:
        axis[0].plot(dataptem[i, :], dataptem[0, :], color='b', alpha=0.1)
    axis[0].plot(np.percentile(dataptem, 10, axis=0), dataptem[0, :], color='g')
    axis[0].plot(np.percentile(dataptem, 50, axis=0), dataptem[0, :], color='r')
    axis[0].plot(np.percentile(dataptem, 90, axis=0), dataptem[0, :], color='g')
    axis[0].set_xlim([1500, 3700])
    axis[0].set_xlabel('$T$ [K]')
    axis[0].set_yscale('log')
    axis[0].set_ylabel('$P$ [bar]')
    axis[0].invert_yaxis()
    axis[0].set_ylim([10., 1e-5])
    ### contribution function
    axistwin = axis[0].twiny()
    axistwin.fill(ctrbtess, presctrb, alpha=0.5, color='grey')
    axistwin.set_xticklabels([])
    ## abundance
    numbcomp = len(listdatacomp)
    indxcomp = np.arange(numbcomp)
    listobjtcolr = sns.color_palette("hls", numbcomp)
    axis[1].set_prop_cycle('color', listobjtcolr)
    listcolr = []
    for k in indxcomp:
        objt, = axis[1].plot(listdatacomp[k][:, 1], listdatacomp[k][:, 0])
        listcolr.append(objt.get_color())

    axis[1].xaxis.tick_top()
    
    arry = np.logspace(-16., 0., 21) # x 0.8
    for k in range(21):
        axis[1].axvline(arry[k], ls='--', alpha=0.1, color='k')
    arry = np.logspace(-5., 1., 11) # y 0.6
    for k in range(11):
        axis[1].axhline(arry[k], ls='--', alpha=0.1, color='k')
    listobjtcolr = sns.color_palette("hls", numbcomp)
    axis[1].set_prop_cycle('color', listobjtcolr)
    for k in indxcomp:
        if k == 0: # CH4
            xpos, ypos = 10**-12.8, 10**-2.3
        elif k == 1: # CO
            xpos, ypos = 10**-2.8, 10**-3.5
        elif k == 2: # FeH
            xpos, ypos = 10**-10.8, 10**-3.5
        elif k == 3: # H+
            xpos, ypos = 10**-12.8, 10**-4.1
        elif k == 4: # H
            xpos, ypos = 10**-1.6, 10**-2
        elif k == 5: # H2
            xpos, ypos = 10**-1.6, 10**-2.6
        elif k == 6: # H20
            xpos, ypos = 10**-8.8, 10**-4.1
        elif k == 7: # H_
            xpos, ypos = 10**-10., 10**0.4
        elif k == 8: # He
            xpos, ypos = 10**-1.6, 10**-4.1
        elif k == 9: # K+
            xpos, ypos = 10**-4.4, 10**-4.8
        elif k == 10: # K
            xpos, ypos = 10**-8.4, 10**-4.8
        elif k == 11: # Nh3
            xpos, ypos = 10**-13.6, 10**-4.1
        elif k == 12: # Na+
            xpos, ypos = 10**-4.4, 10**-3.8
        elif k == 13: # Na
            xpos, ypos = 10**-6, 10**-3.8
        elif k == 14: # TiO
            xpos, ypos = 10**-7.6, 10**-2
        elif k == 15: # VO
            xpos, ypos = 10**-6, 10**-2
        elif k == 16: # e-
            xpos, ypos = 10**-5.6, 10**-0.8
        else:
            xpos = 10**(np.random.rand() * 16. - 16.)
            ypos = 10**(np.random.rand() * 6. - 5.)
        axis[1].text(xpos, ypos, '%s' % listlablcomp[k], color=listcolr[k], size=10, va='center', ha='center')
        
    #liststrgcomp = ['CH4.txt', 'CO.txt', 'FeH.txt', 'H+.txt', 'H.txt', 'H2.txt', 'H2O.txt', 'H_.txt', 'He.txt', 'K+.txt', \
    #                                                    'K.txt', 'NH3.txt', 'Na+.txt', 'Na.txt', 'TiO.txt', 'VO.txt', 'e_.txt']
    #axis[1].set_yticklabels([])
    axis[1].set_xscale('log')
    axis[1].set_xlabel('Volume Mixing Ratio')
    axis[1].set_yscale('log')
    axis[1].set_xlim([1e-16, 1])
    #axis[1].invert_yaxis()
    plt.subplots_adjust(hspace=0., wspace=0., bottom=0.15)
    path = pathimag + 'ptem.pdf'
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()
  
    return


    # type of sampling
    samptype = 'emce'
    
    listaminchi2runs = []
    listsamppararuns = []
    listlistlablpara = []
    listlistlablparafull = []
        
    # get data
    path = pathdata + 'data_preparation/TESS.csv'
    arryfrst = np.loadtxt(path, delimiter=',', skiprows=1)
    
    # phase-fold the data
    arrythrd = np.copy(arryfrst)
    arrythrd[:, 0] = ((arryfrst[:, 0] - gdat.epoc) % gdat.peri) / gdat.peri
    
    # parse the data
    gdat.meanphas = arrythrd[:, 0]
    gdat.data = arrythrd[:, 1]
    gdat.datastdv = arrythrd[:, 2]
    
    # sort the data
    indxsort = np.argsort(gdat.meanphas)
    gdat.meanphas = gdat.meanphas[indxsort]
    gdat.data = gdat.data[indxsort]
    gdat.datastdv = gdat.datastdv[indxsort]

    # mask out the primary transit
    gdat.timetole = 6. / 24.
    indx = np.where(abs(gdat.meanphas - 0.5) < (1. - gdat.timetole / gdat.peri) / 2.)[0]
    if indx.size == gdat.meanphas.size:
        raise Exception('')
    gdat.data = gdat.data[indx]
    gdat.datastdv = gdat.datastdv[indx]
    gdat.meanphas = gdat.meanphas[indx]
    numbphas = gdat.data.size
    indxphas = np.arange(numbphas)
    
    fileoutp = open(pathdata + 'post.csv', 'w')
    fileoutp.write(' & ')
    fileoutp.write(' & ')
    
    # list of models
    listmodltype = ['simp', 'shft']
    
    fileoutp.write('\\\\\n')
    fileoutp.write('\\hline\n')
    fileoutp.write('$\chi^2_{\\nu}$ & ')
    fileoutp.write(' & ')
    fileoutp.write('\\hline\n')
    fileoutp.write('\\hline\n')
    fileoutp.write('\\\\\n')
    
    strgbins = 'unbd'
    strgmask = 'defa'
    for modltype in listmodltype:
        if modltype == 'simp':
            numbpara = 8
        if modltype == 'shft':
            numbpara = 9
        
        if samptype == 'emce':
            numbwalk = 5 * numbpara
            numbsampwalk = 1000
            numbsampburnwalk = numbsampwalk / 10
            if numbsampwalk == 0:
                raise Exception('')
            if numbsampburnwalk == 0:
                raise Exception('')
            numbsamp = numbsampwalk * numbwalk
            numbsampburn = numbsampburnwalk * numbwalk
        
        numbdoff = numbphas - numbpara
        
        listlablpara = ['Constant', 'Nightside emission [ppm]', \
                        'Planetary Modulation (B1 amplitude) [ppm]', 'B2 Amplitude [ppm]', 'B3 Amplitude [ppm]', \
                        'A1 Amplitude [ppm]', 'A2 Amplitude [ppm]', 'A3 Amplitude [ppm]']

        liststrgpara = ['cons', 'deptnigh', 'aco1', 'aco2', 'aco3', 'asi1', 'asi2', 'asi3']
        if modltype == 'shft':
            listlablpara += ['Phase offset [degree]']
            liststrgpara += ['shft']
        listlablparafull = listlablpara[:]
        liststrgparafull = liststrgpara[:]
        listlablparafull += ['Secondary transit depth [ppm]', 'Geometric albedo']
        liststrgparafull += ['deptseco', 'albg']
        
        numbparafull = len(liststrgparafull)
        indxpara = np.arange(numbpara)
        indxparafull = np.arange(numbparafull)
        gdat.limtpara = np.empty((2, numbpara))
        # cons
        gdat.limtpara[0, 0] = 0.9
        gdat.limtpara[1, 0] = 1.1
        # nightside emission [ppm]
        gdat.limtpara[0, 1] = 0.
        gdat.limtpara[1, 1] = 500.
        # amfo
        gdat.limtpara[0, 2] = -1e3
        gdat.limtpara[1, 2] = 0.
        # amfo
        gdat.limtpara[0, 3:8] = -1e3
        gdat.limtpara[1, 3:8] = 1e3
        if modltype == 'shft':
            # phas
            gdat.limtpara[0, 8] = -10.
            gdat.limtpara[1, 8] = 10.
        
        gdat.indxphasseco = np.where(abs(gdat.meanphas - 0.5) < dura / gdat.peri)[0]
        dictfour = [gdat, modltype]
        dicticdf = []
        
        if samptype == 'emce':
            indxwalk = np.arange(numbwalk)
            parainit = []
            for k in indxwalk:
                parainit.append(np.empty(numbpara))
                meannorm = (gdat.limtpara[0, :] + gdat.limtpara[1, :]) / 2.
                stdvnorm = (gdat.limtpara[0, :] - gdat.limtpara[1, :]) / 10.
                parainit[k]  = (scipy.stats.truncnorm.rvs((gdat.limtpara[0, :] - meannorm) / stdvnorm, \
                                                            (gdat.limtpara[1, :] - meannorm) / stdvnorm)) * stdvnorm + meannorm
            objtsamp = emcee.EnsembleSampler(numbwalk, numbpara, retr_lpos_four, args=dictfour)
            parainitburn, prob, state = objtsamp.run_mcmc(parainit, numbsampburnwalk)
            objtsamp.reset()
            objtsamp.run_mcmc(parainitburn, numbsampwalk)
            objtsave = objtsamp
        else:
        
            sampler = dynesty.NestedSampler(retr_llik_four, icdf, numbpara, logl_args=dictfour, ptform_args=dicticdf, bound='single', dlogz=1000.)
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
            postpara = objtsave.flatchain
        else:
            postpara = samppara
        
        sampllik = objtsave.lnprobability.flatten()
        aminchi2 = (-2. * np.amax(sampllik) / numbdoff)
        listaminchi2runs.append(aminchi2)
        
        postparatemp = np.copy(postpara)
        postpara = np.empty((numbsamp, numbparafull))
        postpara[:, :-(numbparafull - numbpara)] = postparatemp
        for k in indxsamp:
            # secondary depth
            postpara[k, -2] = postpara[k, 2] - postpara[k, 3]
            # geometric albedo
            postpara[:, -1] = -(smajaxis**2 / radiplan**2) * (1e-6 * postpara[:, 2])
    
        listsamppararuns.append(postpara)
        listlistlablpara.append(listlablpara)
        listlistlablparafull.append(listlablparafull)
    
        # plot the posterior
        
        ### histogram
        for k in indxparafull:
            figr, axis = plt.subplots()
            if samptype == 'emce':
                axis.hist(postpara[:, k]) 
            else:
                axis.hist(samppara[:, k]) 
            axis.set_xlabel(listlablparafull[k])
            path = pathimag + 'diag/hist_%s_%s_%s_%s.pdf' % (liststrgparafull[k], modltype, strgmask, strgbins)
            plt.tight_layout()
            print('Writing to %s...' % path)
            plt.savefig(path)
            plt.close()
        
        if samptype == 'nest':
            for keys in objtsave:
                if isinstance(objtsave[keys], np.ndarray) and objtsave[keys].size == numbsamp:
                    figr, axis = plt.subplots()
                    axis.plot(indxsamp, objtsave[keys])
                    path = pathimag + '%s_%s.pdf' % (keys, modltype)
                    print('Writing to %s...' % path)
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
            print('Writing to %s...' % path)
            plt.savefig(path)
            plt.close()
        
            chi2 = -2. * objtsave.lnprobability
            print('Posterior-mean chi2: ')
            print(np.mean(chi2))
            print('Posterior-mean chi2 per dof: ')
            print(np.mean(chi2) / numbdoff)
        
        tdpy.mcmc.plot_grid(pathimag, 'post_sinu_%s_%s_%s' % (modltype, strgmask, strgbins), postpara, listlablparafull, plotsize=2.5)
    
        ### sample model phas
        numbsampplot = 100
        indxsampplot = np.random.choice(indxsamp, numbsampplot, replace=False)
        yerr = gdat.datastdv
        
        numbphasfine = 1000
        gdat.meanphasfine = np.linspace(0.1, 0.9, numbphasfine)
        phasmodlfine = np.empty((numbsampplot, numbphasfine))
        indxphasfineseco = np.where(abs(gdat.meanphasfine - 0.5) < dura / gdat.peri)[0]
        for k, indxsampplottemp in enumerate(indxsampplot):
            if samptype == 'emce':
                objttemp = objtsave.flatchain
            else:
                objttemp = samppara
            cons = objttemp[indxsampplottemp, 0]
            deptnigh = objttemp[indxsampplottemp, 1]
            amfo = objttemp[indxsampplottemp, 2:8]
            if modltype == 'shft':
                shft = objttemp[indxsampplottemp, 8]
            else:
                shft = 0.
            phasmodlfine[k, :], deptpmodfine = retr_modl_four(gdat, gdat.meanphasfine, indxphasfineseco, cons, deptnigh, amfo, shft, modltype)
        
        strgextn = '%s_%s_%s' % (modltype, strgmask, strgbins)
        ## log-likelihood
        figr, axis = plt.subplots(figsize=(12, 6))
        axis.errorbar(gdat.phasbind % 1., (gdat.pcurdetrbind - 1) * 1e6, yerr=1e6*gdat.stdvpcurdetrbind, color='k', marker='o', ls='', markersize=1)
        for k, indxsampplottemp in enumerate(indxsampplot):
            axis.plot(gdat.meanphasfine, (phasmodlfine[k, :] - 1) * 1e6, alpha=0.5, color='b')
        axis.set_xlim([0.1, 0.9])
        axis.set_ylim([-400, 1000])
        axis.set_ylabel('Relative Flux - 1 [ppm]')
        axis.set_xlabel('Phase')
        plt.tight_layout()
        path = pathimag + 'pcur_sine_%s.pdf' % strgextn
        print('Writing to %s...' % path)
        plt.savefig(path)
        plt.close()

        for n in indxparafull:
            fileoutp.write('%s & ' % listlablparafull[n])
            ydat = np.median(postpara[:, n])
            uerr = np.percentile(postpara[:, n], 84.) - ydat
            lerr = ydat - np.percentile(postpara[:, n], 16.)
            fileoutp.write('$%.3g \substack{+%.3g \\\\ -%.3g}$' % (ydat, uerr, lerr))
            fileoutp.write(' & ')
            fileoutp.write('\\\\\n')
    fileoutp.write('\\hline\n')
    fileoutp.close()

init()
