# plotting
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
mpl.rc('image', interpolation='none', origin='lower')

import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)

# numpy
import numpy as np
from numpy import *
from numpy.random import *
from numpy.random import choice

# scipy
import scipy as sp
from scipy import ndimage
from scipy.interpolate import *
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss
from scipy import ndimage

# multiprocessing
import multiprocessing as mp

# healpy
import healpy as hp
from healpy.rotator import angdist
from healpy import ang2pix

# pyfits
import pyfits as pf

# utilities
import os, time, sys, datetime, warnings, getpass, glob, inspect

# tdpy
import tdpy.util
import tdpy.mcmc

# profiling
import cProfile, pstats, cPickle


def make_maps():
    
    bind = False
    
    cmnd = 'mkdir -p $FERMI_DATA/exposure/ferm_line'
    os.system(cmnd)
    
    liststrgener = ['ENERGY%d' % (k + 1) for k in indxener]
    liststrgchan = ['CHANNEL%d' % (k + 1) for k in indxener]
    
    numbproc = 20
    minmweek = 11
    maxmweek = 411
    weekpart = (maxmweek - minmweek) / numbproc
    
    listweek = []
    for k in range(numbproc):
        listweek.append(arange(minmweek + k * weekpart, minmweek + (k + 1) * weekpart))
        
    # process pool
    pool = mp.Pool(numbproc)

    # spawn the processes
    pool.map(make_maps_sing, range(numbproc))

    pool.close()
    pool.join()


def make_maps_sing(indxproc):

    for t in listweek[indxproc]:
        
        print inspect.stack()[0][3] + ', proc%d is working on week %d...' % (indxproc, t)
        
        for m in indxevtt:
            
            evnt = os.environ["FERMI_DATA"] + 'weekly/photon/lat_photon_weekly_w%03d_p302_v001.fits' % t
            spac = os.environ["FERMI_DATA"] + 'weekly/spacecraft/lat_spacecraft_weekly_w%03d_p202_v001.fits' % t
            sele = os.environ["FERMI_DATA"] + 'exposure/ferm_line/sele_pass8_evtc128_evtt%03d_week%03d.fits' % (evtt[m], t)
            filt = os.environ["FERMI_DATA"] + 'exposure/ferm_line/filt_pass8_evtc128_evtt%03d_week%03d.fits' % (evtt[m], t)
            live = os.environ["FERMI_DATA"] + 'exposure/ferm_line/live_pass8_evtc128_evtt%03d_week%03d.fits' % (evtt[m], t)
            cnts = os.environ["FERMI_DATA"] + 'exposure/ferm_line/cnts_pass8_evtc128_evtt%03d_week%03d.fits' % (evtt[m], t)
            expo = os.environ["FERMI_DATA"] + 'exposure/ferm_line/expo_pass8_evtc128_evtt%03d_week%03d.fits' % (evtt[m], t)

            cmnd = 'gtselect infile=' + evnt + ' outfile=' + sele + ' ra=INDEF dec=INDEF rad=INDEF tmin=INDEF tmax=INDEF' + \
                        ' emin=%d emax=%d zmax=90 evclass=%d evtype=%d' % (minmener * 1e3, maxmener * 1e3, evtc, evtt[m])
            os.system(cmnd)

            cmnd = 'gtmktime evfile=' + sele + ' scfile=' + spac + ' filter="DATA_QUAL==1 && LAT_CONFIG==1"' + ' outfile=' + filt + ' roicut=no'

            os.system(cmnd)
            
            if bind or (t == 0 and m == 0):
                cntstemp = os.environ["FERMI_DATA"] + 'exposure/ferm_line/cnts.fits'
                cmnd = 'gtbin evfile=' + filt + ' scfile=' + spac + ' outfile=' + cntstemp + \
                            ' ebinalg=LOG emin=%d emax=%d enumbins=%d algorithm=HEALPIX' % (minmener * 1e3, maxmener * 1e3, numbener) + \
                            ' hpx_ordering_scheme=RING coordsys=GAL hpx_order=8 hpx_ebin=yes'
                os.system(cmnd)
            else:
                cntstemp = cnts

            cmnd = 'gtltcube evfile=' + filt + ' scfile=' + spac + ' outfile=' + live + ' dcostheta=0.025 binsz=1'
            os.system(cmnd)

            cmnd = 'gtexpcube2 infile=' + live + ' outfile=' + expo + ' cmap=' + cntstemp + ' irfs=CALDB evtype=%03d bincalc=CENTER' % evtt[m]
            os.system(cmnd)


def prep_maps():
    

    liststrgener = ['ENERGY%d' % (k + 1) for k in indxener]
    liststrgchan = ['CHANNEL%d' % (k + 1) for k in indxener]
    
    minmweek = 11.
    maxmweek = 411.
    ntime = 4
        
    cntsweek = zeros((numbener, numbpixl, numbevtt, ntime))
    expoweek = zeros((numbener, numbpixl, numbevtt, ntime))
    cnts = zeros((numbener, numbpixl, numbevtt))
    expo = zeros((numbener, numbpixl, numbevtt))

    for w in arange(minmweek, maxmweek):

        print inspect.stack()[0][3] + ' is working on week %d...' % w 

        for m in indxevtt:
            path = os.environ["FERMI_DATA"] + 'exposure/ferm_line/expo_pass%d_evtc%03d_evtt' % (reco, evtc)
            path += '%03d_week%03d.fits' % (evtt[m], w)
            
            expoarry = pf.getdata(path, 1)
            for i in indxener:
                expo[i, :, m] = expoarry[liststrgener[i]]

        for m in indxevtt:
            path = os.environ["FERMI_DATA"] + 'exposure/ferm_line/cnts_pass%d_evtc%03d_evtt' % (reco, evtc)
            path += '%03d_week%03d.fits' % (evtt[m], w)
            cntsarry = pf.getdata(path)
            for i in indxener:
                cnts[i, :, m] = cntsarry[liststrgchan[i]]

        t = int(floor(ntime * (w - minmweek) / (maxmweek - minmweek)))
        expoweek[:, :, :, t] += expo
        cntsweek[:, :, :, t] += cnts

    fluxweek = zeros((numbener, numbpixl, numbevtt, ntime))
    indxexpoweek = where(expoweek > 0.) 
    fluxweek[indxexpoweek] = cntsweek[indxexpoweek] / expoweek[indxexpoweek] / apix
    fluxweek /= diffener[:, None, None, None]
    
    expo = sum(expoweek, 3)
    cnts = sum(cntsweek, 3)
    
    flux = zeros((numbener, numbpixl, numbevtt))
    indxexpo = where(expo > 0.)
    flux[indxexpo] = cnts[indxexpo] / expo[indxexpo] / apix
    flux /= diffener[:, None, None]
    
    #for t in range(ntime):
    #    path = os.environ["FERM_LINE_DATA_PATH"] + '/flux_time%d' % t + '.fits'
    #    pf.writeto(path, fluxweek[:, :, :, t], clobber=True)
    #    path = os.environ["FERM_LINE_DATA_PATH"] + '/expo_time%d' % t + '.fits'
    #    pf.writeto(path, expoweek[:, :, :, t], clobber=True)

    path = os.environ["FERM_LINE_DATA_PATH"] + '/flux.fits'
    pf.writeto(path, flux, clobber=True)
    path = os.environ["FERM_LINE_DATA_PATH"] + '/expo.fits'
    pf.writeto(path, expo, clobber=True)


def plot_datacntsmean(gdat):

    figr, axcl = plt.subplots(gdat.numbevtt, 1, figsize=(10, 25))
    for m, axis in enumerate(axcl):
        xdat = gdat.meanenerwndw
        ydat = gdat.datacntsmeanwndw[:, m]
        yerr = sqrt(gdat.datacntsmeanwndw[:, m])
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', elinewidth=2, capsize=20)
        axis.set_xlabel('E [GeV]')
        axis.set_xlim([gdat.minmenerwndw, gdat.maxmenerwndw])
        axis.set_title('EDISP%d' % m)
        axis.set_xscale('log')
        axis.set_yscale('log')
    figr.text(0.05, 0.5, '$N$', ha='center', va='center', rotation=90)
    plt.subplots_adjust(hspace=0.3)
    plt.savefig(gdat.pathplot + 'datameancnts_%.3g.png' % gdat.enercntrwndw)
    plt.close(figr)
          

def retr_datapara(gdat):

    if gdat.modltype == 'nullmod0':
        gdat.numbpara = 2
    if gdat.modltype == 'altrmod0':
        gdat.numbpara = 2 + gdat.numbalmp
        
    if gdat.modltype == 'nullmod1':
        gdat.numbpara = 2 + gdat.numbalmp
    if gdat.modltype == 'altrmod1':
        gdat.numbpara = 2 + 2 * gdat.numbalmp

    dictpara = dict()
    minmpara = zeros(gdat.numbpara)
    maxmpara = zeros(gdat.numbpara)
    namepara = empty(gdat.numbpara, dtype=object)
    scalpara = empty(gdat.numbpara, dtype=object)
    lablpara = empty(gdat.numbpara, dtype=object)
    unitpara = empty(gdat.numbpara, dtype=object)
    varipara = zeros(gdat.numbpara)
    
    dictpara['slopback'] = 0
    namepara[0] = 'slopback'
    minmpara[0] = 1.5
    maxmpara[0] = 4.
    scalpara[0] = 'self'
    lablpara[0] = r'$\alpha_b$'
    unitpara[0] = ''
    varipara[0] = 1e-1

    cntrpara = 1
    if gdat.modltype == 'altrmod0' or gdat.modltype == 'nullmod0' or gdat.modltype == 'altrmod1' or gdat.modltype == 'nullmod1':
            
        dictpara['normback'] = 1
        namepara[1] = 'normback'
        minmpara[1] = 1e-3
        maxmpara[1] = 1e0
        scalpara[1] = 'logt'
        lablpara[1] = r'$A_b$'
        unitpara[1] = ''
        varipara[1] = 1e-1

        cntrpara += 1
    
    if gdat.modltype == 'nullmod1' or gdat.modltype == 'altrmod1':
        defn_almcpara('powrspec', cntrpara, dictpara, namepara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara)
        cntrpara += gdat.numbalmp
        
    if gdat.modltype == 'altrmod0' or gdat.modltype == 'altrmod1':
        defn_almcpara('linespec', cntrpara, dictpara, namepara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara)
        
    strgpara = lablpara + ' ' + unitpara
    datapara = namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara
    
    return datapara


def defn_almcpara(typepara, cntrpara, dictpara, namepara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara):
    
    listsphl, listsphm = hp.Alm.getlm(gdat.maxmsphl)
    listsphl = tile(listsphl, 2)
    listsphm = tile(listsphm, 2)

    temp = []
    for k in range(gdat.numbalmc):
        if listsphm[gdat.numbalmc+k] == 0:
            temp.append(gdat.numbalmc + k)
    listsphl = delete(listsphl, temp)
    listsphm = delete(listsphm, temp)

    maxmindxpara = cntrpara + gdat.numbalmp
        
    for k in range(cntrpara, maxmindxpara):
        dictpara['para%04d' % k] = k
        namepara[k] = 'para%04d' % k
        minmpara[k] = -100.
        maxmpara[k] = 100.
        scalpara[k] = 'self'
        if typepara == 'powrspec':
            strg = '$a^{P}'
        if typepara == 'linespec':
            strg = '$a^{L}'
        if k < cntrpara + gdat.numbalmc:
            lablpara[k] = 'Re{'
        else:
            lablpara[k] = 'Im{' 
        lablpara[k] += '%s_{%d%d}$' % (strg, listsphl[k - cntrpara],   listsphm[k - cntrpara])
        if k < gdat.numbalmc:
            lablpara[k] += '}'
        else:
            lablpara[k] += '}'
        unitpara[k] = ''
        varipara[k] = 1e-5


def retr_aexp(scaldevi, stdv, skew, bias, slop):

    aexp = slop / stdv / sp.special.gamma(1. / slop) * skew / (1. + skew**2)
    indxscalfrwd = where(scaldevi - bias > 0.)
    indxscalback = where(scaldevi - bias <= 0.)
    aexp[indxscalfrwd] *= exp(-(skew[indxscalfrwd] / stdv[indxscalfrwd])**slop[indxscalfrwd] * (scaldevi[indxscalfrwd] - bias[indxscalfrwd])**slop[indxscalfrwd])
    aexp[indxscalback] *= exp(-(1. / skew[indxscalback] / stdv[indxscalback])**slop[indxscalback] * (bias[indxscalback] - scaldevi[indxscalback])**slop[indxscalback])

    return aexp


def retr_fermedfn(gdat, thisener, enercntr):
  
    enerdevi = thisener - enercntr

    thisnumbener = thisener.size

    path = os.environ["PCAT_DATA_PATH"] + '/irfn/edisp_P8R2_SOURCE_V6_EDISP.fits'

    strgformpara = ['F', 'S1', 'S2', 'K1', 'K2', 'BIAS', 'BIAS2', 'PINDEX1', 'PINDEX2']

    listhdun = pf.open(path)

    numbctht = 8
    numbformpara = 9
    numbscalpara = 6
    frac = zeros((thisnumbener, numbctht, gdat.numbevtt))
    stdv = zeros((2, thisnumbener, numbctht, gdat.numbevtt))
    skew = zeros((2, thisnumbener, numbctht, gdat.numbevtt))
    bias = zeros((2, thisnumbener, numbctht, gdat.numbevtt))
    slop = zeros((2, thisnumbener, numbctht, gdat.numbevtt))
    scalpara = zeros((gdat.numbevtt, numbscalpara))
    for m in gdat.indxevtt:
        
        # form parameters
        ## get HDU
        hdun = listhdun['ENERGY DISPERSION_EDISP%d' % m]

        ## get energy and angle axes
        if m == 0:
            minmcthtirfn = hdun.data['CTHETA_LO'].flatten()
            maxmcthtirfn = hdun.data['CTHETA_HI'].flatten()
            cthtirfn = sqrt(minmcthtirfn * maxmcthtirfn)
            minmenerirfn = hdun.data['ENERG_LO'].flatten() * 1e-3 # [GeV]
            maxmenerirfn = hdun.data['ENERG_HI'].flatten() * 1e-3 # [GeV]
            meanenerirfn = sqrt(minmenerirfn * maxmenerirfn)

        ## get form parameters
        for k in range(numbformpara):
            formparatemp = interp1d(meanenerirfn, hdun.data[strgformpara[k]].squeeze(), axis=1)(enercntr)
            if k == 0:
                frac[:, :, m] = formparatemp[None, :]
            if k == 1:
                stdv[0, :, :, m] = formparatemp[None, :]
            if k == 2:
                stdv[1, :, :, m] = formparatemp[None, :]
            if k == 3:
                skew[0, :, :, m] = formparatemp[None, :]
            if k == 4:
                skew[1, :, :, m] = formparatemp[None, :]
            if k == 5:
                bias[0, :, :, m] = formparatemp[None, :]
            if k == 6:
                bias[1, :, :, m] = formparatemp[None, :]
            if k == 7:
                slop[0, :, :, m] = formparatemp[None, :]
            if k == 8:
                slop[1, :, :, m] = formparatemp[None, :]

        # scale parameters
        scalpara[m, :] = listhdun['EDISP_SCALING_PARAMS_EDISP%d' % m].data['EDISPSCALE']

    if False: 
        print 'hey'
        print mean(frac, axis=1)
        for k in range(2):
            print 'stdv'
            print mean(stdv, axis=1)
            print 'skew'
            print mean(skew, axis=1)
            print 'bias'
            print mean(bias, axis=1)
            print 'slop'
            print mean(slop, axis=1)
            print
        print 
        print
    scalfact = scalpara[None, None, :, 0] * log(thisener[:, None, None])**2 + \
               scalpara[None, None, :, 1] * cthtirfn[None, :, None]**2 + \
               scalpara[None, None, :, 2] * log(thisener[:, None, None]) + \
               scalpara[None, None, :, 3] * cthtirfn[None, :, None] + \
               scalpara[None, None, :, 4] * log(thisener[:, None, None]) * cthtirfn[None, :, None] + \
               scalpara[None, None, :, 5]
            
    gdat.scaldevi = enerdevi[:, None, None] / enercntr / scalfact
    
    edfncom0 = frac[:, :, :] * retr_aexp(gdat.scaldevi, stdv[0, :, :, :], skew[0, :, :, :], bias[0, :, :, :], slop[0, :, :, :])
    edfncom1 = (1. - frac[:, :, :]) * retr_aexp(gdat.scaldevi, stdv[1, :, :, :], skew[1, :, :, :], bias[1, :, :, :], slop[1, :, :, :])
    edfn = edfncom0 + edfncom1
    
    edfncom0 = mean(edfncom0, axis=1)
    edfncom1 = mean(edfncom1, axis=1)
    edfn = mean(edfn, axis=1)

    return edfn, edfncom0, edfncom1


def plot_edfnwndw(gdat):

    listlinestyl = ['-', '--', '--']
    listcolr = ['b', 'g', 'r', 'm']
    figr, axis = plt.subplots()
    for m in gdat.indxevtt:
        axis.plot(gdat.meanenerwndw, gdat.edfnwndw[:, m], c=listcolr[m], ls=listlinestyl[0])
        axis.plot(gdat.meanenerwndw, gdat.edfncom0wndw[:, m], c=listcolr[m], ls=listlinestyl[1])
        axis.plot(gdat.meanenerwndw, gdat.edfncom1wndw[:, m], c=listcolr[m], ls=listlinestyl[2])
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_yscale('log')
    plt.xlim([gdat.binsenerwndw[0], gdat.binsenerwndw[-1]])
    plt.savefig(gdat.pathplot + 'edfnwndw_%s.png' % gdat.strgenercntrwndw)
    plt.close(figr) 

    figr, axis = plt.subplots()
    for m in gdat.indxevtt:
        axis.plot(mean(gdat.scaldevi[:, :, m], 1), gdat.edfnwndw[:, m], c=listcolr[m], ls=listlinestyl[0])
        axis.plot(mean(gdat.scaldevi[:, :, m], 1), gdat.edfncom0wndw[:, m], c=listcolr[m], ls=listlinestyl[1])
        axis.plot(mean(gdat.scaldevi[:, :, m], 1), gdat.edfncom1wndw[:, m], c=listcolr[m], ls=listlinestyl[2])
    axis.set_xlabel(r'$x$')
    plt.xlim([-15., 15.])
    plt.ylim([1e-6, 1e0])
    axis.set_yscale('log')
    plt.savefig(gdat.pathplot + 'edfnwndwscal_%s.png' % gdat.strgenercntrwndw)
    plt.close(figr) 


def writ_maps():

    fdfm = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + '/fdfm.fits')
    dataflux = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + '/flux.fits')
    expo = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + '/expo.fits')
        
    listnumbside = [1, 2, 4, 8, 16, 32, 64]
    for numbside in listnumbside:
        
        numbpixl = numbside**2 * 12
        
        fdfmtemp = zeros((numbener, numbpixl, numbevtt))
        datafluxtemp = zeros((numbener, numbpixl, numbevtt))
        expotemp = zeros((numbener, numbpixl, numbevtt))
        
        fdfmtemp = hp.ud_grade(fdfm, numbside)
        for i in indxener:
            for m in indxevtt:
                datafluxtemp[i, :, m] = hp.ud_grade(dataflux[i, :, m], numbside)
                expotemp[i, :, m] = hp.ud_grade(expo[i, :, m], numbside)

        pathfdfm = os.environ["FERM_LINE_DATA_PATH"] + '/fdfm%03d.fits' % numbside
        pathflux = os.environ["FERM_LINE_DATA_PATH"] + '/flux%03d.fits' % numbside
        pathexpo = os.environ["FERM_LINE_DATA_PATH"] + '/expo%03d.fits' % numbside
        pf.writeto(pathfdfm, fdfmtemp, clobber=True)
        pf.writeto(pathflux, datafluxtemp, clobber=True)
        pf.writeto(pathexpo, expotemp, clobber=True)


def retr_linecnts(thisener, enercntr):
    
    linecnts = 1. / sqrt(2. * pi * stdvdisp[None, :]**2) * exp(-0.5 * (thisener[:, None] - enercntr)**2 / stdvdisp[None, :]**2)
    
    return linecnts


def retr_modlcnts(sampvarb):
    
    slopback = sampvarb[0]
    normback = sampvarb[1]
   
    enertemppowr = (gdat.meanenerwndw / gdat.enercntrwndw)**(-slopback)
    gdat.modlcnts[:] = normback * gdat.fdfmfluxwndw * gdat.expowndw * gdat.diffenerwndw[:, None, None] * gdat.apix * enertemppowr[:, None, None]
    cntr = 2
        
    if gdat.modltype == 'nullmod0':
        numbiter = 0
    elif gdat.modltype == 'altrmod0' or gdat.modltype == 'nullmod1':
        numbiter = 1
    else:
        numbiter = 2
        
    for k in range(numbiter):
        
        almcreal = sampvarb[cntr:cntr+gdat.numbalmc]
        gdat.almcimag[gdat.maxmsphl+1:] = sampvarb[cntr+gdat.numbalmc:cntr+gdat.numbalmp]
        almc = almcreal + gdat.almcimag * 1j
        modlcntstemp = hp.alm2map(almc, gdat.numbside, verbose=False)
        cntr += gdat.numbalmp
        if k == 1 or gdat.modltype == 'altrmod0':
            gdat.modlcnts[:] += modlcntstemp[None, :, None] * gdat.edfnwndw[:, None, :]
        else:
            gdat.modlcnts[:] += modlcntstemp[None, :, None] * enertemppowr[:, None, None]
       
        if gdat.verbtype > 1:
            print 'Iter %d' % k
            print 'almc'
            print almc
    if gdat.verbtype > 1:
        print 'sampvarb'
        print sampvarb
        print 'modltype'
        print gdat.modltype
        print 'modlcnts'
        print amin(gdat.modlcnts)
        print amax(gdat.modlcnts)
        print mean(gdat.modlcnts)
        print

    return gdat.modlcnts


def retr_llik(sampvarb, init=False):
    
    modlcntstemp = retr_modlcnts(sampvarb)

    if gdat.lliktype == 'bind':
        llik = sum(gdat.datacntswndw * log(modlcntstemp) - modlcntstemp)
    else:
        llik = sum(exp(-modlcntsintp))
    
    return llik, gdat.sampcalc
    

def plot_cnts(thiscnts, strg):

    cart = tdpy.util.retr_cart(sum(thiscnts[gdat.indxenerwndwplot, :, :], 1))

    figr, axis = plt.subplots()
    
    imag = plt.imshow(cart, origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)

    plt.savefig(gdat.pathplot + '%scnts%s.png' % (strg, gdat.rtag))
    plt.close(figr)


def retr_numbalmc(maxmsphl):
    
    numbalmc = maxmsphl * (maxmsphl + 1) / 2 + maxmsphl + 1
    
    return numbalmc

def retr_numbalmp(maxmsphl, numbalmc):
    
    numbalmp = numbalmc * 2 - maxmsphl - 1
    
    return numbalmp


def plot_sphl(gdat):

    listmaxmsphl = logspace(0., 2., 20)
    listnumbalmc = retr_numbalmc(listmaxmsphl)
    listnumbalmp = retr_numbalmp(listmaxmsphl, listnumbalmc)
    listangl = 180. * sqrt(2. / pi / listmaxmsphl)
    
    figr, axis = plt.subplots()
    
    axis.loglog(listmaxmsphl, listnumbalmc, label='$N_{almc}$')
    axis.loglog(listmaxmsphl, listnumbalmp, label='$N_{almp}$')
    axistwin = axis.twinx()
    axistwin.loglog(listmaxmsphl, listangl, label=r'$\theta$', color='r')

    axis.legend(loc=2)
    axistwin.legend(loc=4)
    
    axis.set_ylabel('$N$')
    axistwin.set_ylabel(r'$\theta$ [degree]')
    axis.set_xlabel('$l$')
    axis.set_title('')
    plt.subplots_adjust(hspace=0.1)
    plt.savefig(gdat.pathplot + 'sphl.png')
   

class globdatastrt(object):
    
    def __init__(self):
        pass


def init( \
         methtype='almc', \
         datatype='mock', \
        ):

    # initialize the global object 
    gdat = globdatastrt()
    
    # event class
    #gdat.reco = 8
    #gdat.evtc = 128
    
    # energy axis
    gdat.minmener = 30. # [GeV]
    gdat.maxmener = 300. # [GeV]
    gdat.resoener = 0.1
    gdat.numbener = int(4. * (gdat.maxmener / gdat.minmener) / gdat.resoener)
    gdat.binsener = logspace(log10(gdat.minmener), log10(gdat.maxmener), gdat.numbener + 1)
    gdat.meanener = sqrt(gdat.binsener[1:] * gdat.binsener[0:-1])
    gdat.diffener = gdat.binsener[1:] - gdat.binsener[0:-1]
    gdat.indxener = arange(gdat.numbener)
    #gdat.indxenerpivt = gdat.numbener / 2
    #gdat.enerpivt = gdat.meanener[gdat.indxenerpivt]

    # event type axis
    gdat.evtt = array([64, 128, 256, 512])
    gdat.numbevtt = gdat.evtt.size
    gdat.indxevtt = arange(gdat.numbevtt)
    
    # healpix setup
    gdat.numbside = 8
    gdat.numbpixl = gdat.numbside**2 * 12
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.apix = 4. * pi / gdat.numbpixl
    
    # base project path
    gdat.pathbase = os.environ["FERM_LINE_DATA_PATH"]
    
    # plot setup
    ## ROI settings
    gdat.exttrofi = [-180., 180., -90., 90.]
    ## prepare plot folder
    gdat.pathplot = gdat.pathbase + '/png/'
    cmnd = 'mkdir -p ' + gdat.pathplot
    os.system(cmnd)
    if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
        cmnd = 'mv ' + gdat.pathplot + '/* /n/pan/www/tansu/png/ferm_line/'
        os.system(cmnd)
    
    # maximum multipole moment
    gdat.maxmsphl = 1

    # number of alm coefficients
    gdat.numbalmc = retr_numbalmc(gdat.maxmsphl)

    # number of parameters
    gdat.numbalmp = retr_numbalmp(gdat.maxmsphl, gdat.numbalmc)
    
    # initial plots
    ## spherical harmonics
    plot_sphl(gdat)
    
    # auxiliary variables to be collected during sampling
    gdat.sampcalc = []
   
    gdat.lliktype = 'bind'

    # energy window
    ## number of energy bins in the sliding window
    gdat.numbenerwndw = 10 * int(0.03 * gdat.numbener) + 1
    gdat.indxenerwndwplot = gdat.numbenerwndw / 2
    ## hald of the sliding window
    gdat.numbenerwndwside = (gdat.numbenerwndw - 1) / 2
    
    # exposure
    
    # temp
    strgexpo = '/expo%03d.fits' % gdat.numbside   
    gdat.expo = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + strgexpo)
    expotemp = copy(gdat.expo)
    gdat.expo = empty((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    # temp
    gdat.expo[:] = 1000. * expotemp[0, :, :][None, :, :]

    # diffuse model
    fdfmfluxtemp = tdpy.util.retr_fdfm(gdat.binsener, numbside=gdat.numbside)  
    gdat.fdfmflux = empty((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for m in gdat.indxevtt:
        gdat.fdfmflux[:, :, m] = fdfmfluxtemp
    
    # temp
    #path = os.environ["FERM_LINE_DATA_PATH"] + '/fdfm%03d.fits' % gdat.numbside
    #gdat.fdfmflux = pf.getdata(path)
    
    gdat.fdfmcnts = gdat.fdfmflux * gdat.expo * gdat.diffener[:, None, None] * gdat.apix
    
    # common MCMC settings 
    gdat.verbtype = 1
    gdat.factthin = 1
    gdat.numbproc = 1
    gdat.optiprop = True
    
    # convenience data structure
    ## alm coefficients
    gdat.almcimag = zeros(gdat.numbalmc)
    # model counts
    gdat.modlcnts = empty((gdat.numbenerwndw, gdat.numbpixl, gdat.numbevtt))
    
    # energy axis label
    gdat.strgenercntr = ['%.3g' % gdat.meanener[i] for i in gdat.indxener]

    # center energies
    gdat.numbenercntr = 2
    gdat.listindxenercntr = linspace(gdat.numbenerwndw, gdat.numbener - gdat.numbenerwndw, gdat.numbenercntr).astype(int)

    # get data counts
    if datatype == 'mock':
        gdat.datacnts = zeros_like(gdat.expo)
        ## temp -- FDM morphology should be energy dependent
        mockenertemppowr = gdat.meanener**(-2.6) * 1e5
        gdat.fdfmcnts = gdat.fdfmflux * gdat.expo * gdat.diffener[:, None, None] * gdat.apix * mockenertemppowr[:, None, None]
        gdat.mockcnts = gdat.fdfmcnts
        for i in gdat.indxener:
            for j in gdat.indxpixl:
                for m in gdat.indxevtt:
                    gdat.datacnts[i, j, m] = poisson(gdat.mockcnts[i, j, m])
    else:
        strgflux = '/flux%03d.fits' % gdat.numbside
        gdat.dataflux = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + strgflux)
        gdat.datacnts = gdat.dataflux * gdat.expo * gdat.diffener[:, None, None] * gdat.apix

    # plot data counts
    gdat.datacntsmean = sum(gdat.datacnts, 1)
    
    # model types
    gdat.listmodltype = ['altrmod0']
    #gdat.listmodltype = ['nullmod0', 'altrmod0']
    #listmodltype = ['nullmod1', 'altrmod1']
    
    gdat.bayefact = empty(gdat.numbenercntr)
    print 'hey'
    print 'gdat.bayefact'
    print gdat.bayefact

    # temp
    gdat.listindxenercntr = gdat.listindxenercntr[gdat.numbenercntr/2-1:gdat.numbenercntr/2+2]
    gdat.listenercntr = gdat.meanener[gdat.listindxenercntr]
    gdat.numbenercntr = gdat.listenercntr.size
    gdat.listlevi = empty(2)
    
    for gdat.thisindxenercntr, gdat.enercntrwndw in enumerate(gdat.listenercntr):

        # get window variables
        gdat.indxenercntrwndw = gdat.listindxenercntr[gdat.thisindxenercntr]
        gdat.minmindxenerwndw = gdat.indxenercntrwndw - gdat.numbenerwndwside
        gdat.maxmindxenerwndw = gdat.indxenercntrwndw + gdat.numbenerwndwside
        gdat.indxenerwndw = arange(gdat.minmindxenerwndw, gdat.maxmindxenerwndw + 1)

        indxmesh = meshgrid(gdat.indxenerwndw, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        gdat.datacntswndw = gdat.datacnts[indxmesh]
        gdat.expowndw = gdat.expo[indxmesh]
        gdat.fdfmfluxwndw = gdat.fdfmflux[indxmesh]
        
        indxmesh = meshgrid(gdat.indxenerwndw, gdat.indxevtt, indexing='ij')
        gdat.datacntsmeanwndw = gdat.datacntsmean[indxmesh]

        gdat.enercntrwndw = gdat.meanener[gdat.indxenercntrwndw]
        
        gdat.minmenerwndw = gdat.binsener[gdat.minmindxenerwndw]
        gdat.maxmenerwndw = gdat.binsener[gdat.maxmindxenerwndw+1]
        
        gdat.meanenerwndw = gdat.meanener[gdat.minmindxenerwndw:gdat.maxmindxenerwndw+1]
        gdat.binsenerwndw = gdat.binsener[gdat.minmindxenerwndw:gdat.maxmindxenerwndw+2]
        gdat.diffenerwndw = gdat.diffener[gdat.minmindxenerwndw:gdat.maxmindxenerwndw+1]

        gdat.edfnwndw, gdat.edfncom0wndw, gdat.edfncom1wndw = retr_fermedfn(gdat, gdat.meanenerwndw, gdat.enercntrwndw)
        
        gdat.strgenercntrwndw = gdat.strgenercntr[gdat.indxenercntrwndw]
   
        plot_edfnwndw(gdat)
        plot_datacntsmean(gdat)

        if methtype == 'almc':
            almc(gdat)
        if methtype == 'diff':
            diff(gdat)

    figr, axis = plt.subplots()
    print 'hey'
    print 'figure'
    print 'gdat.listenercntr'
    print gdat.listenercntr
    print 'gdat.bayefact'
    print gdat.bayefact
    axis.plot(gdat.listenercntr, gdat.bayefact)
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_ylabel('BF')
    plt.savefig(gdat.pathplot + 'bayefact.png')
    plt.close(figr) 
        

def diff(gdattemp):
    global gdat
    gdat = gdattemp

    for k, gdat.modltype in enumerate(gdat.listmodltype):

        for j in indxpixl:
            
            sampbund = tdpy.mcmc.init(numbproc, numbswep, retr_llik, datapara, thissamp=thissamp, numbburn=numbburn, \
                factthin=factthin, optiprop=optiprop, verbtype=verbtype, pathbase=pathbase, rtag=gdat.rtag, numbplotside=numbplotside)

    diffmaps = flux


def almc(gdattemp):
    global gdat
    gdat = gdattemp

    for k, gdat.modltype in enumerate(gdat.listmodltype):

        gdat.rtag = '%s_%s_%02d' % (gdat.modltype, gdat.strgenercntrwndw, gdat.maxmsphl)
                
        datapara = retr_datapara(gdat)
        namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varindxpara, dictpara = datapara
        gdat.numbpara = len(lablpara)

        gdat.numbswep = 100 * gdat.numbpara
        gdat.plotperd = gdat.numbswep / 10
        gdat.numbburn = gdat.numbswep / 10
        gdat.numbsamp = tdpy.mcmc.retr_numbsamp(gdat.numbswep, gdat.numbburn, gdat.factthin)

        pathmedisamp = os.environ["FERM_LINE_DATA_PATH"] + '/medisamp%s.fits' % gdat.rtag
        thissamp = empty((gdat.numbproc, gdat.numbpara))
        if os.path.isfile(pathmedisamp):
            print 'Loading initial sample from the previous run.'
            thissampvarb = pf.getdata(pathmedisamp)
            thissamp[:] = tdpy.mcmc.cdfn_samp(thissampvarb, datapara)[None, :]
        else:
            if gdat.modltype == 'nullmod0' or gdat.modltype == 'altrmod0' or gdat.modltype == 'nullmod1' or gdat.modltype == 'altrmod1':
                thissamp[:, 0:2] = rand(2 * gdat.numbproc)
                thissamp[:, 2:] = 0.5
                
        if gdat.numbpara > 10:
            numbplotside = 10
        else:
            numbplotside = gdat.numbpara
        sampbund = tdpy.mcmc.init(gdat.numbproc, gdat.numbswep, retr_llik, datapara, thissamp=thissamp, numbburn=gdat.numbburn, \
            factthin=gdat.factthin, optiprop=gdat.optiprop, verbtype=gdat.verbtype, pathbase=gdat.pathbase, rtag=gdat.rtag, numbplotside=numbplotside)

        listsampvarb = sampbund[0]
        listsamp = sampbund[1]
        listsampcalc = sampbund[2]
        listllik = sampbund[3]
        listaccp = sampbund[4]
        listjsampvari = sampbund[5]
        propeffi = sampbund[6]
        levi = sampbund[7]
        info = sampbund[8]
            
        gdat.listlevi[k] = levi
    
        medisampvarb = percentile(listsampvarb, 50., axis=0)
        pf.writeto(pathmedisamp, medisampvarb, clobber=True)
        
        gdat.modlcntswndw = retr_modlcnts(medisampvarb)
        gdat.resicntswndw = gdat.datacntswndw - gdat.modlcntswndw
        plot_cnts(gdat.datacntswndw, 'datacnts')
        plot_cnts(gdat.modlcntswndw, 'modlcnts')
        plot_cnts(gdat.resicntswndw, 'resicnts')

    print 'hey'
    print 'gdat.thisindxenercntr'
    print gdat.thisindxenercntr
    print 'gdat.bayefact'
    print gdat.bayefact
    print 
    gdat.bayefact[gdat.thisindxenercntr] = exp(gdat.listlevi[1] - gdat.listlevi[0])
    

if __name__ == '__main__':
    
    pass

    #writ_maps()
    #make_maps()
    #prep_maps()
    init(methtype='almc')
