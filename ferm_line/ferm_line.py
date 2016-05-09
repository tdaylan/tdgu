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


def retr_axes():

    reco = 8
    evtc = 128
    
    numbener = 20
    minmener = 55. # [GeV]
    maxmener = 65. # [GeV]
    binsener = logspace(log10(minmener), log10(maxmener), numbener + 1)
    meanener = sqrt(binsener[1:] * binsener[0:-1])
    diffener = binsener[1:] - binsener[0:-1]
    indxener = arange(numbener)
    
    global indxevtt
    evtt = array([64, 128, 256, 512])
    numbevtt = evtt.size
    indxevtt = arange(numbevtt)
    
    numbside = 8
    numbpixl = numbside**2 * 12
    apix = 4. * pi / numbpixl
    
    return reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix


def make_maps():
    
    global bind
    bind = False
    
    cmnd = 'mkdir -p $FERMI_DATA/exposure/ferm_line'
    os.system(cmnd)
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix = retr_axes()

    global liststrgener, liststrgchan
    liststrgener = ['ENERGY%d' % (k + 1) for k in indxener]
    liststrgchan = ['CHANNEL%d' % (k + 1) for k in indxener]
    
    numbproc = 20
    minmweek = 11
    maxmweek = 411
    weekpart = (maxmweek - minmweek) / numbproc
    
    global listweek
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

            cmnd = 'gtselect infile=' + evnt + ' outfile=' + sele +                 ' ra=INDEF dec=INDEF rad=INDEF tmin=INDEF tmax=INDEF' +                 ' emin=%d emax=%d zmax=90 evclass=%d evtype=%d' % (minmener * 1e3, maxmener * 1e3, evtc, evtt[m])
            os.system(cmnd)

            cmnd = 'gtmktime evfile=' + sele + ' scfile=' + spac + ' filter="DATA_QUAL==1 && LAT_CONFIG==1"' +                 ' outfile=' + filt + ' roicut=no'

            os.system(cmnd)
            
            if bind or (t == 0 and m == 0):
                cntstemp = os.environ["FERMI_DATA"] + 'exposure/ferm_line/cnts.fits'
                cmnd = 'gtbin evfile=' + filt + ' scfile=' + spac + ' outfile=' + cntstemp +                     ' ebinalg=LOG emin=%d emax=%d enumbins=%d algorithm=HEALPIX'                     % (minmener * 1e3, maxmener * 1e3, numbener) +                     ' hpx_ordering_scheme=RING coordsys=GAL hpx_order=8 hpx_ebin=yes'
                os.system(cmnd)
            else:
                cntstemp = cnts

            cmnd = 'gtltcube evfile=' + filt + ' scfile=' + spac + ' outfile=' + live +                 ' dcostheta=0.025 binsz=1'
            os.system(cmnd)

            cmnd = 'gtexpcube2 infile=' + live + ' outfile=' + expo + ' cmap=' + cntstemp +                 ' irfs=CALDB evtype=%03d bincalc=CENTER' % evtt[m]
            os.system(cmnd)


def prep_maps():
    

    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix = retr_axes()

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


def plot_datacntsmean(thisdatacntsmean, rtag):

    figr, axcl = plt.subplots(numbevtt, 1, figsize=(10, 20))
    for m, axis in enumerate(axcl):
        xdat = meanener
        ydat = thisdatacntsmean[:, m]
        yerr = sqrt(thisdatacntsmean[:, m])
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', elinewidth=2, capsize=20)
        axis.set_xlabel('E [GeV]')
        axis.set_xlim([minmener, maxmener])
        axis.set_title('EDISP%d' % m)
    figr.text(0.05, 0.5, '$N$', ha='center', va='center', rotation=90)
    plt.subplots_adjust(hspace=0.1)
    plt.savefig(plotpath + 'datameancnts' + rtag + '.png')
    plt.close(figr)
          

def retr_datapara():


    if modltype == 'nullmod0':
        numbpara = 2
    if modltype == 'altrmod0':
        numbpara = 2 + numbalmp
        
    if modltype == 'nullmod1':
        numbpara = 2 + numbalmp
    if modltype == 'altrmod1':
        numbpara = 2 + 2 * numbalmp

    if modltype == 'nullmod2':
        numbpara = 1 + numbalmp
    if modltype == 'altrmod2':
        numbpara = 1 + 2 * numbalmp
        
    dictpara = dict()
    minmpara = zeros(numbpara)
    maxmpara = zeros(numbpara)
    namepara = empty(numbpara, dtype=object)
    scalpara = empty(numbpara, dtype=object)
    lablpara = empty(numbpara, dtype=object)
    unitpara = empty(numbpara, dtype=object)
    varipara = zeros(numbpara)
    
    dictpara['slopback'] = 0
    namepara[0] = 'slopback'
    minmpara[0] = 1.5
    maxmpara[0] = 4.
    scalpara[0] = 'self'
    lablpara[0] = r'$\alpha_b$'
    unitpara[0] = ''
    varipara[0] = 1e-1

    cntrpara = 1
    if modltype == 'altrmod0' or modltype == 'nullmod0' or modltype == 'altrmod1' or modltype == 'nullmod1':
            
        dictpara['normback'] = 1
        namepara[1] = 'normback'
        minmpara[1] = 1e-1
        maxmpara[1] = 1e1
        scalpara[1] = 'logt'
        lablpara[1] = r'$A_b$'
        unitpara[1] = ''
        varipara[1] = 1e-2

        cntrpara += 1
    
    if modltype == 'nullmod1' or modltype == 'altrmod1' or modltype == 'nullmod2' or modltype == 'altrmod2':
        defn_almcpara('powrspec', cntrpara, dictpara, namepara,                       minmpara, maxmpara, scalpara, lablpara, unitpara, varipara)
        cntrpara += numbalmp
        
    if modltype == 'altrmod0' or modltype == 'altrmod1' or modltype == 'altrmod2':
        defn_almcpara('linespec', cntrpara, dictpara, namepara,                       minmpara, maxmpara, scalpara, lablpara, unitpara, varipara)
        
    strgpara = lablpara + ' ' + unitpara
    datapara = namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara
    
    return datapara


def defn_almcpara(typepara, cntrpara, dictpara, namepara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara):
    
    global listsphm
    listsphl, listsphm = hp.Alm.getlm(maxmsphl)
    listsphl = tile(listsphl, 2)
    listsphm = tile(listsphm, 2)

    temp = []
    for k in range(numbalmc):
        if listsphm[numbalmc+k] == 0:
            temp.append(numbalmc + k)
    listsphl = delete(listsphl, temp)
    listsphm = delete(listsphm, temp)

    maxmindxpara = cntrpara + numbalmp
        
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
        if k < cntrpara + numbalmc:
            lablpara[k] = 'Re{'
        else:
            lablpara[k] = 'Im{' 
        lablpara[k] += '%s_{%d%d}$' % (strg, listsphl[k - cntrpara],   listsphm[k - cntrpara])
        if k < numbalmc:
            lablpara[k] += '}'
        else:
            lablpara[k] += '}'
        unitpara[k] = ''
        varipara[k] = 1e-2


def retr_aexp(scaldevi, stdv, skew, bias, slop):

    
    aexp = slop / stdv / sp.special.gamma(1. / slop) * skew / (1. + skew**2)
    indxscalfrwd = where(scaldevi - bias > 0.)
    indxscalback = where(scaldevi - bias <= 0.)
    aexp[indxscalfrwd] *= exp(-(skew[indxscalfrwd] / stdv[indxscalfrwd])**slop[indxscalfrwd] * (scaldevi[indxscalfrwd] - bias[indxscalfrwd])**slop[indxscalfrwd])
    aexp[indxscalback] *= exp(-(1. / skew[indxscalback] / stdv[indxscalback])**slop[indxscalback] * (bias[indxscalback] - scaldevi[indxscalback])**slop[indxscalback])

    return aexp


def retr_fermedfn(thisener, cntrener=None):
  
    if cntrener == None:
        cntrener = sqrt(thisener[0] * thisener[-1])
    enerdevi = thisener - cntrener

    numbevtt = 4
    indxevtt = arange(numbevtt)
    
    thisnumbener = thisener.size

    ctht = linspace(0.7, 1., 20)

    path = os.environ["PCAT_DATA_PATH"] + '/irf/edisp_P8R2_SOURCE_V6_EDISP.fits'

    #strgformpara = ['F', 'S1', 'K1', 'BIAS', 'PINDEX1', 'S2', 'K2', 'BIAS2', 'PINDEX2']
    strgformpara = ['F', 'S1', 'S2', 'K1', 'K2', 'BIAS', 'BIAS2', 'PINDEX1', 'PINDEX2']

    listhdun = pf.open(path)

    numbctht = 8
    numbformpara = 9
    numbscalpara = 6
    frac = zeros((thisnumbener, numbctht, numbevtt))
    stdv = zeros((2, thisnumbener, numbctht, numbevtt))
    skew = zeros((2, thisnumbener, numbctht, numbevtt))
    bias = zeros((2, thisnumbener, numbctht, numbevtt))
    slop = zeros((2, thisnumbener, numbctht, numbevtt))
    scalpara = zeros((numbevtt, numbscalpara))
    for m in indxevtt:
        
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
            formparatemp = interp1d(meanenerirfn, hdun.data[strgformpara[k]].squeeze(), axis=1)(thisener).T
            if k == 0:
                frac[:, :, m] = formparatemp
            if k == 1:
                stdv[0, :, :, m] = formparatemp
            if k == 2:
                stdv[1, :, :, m] = formparatemp
            if k == 3:
                skew[0, :, :, m] = formparatemp
            if k == 4:
                skew[1, :, :, m] = formparatemp
            if k == 5:
                bias[0, :, :, m] = formparatemp
            if k == 6:
                bias[1, :, :, m] = formparatemp
            if k == 7:
                slop[0, :, :, m] = formparatemp
            if k == 8:
                slop[1, :, :, m] = formparatemp

        # scale parameters
        scalpara[m, :] = listhdun['EDISP_SCALING_PARAMS_EDISP%d' % m].data['EDISPSCALE']

    scalfact = scalpara[None, None, :, 0] * log(thisener[:, None, None])**2 + \
               scalpara[None, None, :, 1] * cthtirfn[None, :, None]**2 + \
               scalpara[None, None, :, 2] * log(thisener[:, None, None]) + \
               scalpara[None, None, :, 3] * cthtirfn[None, :, None] + \
               scalpara[None, None, :, 4] * log(thisener[:, None, None]) * cthtirfn[None, :, None] + \
               scalpara[None, None, :, 5]
            
    scaldevi = enerdevi[:, None, None] / cntrener / scalfact
    
    edfn = frac * retr_aexp(scaldevi, stdv[0, :, :, :], skew[0, :, :, :], bias[0, :, :, :], slop[0, :, :, :]) + \
                    (1. - frac) * retr_aexp(scaldevi, stdv[1, :, :, :], skew[1, :, :, :], bias[1, :, :, :], slop[1, :, :, :])
  
    edfn = mean(edfn, axis=1)

    return edfn


def plot_fermedfn():

    cntrener = 60.
    meanener = linspace(55., 65., 10)
    fermedfn = retr_fermedfn(meanener, cntrener=cntrener)

    figr, axis = plt.subplots()
    axis.plot(meanener, fermedfn)
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    plt.savefig(plotpath + 'fermedfn_%.3g.png' % cntrener)
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
    
    linecnts = 1. / sqrt(2. * pi * stdvdisp[None, :]**2) *         exp(-0.5 * (thisener[:, None] - enercntr)**2 / stdvdisp[None, :]**2)
    
    return linecnts


def retr_modlcnts(sampvarb):
    
    global almcimag
    
    modlcnts = zeros((numbener, numbpixl, numbevtt))
    
    slopback = sampvarb[0]
    enertemppowr = (meanener / enercntr)**(-slopback)
    
    cntr = 1
    
    if modltype != 'nullmod2' and modltype != 'altrmod2':
        
        normback = sampvarb[1]
        modlcnts = normback * fdfmcnts * enertemppowr[:, None, None]
        cntr += 1
        
    if modltype == 'nullmod0':
        numbiter = 0
    elif modltype == 'altrmod0' or modltype == 'nullmod1' or modltype == 'nullmod2':
        numbiter = 1
    else:
        numbiter = 2
        
    for k in range(numbiter):

        almcreal = sampvarb[cntr:cntr+numbalmc]
        almcimag[maxmsphl+1:] = sampvarb[cntr+numbalmc:cntr+numbalmp]
        almc = almcreal + almcimag * 1j
        modlcntstemp[:] = hp.alm2map(almc, numbside, verbose=False)[None, :, None]
        cntr += numbalmp
        if k == 0:
            modlcnts += modlcntstemp * enertempline[:, None, :]
        else:
            modlcnts += modlcntstemp * enertemppowr[:, None, None]
    
    if lliktype == 'bind':
        llik = sum(datacnts * log(modlcnts) - modlcnts)
    else:
        llik = sum(exp(-modlcntsintp))
    
    return modlcnts


def retr_llik(sampvarb, init=False):
    
    modlcnts = retr_modlcnts(sampvarb)
   
    if lliktype == 'bind':
        llik = sum(datacnts * log(modlcnts) - modlcnts)
    else:
        llik = sum(exp(-modlcntsintp))
    
    return llik, sampcalc
    

def plot_enertempline():
    
    figr, axis = plt.subplots()
    axis.plot(meanener, enertempline)
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    plt.savefig(plotpath + 'enertempline_' + strgenercntr + '.png')
    plt.close(figr) 


def plot_cnts(thiscnts, strg):

    cart = tdpy.util.retr_cart(sum(thiscnts[indxenercntr, :, :], 1))

    figr, axis = plt.subplots()
    
    imag = plt.imshow(cart, origin='lower', cmap='Reds', extent=extt)
    plt.colorbar(imag, fraction=0.05)

    plt.savefig(plotpath + '%scnts%s.png' % (strg, rtag))
    plt.close(figr)


def retr_numbalmc(maxmsphl):
    
    numbalmc = maxmsphl * (maxmsphl + 1) / 2 + maxmsphl + 1
    
    return numbalmc

def retr_numbalmp(maxmsphl, numbalmc):
    
    numbalmp = numbalmc * 2 - maxmsphl - 1
    
    return numbalmp


def plot_sphl():

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
    plt.savefig(plotpath + 'sphl.png')

    
def init():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, numbside, numbpixl, apix = retr_axes()

    global plotpath, plotpath
    if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
        plotfold = '/n/pan/www/tansu/png/ferm_line/'
    else:
        plotfold = os.environ["FERM_LINE_DATA_PATH"] + '/png/'
    plotpath = plotfold
    cmnd = 'mkdir -p ' + plotpath
    os.system(cmnd)
    
    global extt
    extt = [-180., 180., -90., 90.]
    
    global numbalmc, maxmsphl, numbalmp
    maxmsphl = 4
    numbalmc = retr_numbalmc(maxmsphl)
    numbalmp = retr_numbalmp(maxmsphl, numbalmc)
    
    plot_sphl()
    
    global sampcalc
    sampcalc = []
   
    global lliktype
    lliktype = 'bind'

    global modlcntstemp
    modlcntstemp = empty((numbener, numbpixl, numbevtt))
    
    global datacnts, datacntsmean

    global enerpivt, indxenerpivt
    indxenerpivt = numbener / 2
    enerpivt = meanener[indxenerpivt]
    
    strgflux = '/flux%03d.fits' % numbside
    strgexpo = '/expo%03d.fits' % numbside   
    dataflux = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + strgflux)
    expo = pf.getdata(os.environ["FERM_LINE_DATA_PATH"] + strgexpo)
    
    global fdfmflux, fdfmcnts
    path = os.environ["FERM_LINE_DATA_PATH"] + '/fdfm%03d.fits' % numbside
    fdfmflux = pf.getdata(path)
    fdfmcnts = fdfmflux[None, :, None] * expo * diffener[:, None, None] * apix
    datacnts = dataflux * expo * diffener[:, None, None] * apix
    datacntsmean = sum(datacnts, 1)
    plot_datacntsmean(datacntsmean, '_init')
    
    global stdvdisp
    stdvdisp = array([7., 5., 3., 1.])
   
    # common MCMC settings 
    verbtype = 1
    factthin = 1
    numbproc = 2
    optiprop = True
    
    global almcimag
    almcimag = zeros(numbalmc)
   

    
    plot_fermedfn()

    global enercntr, enertempline, modltype, strgenercntr, rtag, indxenercntr
    listindxenercntr = array([numbener / 2]) # array([58., 60., 62.])
    listenercntr = meanener[listindxenercntr]
    numbenercntr = listenercntr.size
    listlevi = zeros((2, numbenercntr))
    for n, enercntr in enumerate(listenercntr):

        indxenercntr = listindxenercntr[n]
        
        enertempline = retr_linecnts(meanener, enercntr)

        strgenercntr = '%04d' % enercntr
        
        plot_enertempline()
  
        listmodltype = ['nullmod2', 'altrmod2']
        for k, modltype in enumerate(listmodltype):

            rtag = '_%s_%s_%02d' % (modltype, strgenercntr, maxmsphl)
                    
            datapara = retr_datapara()
            namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varindxpara, dictpara = datapara
            numbpara = len(lablpara)

            numbswep = 10000 * numbpara
            plotperd = numbswep / 10
            numbburn = numbswep / 10
            numbsamp = tdpy.mcmc.retr_numbsamp(numbswep, numbburn, factthin)

            pathmedisamp = os.environ["FERM_LINE_DATA_PATH"] + '/medisamp%s.fits' % rtag
            thissamp = empty((numbproc, numbpara))
            if os.path.isfile(pathmedisamp):
                print 'Loading initial sample from the previous run.'
                thissampvarb = pf.getdata(pathmedisamp)
                thissamp[:] = tdpy.mcmc.cdfn_samp(thissampvarb, datapara)[None, :]
            else:
                if modltype == 'nullmod0' or modltype == 'altrmod0' or modltype == 'nullmod1' or modltype == 'altrmod1':
                    thissamp[:, 0:2] = rand(2 * numbproc)
                    thissamp[:, 2:] = 0.5
                if modltype == 'nullmod2' or modltype == 'altrmod2':
                    thissamp[:, 0] = rand(numbproc)
                    thissamp[:, 1] = 1.
                    thissamp[:, 2:] = 0.5
                    
            if numbpara > 10:
                numbplotside = 10
            else:
                numbplotside = numbpara
            
            
            sampbund = tdpy.mcmc.init(numbproc, numbswep, retr_llik, datapara, thissamp=thissamp, numbburn=numbburn, \
                factthin=factthin, optiprop=optiprop, verbtype=verbtype, plotpath=plotpath, rtag=rtag, numbplotside=numbplotside)

            listsampvarb = sampbund[0]
            listsamp = sampbund[1]
            listsampcalc = sampbund[2]
            listllik = sampbund[3]
            listaccp = sampbund[4]
            listjsampvari = sampbund[5]
            propeffi = sampbund[6]
            levi = sampbund[7]
            info = sampbund[8]
                
            listlevi[k, n] = levi
        
            medisampvarb = percentile(listsampvarb, 50., axis=0)
            pf.writeto(pathmedisamp, medisampvarb, clobber=True)
            
            medimodlcnts = retr_modlcnts(medisampvarb)
            plot_cnts(medimodlcnts, 'medimodl')
            plot_cnts(datacnts - medimodlcnts, 'mediresi')

    print listlevi.T
    
    bayefact = exp(listlevi[1, :] - listlevi[0, :])
    
    figr, axis = plt.subplots()
    axis.plot(listenercntr, bayefact)
    axis.set_xlabel(r'$E_\gamma$ [GeV]')
    axis.set_ylabel('BF')
    plt.savefig(plotpath + 'bayefact.png')
    plt.close(figr) 
        

if __name__ == '__main__':
    
    pass

    #writ_maps()
    #make_maps()
    #prep_maps()
    init()
