import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
mpl.rc('image', interpolation='none', origin='lower')

# numpy
import numpy as np
from numpy import *
from numpy.random import *
from numpy.random import choice

import matplotlib.pyplot as plt

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

import sympy, functools

def retr_axes():

    reco = 8
    evtc = 128
    
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    meanener = sqrt(binsener[1:] * binsener[0:-1])
    diffener = binsener[1:] - binsener[0:-1]
    numbener = meanener.size
    indxener = arange(numbener)
    minmener = amin(binsener)
    maxmener = amax(binsener)
    
    global indxevtt
    evtt = array([4, 16, 32, 64])
    numbevtt = evtt.size
    indxevtt = arange(numbevtt)
    
    nside = 256
    numbpixl = nside**2 * 12
    apix = 4. * pi / numbpixl
    
    return reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix


def make_maps_pcat():
    
    #reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    gdat = tdpy.util.gdatstrt()
    
    numbtime = 4
    gdat.rtag = ['cmp%d' % (numbtime - t - 1) for t in range(numbtime)] + ['full']
    gdat.reco = [7 for t in range(numbtime)] + [8]
    gdat.evtc = [2 for t in range(numbtime)] + [128]
    gdat.strgtime = ['tmin=239155201 tmax=364953603' for t in range(numbtime)] + ['tmin=INDEF tmax=INDEF']
    gdat.weekinit = [9 for t in range(numbtime)] + [11]
    gdat.weekfinl = [218 for t in range(numbtime)] + [411]
    gdat.listtimefrac = [0.25, 0.50, 0.75, 1.] + [1.]
    gdat.photpath = ['p7v6c' for t in range(numbtime)] + ['photon']
    gdat.strgregi = [' ra=INDEF dec=INDEF rad=INDEF ' for n in range(numbtime + 1)]
    gdat.strgregi = [256 for n in range(numbtime + 1)]

    gdat.evtt = [4, 8, 16, 32]
    gdat.test = True

    make_maps_main(gdat)


def make_maps_main(gdat):
    
    numbproc = len(gdat.rtag)
    indxproc = arange(numbproc)
    
    # process pool
    pool = mp.Pool(numbproc)

    # spawn the processes
    make_maps_part = functools.partial(make_maps_work, gdat)
    pool.map(make_maps_part, indxproc)
    pool.close()
    pool.join()


def make_maps_work(gdat, indxprocwork):

    cmnd = 'mkdir -p $FERMI_DATA/exposure/%s' % gdat.rtag[indxprocwork]
    os.system(cmnd)

    # make file lists
    infl = '$PCAT_DATA_PATH/phot_%s.txt' % gdat.rtag[indxprocwork]
    spac = '$PCAT_DATA_PATH/spac_%s.txt' % gdat.rtag[indxprocwork]
        
    numbweek = (gdat.weekfinl[indxprocwork] - gdat.weekinit[indxprocwork]) * gdat.listtimefrac[indxprocwork]
    listweek = floor(linspace(gdat.weekinit[indxprocwork], gdat.weekfinl[indxprocwork] - 1, numbweek)).astype(int)
    cmnd = 'rm ' + infl
    os.system(cmnd)
    cmnd = 'rm ' + spac
    os.system(cmnd)
    for week in listweek:
        cmnd = 'ls -d -1 $FERMI_DATA/weekly/spacecraft/*_w%03d_* >> ' % week + spac
        os.system(cmnd)
        cmnd = 'ls -d -1 $FERMI_DATA/weekly/%s/*_w%03d_* >> ' % (gdat.photpath[indxprocwork], week) + infl
        os.system(cmnd)
    for m in indxevtt:

        if gdat.reco[indxprocwork] == 7:
            if m == 3:
                thisevtt = 1
                thisevttdepr = 0
            elif m == 2:
                thisevtt = 2
                thisevttdepr = 1
            else:
                continue
            strgpsfn = 'convtype=%d' % thisevttdepr

        if gdat.reco[indxprocwork] == 8:
            thisevtt = gdat.evtt[m]
            strgpsfn = 'evtype=%d' % thisevtt
         
        sele = '$PCAT_DATA_PATH/sele_evtt%03d_%s.fits' % (thisevtt, gdat.rtag[indxprocwork])
        filt = '$PCAT_DATA_PATH/filt_evtt%03d_%s.fits' % (thisevtt, gdat.rtag[indxprocwork])
        live = '$PCAT_DATA_PATH/live_evtt%03d_%s.fits' % (thisevtt, gdat.rtag[indxprocwork])
        cnts = '$PCAT_DATA_PATH/cnts_evtt%03d_%s.fits' % (thisevtt, gdat.rtag[indxprocwork])
        expo = '$PCAT_DATA_PATH/expo_evtt%03d_%s.fits' % (thisevtt, gdat.rtag[indxprocwork])

        cmnd = 'gtselect infile=' + infl + ' outfile=' + sele + gdat.strgregi[indxprocwork] + \
            gdat.strgtime[indxprocwork] + ' emin=100 emax=100000 zmax=90 evclass=%d %s' % (gdat.evtc[indxprocwork], strgpsfn)
        if gdat.test:
            print cmnd
            print ''
        else:
            os.system(cmnd)

        cmnd = 'gtmktime evfile=' + sele + ' scfile=' + spac + ' filter="DATA_QUAL==1 && LAT_CONFIG==1"' + ' outfile=' + filt + ' roicut=no'
        if gdat.test:
            print cmnd
            print ''
        else:
            os.system(cmnd)

        cmnd = 'gtbin evfile=' + filt + ' scfile=NONE outfile=' + cnts + \
            ' ebinalg=FILE ebinfile=/n/fink1/fermi/exposure/gcps_time/gtbndefn.fits ' + \
            'algorithm=HEALPIX hpx_ordering_scheme=RING coordsys=GAL hpx_order=%d hpx_ebin=yes' % log2(gdat.numbside[indxprocwork])
        if gdat.test:
            print cmnd
            print ''
        else:
            os.system(cmnd)

        cmnd = 'gtltcube evfile=' + filt + ' scfile=' + spac + ' outfile=' + live + ' dcostheta=0.025 binsz=1'
        if gdat.test:
            print cmnd
            print ''
        else:
            os.system(cmnd)

        cmnd = 'gtexpcube2 infile=' + live + ' cmap=' + cnts + ' outfile=' + expo + ' irfs=CALDB evtype=%03d bincalc=CENTER' % thisevtt
        if gdat.test:
            print cmnd
            print ''
        else:
            os.system(cmnd)


def prep_maps():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    liststrgener = ['ENERGY1', 'ENERGY2', 'ENERGY3', 'ENERGY4', 'ENERGY5']
    liststrgchan = ['CHANNEL1', 'CHANNEL2', 'CHANNEL3', 'CHANNEL4', 'CHANNEL5']
    listdatatype = ['cmp0', 'cmp1', 'cmp2', 'cmp3', 'full']
    listregitype = ['igal', 'ngal']
    
    cnts = zeros((numbener, numbpixl, numbevtt))
    expo = zeros((numbener, numbpixl, numbevtt))
    flux = zeros((numbener, numbpixl, numbevtt))
    
    for k, datatype in enumerate(listdatatype):
    
        for m in indxevtt:

            if datatype != 'full':
                if m < 2:
                    continue
                elif m == 2:
                    thisevtt = 2
                elif m == 3:
                    thisevtt = 1
            else:
                thisevtt = evtt[m]

            path = os.environ["PCAT_DATA_PATH"] + '/expo_evtt%03d_%s.fits' % (thisevtt, datatype)
            expoarry = pf.getdata(path, 1)
            for i in indxener:
                expo[i, :, m] = expoarry[liststrgener[i]]

            path = os.environ["PCAT_DATA_PATH"] + '/cnts_evtt%03d_%s.fits' % (thisevtt, datatype)
            cntsarry = pf.getdata(path)
            for i in indxener:
                cnts[i, :, m] = cntsarry[liststrgchan[i]]

        indxexpo = where(expo > 0.) 
        flux[indxexpo] = cnts[indxexpo] / expo[indxexpo] / apix
        flux /= diffener[:, None, None]

        for regitype in listregitype:
            if regitype == 'ngal':
                for i in indxener:
                    for m in indxevtt:
                        
                        if datatype != 'full':
                            if m < 2:
                                continue
                            elif m == 2:
                                thisevtt = 2
                            elif m == 3:
                                thisevtt = 1
                        else:
                            thisevtt = evtt[m]

                        almc = hp.map2alm(flux[i, :, m])
                        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                        flux[i, :, m] = hp.alm2map(almc, nside)

                        almc = hp.map2alm(expo[i, :, m])
                        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                        expo[i, :, m] = hp.alm2map(almc, nside)

            path = os.environ["PCAT_DATA_PATH"] + '/fermexpo_%s_%s.fits' % (datatype, regitype)
            pf.writeto(path, expo, clobber=True)

            path = os.environ["PCAT_DATA_PATH"] + '/fermflux_%s_%s.fits' % (datatype, regitype)
            pf.writeto(path, flux, clobber=True)


def writ_isot():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    # isotropic background
    path = os.environ["PCAT_DATA_PATH"] + '/iso_P8R2_ULTRACLEAN_V6_v06.txt'
    isotdata = loadtxt(path)
    enerisot = isotdata[:, 0] * 1e-3 # [GeV]
    isotflux = isotdata[:, 1] * 1e3 # [1/cm^2/s/sr/GeV]
    isotfluxheal = empty((numbener, numbpixl, numbevtt))
    nsampbins = 10
    enersamp = logspace(log10(amin(binsener)), log10(amax(binsener)), nsampbins * numbener)
    isotflux = interpolate.interp1d(enerisot, isotflux)(enersamp)
    for i in range(numbener):
        isotfluxheal[i, :, :] = trapz(isotflux[i*nsampbins:(i+1)*nsampbins], enersamp[i*nsampbins:(i+1)*nsampbins]) / diffener[i]
        
    path = os.environ["PCAT_DATA_PATH"] + '/fermisotflux.fits'
    pf.writeto(path, isotfluxheal, clobber=True)

    nfwpfluxtemp = tdpy.util.retr_nfwp(1., nside, norm=5.)
    nfwpspec = ones(numbener)
    nfwpfluxheal = zeros((numbener, numbpixl, numbevtt))
    for i in indxener:
        for m in indxevtt:
            nfwpfluxheal[i, :, m] = nfwpspec[i] * nfwpfluxtemp
    path4 = os.environ["PCAT_DATA_PATH"] + '/nfwpflux.fits'
    pf.writeto(path4, nfwpfluxheal, clobber=True)


def writ_fdfm_doug():
    
    numbevtt = 4
    nside = 256
    numbener = 5
    numbpixl = nside**2 * 12
    fermfdfmfluxigal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxngal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxigal[0, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff0.fits')[:, None]
    fermfdfmfluxigal[1, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff1.fits')[:, None]
    fermfdfmfluxigal[2, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff2.fits')[:, None]
    fermfdfmfluxigal[3, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff3.fits')[:, None]
    fermfdfmfluxigal[4, :, :] = pf.getdata(os.environ["PCAT_DATA_PATH"] + '/diff4.fits')[:, None]
    
    path =os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    pf.writeto(path, fermfdfmfluxigal, clobber=True)
    
    for i in range(numbener):
        
        almc = hp.map2alm(fermfdfmfluxigal[i, :, 0])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        fermfdfmfluxngal[i, :, :] = hp.alm2map(almc, nside)[:, None]
        
    path =os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_ngal.fits'
    pf.writeto(path, fermfdfmfluxngal, clobber=True)
        

def writ_fdfm():
    
    nside = 256
    numbpixl = 12 * nside**2
    numbener = 3
    numbevtt = 4

    binsener = array([0.3, 1., 3., 10.])
    
    fermfdfmfluxigaltemp = tdpy.util.retr_fdfm(binsener, nside)

    fermfdfmfluxngaltemp = zeros((numbener, numbpixl))
    for i in range(numbener):
        almc = hp.map2alm(fermfdfmfluxigaltemp[i, :])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        fermfdfmfluxngaltemp[i, :] = hp.alm2map(almc, nside)

    fermfdfmfluxngal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxigal = zeros((numbener, numbpixl, numbevtt))
    for m in range(numbevtt):
        fermfdfmfluxigal[:, :, m] = fermfdfmfluxigaltemp
        fermfdfmfluxngal[:, :, m] = fermfdfmfluxngaltemp

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    pf.writeto(path, fermfdfmfluxigal, clobber=True)

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_ngal.fits'
    pf.writeto(path, fermfdfmfluxngal, clobber=True)


def prep_dust():

    minmlgal = -20.
    maxmlgal = 20.
    minmbgal = -20.
    maxmbgal = 20.

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    fdfmflux = pf.getdata(path)

    path = os.environ["PCAT_DATA_PATH"] + '/lambda_sfd_ebv.fits'
    dustigal = pf.getdata(path)['TEMPERATURE']
    numbside = int(sqrt(dustigal.size / 12))
    dustigal = hp.reorder(dustigal, n2r=True)

    dustigal = dustigal[None, :, None] * mean(fdfmflux, 1)[:, None, :] / mean(dustigal)

    path = os.environ["PCAT_DATA_PATH"] + '/fermdustflux_igal.fits'
    pf.writeto(path, dustigal, clobber=True)

    dustngal = empty_like(dustigal)
    numbener = 3
    for i in range(numbener):
        almc = hp.map2alm(dustigal[i, :, 0])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        dustngal[i, :, :] = hp.alm2map(almc, numbside)[:, None]

    path = os.environ["PCAT_DATA_PATH"] + '/fermdustflux_igal.fits'
    pf.writeto(path, dustigal, clobber=True)
    
    path = os.environ["PCAT_DATA_PATH"] + '/fermdustflux_ngal.fits'
    pf.writeto(path, dustngal, clobber=True)
    
    
def retr_jcbn():
    
    lgl0, lgla, bgl0, bgla, flx0, flxa, snd0, snda = sympy.symbols('lgl0 lgla bgl0 bgla flx0 flxa snd0 snda')
    matr = sympy.Matrix([[1, 1 - 1 / (flx0 + flxa), 0,               0, algl * aflx / flx0**2, -algl / flx0], \
                         [1,    -1 / (flx0 + flxa), 0,               0, algl * aflx / flx0**2, -algl / flx0], \
                         [0,                     0, 1, 1 - aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0], \
                         [0,                     0, 1,    -aflx / flx0, abgl * aflx / flx0**2, -abgl / flx0], \
                         [0,                     0, 0,               0,                     1,            0], \
                         [0,                     0, 0,               0,                    -1,            1]])

    jcbn = matr.det()
    return jcbn


def plot_maps():
    
    global numbpixl
    numbside = 256
    lgalheal, bgalheal, numbpixl, apix = tdpy.util.retr_healgrid(numbside)
    
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    meanener = sqrt(binsener[1:] * binsener[:-1])
    global minmlgal, maxmlgal, minmbgal, maxmbgal, extt
    minmlgal = -20.
    maxmlgal = 20.
    minmbgal = -20.
    maxmbgal = 20.

    path = os.environ["PCAT_DATA_PATH"] + '/fermflux_cmp0_ngal.fits'
    cmp0 = pf.getdata(path)
    path = os.environ["PCAT_DATA_PATH"] + '/fermflux_comp_ngal.fits'
    comp = pf.getdata(path)
    diff = cmp0 - comp
    for i in range(meanener.size):
        for m in range(4):
            
            path = '/n/pan/www/tansu/imag/pcat/cmp0_ene%d_evt%d.pdf' % (i, m)
            tdpy.util.plot_heal(cmp0[i, :, m], path=path, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
            
            path = '/n/pan/www/tansu/imag/pcat/comp_ene%d_evt%d.pdf' % (i, m)
            tdpy.util.plot_heal(comp[i, :, m], path=path, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)

            path = '/n/pan/www/tansu/imag/pcat/diff_ene%d_evt%d.pdf' % (i, m)
            tdpy.util.plot_heal(diff[i, :, m], path=path, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)




#writ_isot()
#writ_fdfm()
#writ_fdfm_doug()
#plot_maps()
#make_maps_pcat()
#prep_maps()
prep_dust()
#plot_maps()
