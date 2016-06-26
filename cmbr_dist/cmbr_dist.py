import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from numpy import *
from numpy.random import *

import scipy as sp
from scipy.optimize import optimize
from scipy.special import erf, zetac, gamma
from scipy.interpolate import *
from scipy.integrate import odeint  

import os
import pyfits as pf
import healpy as hp

import tdpy.util
import tdpy.mcmc

import emcee

import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)

def plot_intr(eventype='norm'):
    
    with plt.xkcd():

        def sigd(varb, varbcntr):
            sigd = 1. / (1 + exp(-0.1 * (varb - varbcntr)))
            return sigd / amax(sigd)

        fig, ax = plt.subplots(figsize=(14, 6))

        redsxkcd = arange(1000)
        rele = sigd(redsxkcd, 300.) - sigd(redsxkcd, 500.) + sigd(redsxkcd, 950.)

        ax.plot(redsxkcd, rele)

        ax.text(200, -0.1, "z=One Million")
        ax.text(450, -0.1, "z=One Thousand")
        ax.text(650, -0.1, "z=One Hundred")
        ax.text(900, -0.1, "04/01/2016")
        ax.set_xticks([])
        if eventype == 'scra':
            ax.annotate('Nucleosynthesis', xy=[950, 0.3], xytext=[20, 0.7], rotation=45)
            ax.annotate('Reheating', xy=[950, 0.3], xytext=[120, 0.7], rotation=45)
            ax.annotate('Inflation', xy=[950, 0.3], xytext=[200, 0.7], rotation=45)
            ax.annotate('CMB Spectral Distortions', xy=[950, 0.3], xytext=[300, 0.45], rotation=45)
            ax.annotate('Reionization', xy=[550, 0.3], xytext=[550, 0.7], rotation=45)
            ax.annotate('Structure formation', xy=[750, 0.3], xytext=[700, 0.7], rotation=45)
            ax.annotate('Recombination', xy=[750, 0.3], xytext=[650, 0.7], rotation=45)
            ax.set_title("History of the Universe")
        else:
            ax.annotate('Inflaton', xy=[950, 0.3], xytext=[20, 0.7], rotation=45)
            ax.annotate('Reheating', xy=[950, 0.3], xytext=[70, 0.7], rotation=45)
            ax.annotate('Nucleosynthesis', xy=[950, 0.3], xytext=[150, 0.7], rotation=45)
            ax.annotate('CMB spectral distortions', xy=[950, 0.3], xytext=[300, 0.45], rotation=45)
            ax.annotate('Dark ages', xy=[550, 0.3], xytext=[550, 0.7], rotation=45)
            ax.annotate('Structure formation', xy=[750, 0.3], xytext=[700, 0.7], rotation=45)
            ax.annotate('Reionization', xy=[850, 0.3], xytext=[850, 0.7], rotation=45)
            ax.set_title("History of the Universe")

        ax.set_xlabel('Cosmic Time')
        ax.set_ylabel("Relevance to Today's Talk")

        plt.savefig(os.environ["CMBR_DIST_DATA_PATH"] + '/png/talkintr_%s.png' % eventype)
        plt.close()


def plot_grnf():
    
    fig, ax = plt.subplots(2,1, sharex='all')
    fig.suptitle('CMB distortion visibility', fontsize=25)
    ax[0].plot(gdat.meanreds, gdat.vify, label=r'$\mathcal{J}_y(z_i) = \left(1 + \left(\frac{1+z_i}{6\times10^4}\right)\right)^{-1}$')
    ax[0].plot(gdat.meanreds, gdat.vift, label=r'$\mathcal{J_B}(z_i) = \exp\left(-\frac{z_i}{2\times10^6}\right)^{2.5}$')
    ax[0].plot(gdat.meanreds, gdat.vifm, label=r'$\mathcal{J}_\mu(z_i) = 1 - \exp\left(-\frac{1+z_i}{6\times10^4}\right)^{1.88}$')
    ax[0].plot(gdat.meanreds, gdat.vifm * gdat.vift, label=r'$\mathcal{J}_B(z_i)\mathcal{J}_\mu(z_i)$')
    ax[0].set_xscale('log')
    ax[1].set_xlabel('$z$')
    ax[0].set_ylabel(r'$\mathcal{J}(z)$')
    ax[0].legend(fontsize=16, loc='center left', bbox_to_anchor=[0.,0.5])
    ax[1].plot(gdat.meanreds, gdat.vifm * gdat.vift + gdat.vify - gdat.vift, 'y')
    ax[1].set_xscale('log')
    ax[1].set_ylabel('$\mathcal{J}_B\mathcal{J}_\mu + \mathcal{J}_y - \mathcal{J}_B$')
    plt.savefig(pathplot + 'visifunc.png')
    plt.close()

    with sns.color_palette("Blues", njreds):
        fig, ax = plt.subplots(2, 1, sharex='all')
        for c in range(njreds):
            ax[0].plot(gdat.freqmodl * 1e-9, 1e-6 * fluxfdisgrenconc[:,jreds[c]], label='$z$ = %3.1g' % gdat.meanreds[jreds[c]])
        ax[0].set_xscale('log')
        ax[0].set_title('Full distortion')
        for c in range(njreds):
            ax[1].plot(gdat.freqmodl * 1e-9, 1e-6 * fluxodisgrenconc[:,jreds[c]], label='$z$ = %3.1g' % gdat.meanreds[jreds[c]])
        ax[1].set_xscale('log')
        ax[1].set_xlabel(r'$\nu$ [GHz]')
        ax[1].set_title('Observable distortion')
        ax[1].legend(bbox_to_anchor=[1.1, 1.2], loc='center')
        fig.text(0.05, 0.5, r'$\Delta I_\nu^G$ [MJy/sr]', va='center', rotation='vertical', fontsize=18)
        strg = r'$\Delta I_\nu^G = \mathcal{J}_\mu \mathcal{J}_B I^M(\nu) + \mathcal{J}_y I^Y(\nu) + (1-\mathcal{J}_B) I^T(\nu)$'
        fig.text(0.15, 0.85, strg, va='center', fontsize=18)
        strg = r'$\Delta I_\nu^G = \mathcal{J}_\mu \mathcal{J}_B I^M(\nu) + \mathcal{J}_y I^Y(\nu)$'
        fig.text(0.15, 0.42, strg, va='center', fontsize=18)

        plt.savefig(pathplot + 'fulldist.png')
        plt.close()


def plot_timescal():
    
    timedcom = 1e13 * (gdat.meanreds / 9e3)**(-5.) 
    timebrem = 1e2 * (gdat.meanreds / 1e7)**(-2.5)

    timecomp = 1e-10 * (gdat.meanreds / 6e5)**(-4.)
    timecoul = 1e-3 * (gdat.meanreds / 1e3)**(-1.5)

    timebose = 4e-7 * (gdat.meanreds / 1e7)**(-4.)

    timether0 = 5e12 * (gdat.meanreds / 1e3)**(-4.)
    timether1 = 5e-3 * (gdat.meanreds / 1e7)**(-5)
    timether = minimum(timether0, timether1)

    fig, ax = plt.subplots()
    fig.suptitle('Time scales relevant to Comptonization in the early Universe', fontsize=20)

    #ax.loglog(gdat.meanreds, timebrem, label='Bremsstrahlung')
    #ax.loglog(gdat.meanreds, timedcom, label='Double Compton')
    #ax.loglog(gdat.meanreds, timecomp, label='Compton')
    ax.loglog(gdat.meanreds, timecoul, label='Coulomb')
    ax.loglog(gdat.meanreds, timebose, label='Compton')
    ax.loglog(gdat.meanreds, timether, label='Thermalization')
    ax.loglog(gdat.meanreds, gdat.timehubb / yr2s, color='black', label='Hubble')

    line0 = ax.axvline(5e4, ls='--', color='grey')
    line1 = ax.axvline(2e6, ls='-.', color='grey')

    leg0 = plt.legend([line0, line1], ['Compton Surface', 'Blackbody Surface'], loc='center', bbox_to_anchor=[0.2, 0.1])
    leg1 = ax.legend(loc='center', bbox_to_anchor=[0.65, 0.83])
    ax.add_artist(leg0)

    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$\tau$ [yr]')

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amax(gdat.time[:-1]), amin(gdat.time[:-1])])
    ax_.set_xlabel('$t$ [yr]')

    fig.subplots_adjust(top=0.85)
    plt.savefig(pathplot + 'timescal.png')
    plt.close()


def plot_dist(dist):
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_xlabel(r'$\nu$ [GHz]')
    ax.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')


def plot_pure_dist():
    fig, ax = plt.subplots()
    fig.suptitle('Pure spectral distortions', fontsize=20)
    ax.plot(gdat.freqmodl * 1e-9, fluxcmbrconc * 1e-6, color='black', label=r'$I^B_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{1}{e^x-1}$')
    ax.plot(gdat.freqmodl * 1e-9, fluxtdisgrenconc * 1e-6, label=r'$I^T_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{xe^x}{(e^x-1)^2}$')
    ax.plot(gdat.freqmodl * 1e-9, fluxydisgrenconc * 1e-6, label=r'$I^Y_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{xe^x}{(e^x-1)^2}(x\coth (x/2) - 4)$')
    ax.plot(gdat.freqmodl * 1e-9, fluxmdisgrenconc * 1e-6, label=r'$I^M_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{e^x}{(e^x-1)^2}\left(x/2.2-1\right)$')
    plt.legend(loc=2, fontsize=20)
    ax.set_xscale('log')
    ax.set_xlabel(r'$\nu$ [GHz]')
    ax.set_ylabel(r'$I_\nu$ [MJy/sr]')
    ax_ = ax.twinx()
    ax_.set_ylim(array(ax.get_ylim()) * 1e-20)
    ax_.set_ylabel('[J/m$^2$/s/sr/Hz]')
    ax_ = ax.twiny()
    ax_.set_xlim([amin(sfrqconc), amax(sfrqconc)])
    ax_.set_xlabel('$x$')
    ax_.set_xscale('log')
    plt.savefig(pathplot + 'puredist.png')
    plt.close()


def plot_heat():
    
    heatstrg = ['Silk damping', r's-wave annihilation, $<\sigma v> = 3 \times 10^{-26}$ cm$^3$/s', \
                   r'p-wave annihilation (S-enhanced), $<\sigma v> = 3 \times 10^{-31}$ cm$^3$/s', \
                   r'p-wave annihilation, $<\sigma v> = 3 \times 10^{-36}$ cm$^3$/s', r'Decay, $\tau = 1$ year']
    heatrtag = ['silk', 'swav', 'pwrl', 'pwnr', 'deca']

    csecdmatswav = 3e-20
    csecdmatpwrl = 3e-25
    csecdmatpwnr = 3e-30

    massdmat = 1e11 # [eV]
    timedeca = 1. # [yr]
    ratedeca = 1. / timedeca
    redsdeca = interp1d(gdat.time[::-1], gdat.meanreds[::-1])(timedeca)
    ndendmat = edendmat / massdmat

    ninjetype = 5
    heat = zeros((gdat.nreds, ninjetype))
    heat[:, 0] = 0.1 * edendmat / massdmat * ratedeca
    heat[:, 1] = ndendmat**2 * csecdmatswav
    heat[:, 2] = ndendmat**2 * csecdmatpwrl * (1. + gdat.meanreds)
    heat[:, 3] = ndendmat**2 * csecdmatpwnr * (1. + gdat.meanreds)**2
    heat[:, 4] = 0.05 * ndendmat * ratedeca * exp(-(redsdeca / gdat.meanreds)**2)

    for k in range(ninjetype):
        heat[:, k] *= gdat.timehubb / (1. + gdat.meanreds) / gdat.edenradi

    fig, ax = plt.subplots()
    for k in range(1, ninjetype):
        ax.loglog(gdat.meanreds, gdat.meanreds * heat[:, k], label=heatstrg[k])
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$d\kappa/d\ln z$')
    ax.set_ylim([1e-12, 1e-6])
    ax.legend(loc=2)

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amax(gdat.time[:-1]), amin(gdat.time[:-1])])
    ax_.set_xlabel('$t$ [yr]')

    plt.savefig(pathplot + 'heatrate.png')
    plt.close()
    
    edendmatextd = omegdmat * edencrit * (1. + gdat.meanredsextd)**3
    gdat.edenradiextd = omegradi * edencrit * (1. + gdat.meanredsextd)**4
    ndendmatextd = edendmatextd / massdmat
    gdat.timehubbextd = 4.5e17 * (omegmatt * (1. + gdat.meanredsextd)**3 + omegradi * (1. + gdat.meanredsextd)**4.)**(-0.5)
    timeextd = zeros_like(gdat.meanredsextd)
    for c in range(gdat.nreds-1):
        timeextd[c] = trapz(gdat.timehubbextd[c:] / (1. + gdat.meanredsextd[c:]), gdat.meanredsextd[c:])

    heatextd = zeros((gdat.nreds, ninjetype))
    heatextd[:, 1] = ndendmatextd**2 * csecdmatswav
    heatextd[:, 2] = ndendmatextd**2 * csecdmatpwrl * (1. + gdat.meanredsextd)
    heatextd[:, 3] = ndendmatextd**2 * csecdmatpwnr * (1. + gdat.meanredsextd)**2
    for k in range(1, 4):
        heatextd[:, k] *= gdat.timehubbextd / (1. + gdat.meanredsextd) / gdat.edenradiextd

    fig, ax = plt.subplots()
    for k in range(1, 4):
        ax.loglog(gdat.meanredsextd, gdat.meanredsextd * heatextd[:, k], label=heatstrg[k])
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$d\kappa/d\ln z$')
    #ax.set_ylim([1e-12, 1e-1])
    ax.legend(loc=2)

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amax(timeextd[:-1]), amin(timeextd[:-1])])
    ax_.set_xlabel('$t$ [s]')
    ax_.axvline(1., ls='--', color='grey')
    ax_.annotate('BBN', xy=[1., 1e-12], xytext=[1e2, 1e-12], arrowprops=dict(arrowstyle="->"), fontsize=20)
    plt.savefig(pathplot + 'heatextd.png')
    plt.close()


def plot_deca():
    
    ndeca = 40
    timedecalist = logspace(0., 5., ndeca)

    jdeca = array([0, 9, 19, 29])
    njdeca = 4
    
    dist = zeros((gdat.nreds, ndeca))
    heat = zeros((nreds, ndeca))
    redsdeca = zeros(ndeca)
    for k in range(ndeca):
        disttemp, heattemp, redsdecatemp = retr_fluxdeca(gdat, fluxodisgrenconc, ampldeca, timedecalist[k])
        heatdeca[:, k] = heattemp
        heatdeca[:, k] = dist
        redsdeca[k] = redsdecatemp


    fig, ax = plt.subplots()
    for k in range(njdeca):
        ax.loglog(gdat.meanreds, gdat.meanreds * heatdeca[:, jdeca[k]], label=r'$\tau = %d$ yr' % timedecalist[jdeca[k]])
        ax.set_xlabel(r'$z$')
        ax.set_ylabel(r'$d\kappa/d\ln z$')
        ax.set_ylim([1e-10, 1e-6])
        ax.legend(loc=2, ncol=4)

        ax_ = ax.twiny()
        ax_.set_xscale('log')
        ax_.set_xlim([amax(gdat.time[:-1]), amin(gdat.time[:-1])])
        ax_.set_xlabel('$t$ [yr]')

    plt.savefig(pathplot + 'heatdeca.png')
    plt.close()
    

    timedecaplot = logspace(-4., 6, 10)
    with sns.color_palette("Blues", njreds): 
        fig, ax = plt.subplots()    
        for k in range(timedecaplot.size):
            ax.loglog(gdat.freqmodl * 1e-9, abs(retr_fluxdeca(gdat, fluxodisgren, 1., timedecaplot[k])), label='$t = %.3g$' % timedecaplot[k])
        ax.set_xscale('log')
        ax.set_xlabel(r'$\nu$ [GHz]')
        ax.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')
        ax.legend(loc=4)
        plt.close()


def plot_sampdist():
    diffdistdiffreds = heat[None, :, :] * fluxodisgren[:, :, None]
    dist = zeros((nfreq, ninjetype))
    for k in range(ninjetype):
        dist[:, k] = trapz(diffdistdiffreds[:, :, k], gdat.meanreds, axis=1)

    fig, ax = plt.subplots()
    for k in range(1, ninjetype):
        ax.plot(gdat.freqmodl * 1e-9, dist[:, k], label=heatstrg[k])
    ax.set_xscale('log')
    ax.set_xlabel(r'$\nu$ [GHz]')
    ax.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')
    ax.legend(loc=2)
    ax.axhline(5., ls='--', color='grey')
    ax.axhline(-5., ls='--', color='grey')
    ax.fill_between(gdat.freqmodl * 1e-9, ones_like(gdat.freqmodl) * 5., ones_like(gdat.freqmodl) * -5., color='grey')
    ax.text(2, 10, r'PIXIE 1-$\sigma$ sensitivity', color='grey', fontsize=20)
    plt.savefig(pathplot + 'totldist.png')
    plt.close()
    
    
    diffdistdifflred = heat[None, :, :] * fluxodisgren[:, :, None] * gdat.meanreds[None, :, None]
    for k in range(1, 2):
        ylim = [amin(diffdistdifflred[:, :, k]), amax(diffdistdifflred[:, :, k])]
        for c in range(njreds):
            fig, ax = plt.subplots()
            ax.plot(gdat.freqmodl * 1e-9, diffdistdifflred[:, numbreds-1-jreds[c], k])
            text = plt.figtext(0.2, 0.8, '$z_i$ = %.3g' % gdat.meanreds[nreds-jreds[c]-1], fontsize=20)

            ax.set_xscale('log')
            ax.set_xlabel(r'$\nu$ [GHz]')
            ax.set_ylabel(r'$d\Delta I_\nu/ d\ln z$ [Jy/sr]')
            ax.set_ylim(ylim)
            ax.set_title(heatstrg[k])

            plt.legend(loc=2, ncol=2)
            plt.savefig(pathplot + 'diffdistdifflred_' + heatrtag[k] + '_%d.png' % jreds[c])
            plt.close(fig)
        

def retr_fluxdeca(gdat, fluxodisgren, ampldeca, timedeca):

    ratedeca = 1. / timedeca
    redsdeca = interp1d(gdat.time, gdat.meanreds)(timedeca)
    
    heatdeca = ampldeca * gdat.ndenbmat * ratedeca * exp(-timedeca / gdat.time) * gdat.timehubb / (1. + gdat.meanreds) / gdat.edenradi

    difffluxdecadiffreds = heatdeca[None, :] * fluxodisgren
    fluxdeca = trapz(difffluxdecadiffreds, gdat.meanreds, axis=1)

    return fluxdeca


def plot_llik():
    lliktopo = zeros((numbbins + 1, numbbins + 1, numbbins + 1, numbbins + 1))
    thissamp = zeros(numbpara)
    for a in range(numbbins + 1):
        for b in range(numbbins + 1):
            for c in range(numbbins + 1):
                for d in range(numbbins + 1):
                    thissamp[0] = binstimedeca[a]
                    thissamp[1] = binsampldeca[b]
                    thissamp[2] = binstempcmbr[c]
                    thissamp[3] = binsydisampl[d]
                    lliktopo[a, b, c, d] = retr_lpos(thissamp)
    
    lliktopo = trapz(trapz(lliktopo, binsydisampl, 3), binstempcmbr, 2)
    
    fig, ax = plt.subplots()
    imag = ax.imshow(lliktopo, origin='lower', interpolation='none', cmap='Reds', extent=[minmtimedeca, maxmtimedeca, minmampldeca, maxmampldeca])
    
    sp.stats.chi2.ppf(1 - 0.05, 1) / 2.
    
    levl = zeros(2)
    levl[0] = array([amax(lliktopo) - sp.stats.chi2.ppf(1 - 0.37, 1) / 2.])
    levl[1] = array([amax(lliktopo) - sp.stats.chi2.ppf(1 - 0.05, 1) / 2.])
    cont = ax.contour(binstimedeca, binsampldeca, lliktopo, origin='lower', color='b', levels=levl)
    
    ax.set_ylabel("$f_X$ [eV]")
    ax.set_xlabel(r"$\tau_X$ [year]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    #plt.colorbar(imag, axis=axis)
    plt.close() 


def retr_occp(sfrq, dist=False):

    occpplnk = 1. / (exp(sfrq) - 1.)
    
    if dist:
        occptdis = sfrq * exp(sfrq) / (exp(sfrq) - 1.)**2
        occpydis = sfrq * exp(sfrq) / (exp(sfrq) - 1.)**2 * (sfrq / tanh(sfrq / 2.) - 4.)
        occpmdis = exp(sfrq) / (exp(sfrq) - 1.)**2 * (sfrq / 2.2 - 1)
        return occpplnk, occptdis, occpydis, occpmdis
    else:
        return occpplnk


def retr_fluxcmbr(gdat, thisfreq, tempcmbr):
    
    sfrq = gdat.plnkcons * thisfreq / gdat.boltcons / tempcmbr
    
    occpplnk = retr_occp(sfrq)

    fluxcmbr = 2. * gdat.plnkcons * thisfreq**3 / gdat.velolght**2 * occpplnk * 1e26
    
    return fluxcmbr


def retr_fluxgren(gdat, thisfreq, tempcmbr, disttype):
    
    sfrq = gdat.plnkcons * thisfreq / gdat.boltcons / tempcmbr

    occpplnk, occptdis, occpydis, occpmdis = retr_occp(sfrq, dist=True)
    
    fluxtdisgren = 0.25 * 2. * gdat.plnkcons * thisfreq**3 / gdat.velolght**2 * occptdis * 1e26
    fluxydisgren = 0.25 * 2. * gdat.plnkcons * thisfreq**3 / gdat.velolght**2 * occpydis * 1e26
    fluxmdisgren = 1.41 * 2. * gdat.plnkcons * thisfreq**3 / gdat.velolght**2 * occpmdis * 1e26
    
    fluxfdisgren = gdat.vifm[None, :] * gdat.vift[None, :] * fluxmdisgren[:, None] + gdat.vify[None, :] * fluxydisgren[:, None] + (1. - gdat.vift[None, :]) * fluxtdisgren[:, None]
    fluxodisgren = gdat.vifm[None, :] * gdat.vift[None, :] * fluxmdisgren[:, None] + gdat.vify[None, :] * fluxydisgren[:, None]
    
    if disttype == 'full':
        return fluxtdisgren, fluxydisgren, fluxmdisgren, fluxfdisgren, fluxodisgren
    elif disttype == 'ydisodis':
        return fluxydisgren, fluxodisgren


def retr_fluxsync(thisfreq, syncnorm, syncindx):
    
    fluxsync = syncnorm * (thisfreq / 1e9)**syncindx
    
    return fluxsync
  
    
def retr_fluxfree(gdat, thisfreq, emmefree, tempfree):
    
    thisplnkfunc = retr_plnkfunc(gdat, thisfreq, tempfree) * 1e26
    
    #print 'thisplnkfunc'
    #print thisplnkfunc
    #print 'emmefree / thisfreq**2.1 / tempfree**1.5'
    #print emmefree / thisfreq**2.1 / tempfree**1.5
    
    gaunfact = log(4.955e2 * (thisfreq / 1e9)**(-1.)) + 1.5 * log(tempfree)
    odepfree = 3.014e-2 * (tempfree)**(-1.5) * (thisfreq / 1e9)**(-2.) * emmefree * gaunfact
    fluxfree = (1. - exp(-odepfree)) * thisplnkfunc
    
    return fluxfree


def retr_plnkfunc(gdat, thisfreq, temp):
    
    thissfrq = gdat.plnkcons * thisfreq / gdat.boltcons / temp
    
    plnk = 2. * gdat.plnkcons * thisfreq**3 / gdat.velolght**2 / (exp(thissfrq) - 1.)
    
    return plnk


def plot_temp():
    
    #def difftempdmatdiffreds(tempdmat, thisreds):
    #    thishubb = interp1d(gdat.meanreds, hubb)(thisreds)
    #    difftempdmatdiffreds = -2. * thishubb * tempdmat# + gamm * b * (tempbmat - tempdmat)
    #    return difftempdmatdiffreds
    #inittempdmat = array([3000.])
    #tempdmat = odeint(difftempdmatdiffreds, inittempdmat, gdat.meanreds)

    tempmatt = gdat.tempcmbrconc * (1. + gdat.meanreds) / (1. + 119.  / (1. + gdat.meanreds) / (1. + ((1. + gdat.meanreds) / 115.)**1.5))
    tempdmatcold = gdat.tempcmbrconc * (1. + gdat.meanreds) / (1. + 1e9  / (1. + gdat.meanreds) / (1. + ((1. + gdat.meanreds) / 1e9)**2.5))
    tempdmatwarm = gdat.tempcmbrconc * (1. + gdat.meanreds) / (1. + 5e6  / (1. + gdat.meanreds) / (1. + ((1. + gdat.meanreds) / 5e6)**2.5))
    tempcmbr = gdat.tempcmbrconc * (1. + gdat.meanreds)

    fig, ax = plt.subplots()
    ax.loglog(gdat.meanreds, tempcmbr, label='CMB')
    ax.loglog(gdat.meanreds, tempdmatwarm, label='WDM')
    ax.loglog(gdat.meanreds, tempdmatcold, label='CDM')
    ax.loglog(gdat.meanreds, tempmatt, label='Baryonic matter')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'T(z) [K]')

    enerradilate = gdat.tempcmbrconc * (1. + gdat.meanreds) * gdat.boltconsnatu

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amin(enerradilate), amax(enerradilate)])
    ax_.set_xlabel(r'$E_\gamma$ [eV]')

    ax.legend(loc=2)
    plt.savefig(pathplot + 'tempevol.png')
    plt.close()
    

def plot_silkscal():
    
    wnumsilk = (5.92e10)**(-0.5) * (1. + gdat.meanreds)**1.5
    wlensilk = 2. * pi / wnumsilk

    wlenahor = (1. + gdat.meanreds) * gdat.time * yr2s / mp2m * gdat.velolght / sqrt(3. * (1. + edenbmat / gdat.edenradi))
    wnumahor = 2. * pi / wlenahor

    wlenphor = (1. + gdat.meanreds) * gdat.time * yr2s / mp2m * gdat.velolght
    wnumphor = 2. * pi / wlenphor

    fig, ax = plt.subplots()
    ax.loglog(gdat.meanreds, wnumsilk, label='Dissipation scale')
    ax.set_ylabel('$k$ [Mpc$^{-1}$]')
    ax.set_xlabel('$z$')
    ax.axvline(5e4, ls='--', color='black')
    ax.axvline(2e6, ls='--', color='black')
    ax.axhline(interp1d(gdat.meanreds, wnumsilk)(5e4), ls='--', color='black')
    ax.axhline(interp1d(gdat.meanreds, wnumsilk)(2e6), ls='--', color='black')

    plt.figtext(0.19, 0.2, 'y-type distortion', fontsize=20)
    plt.figtext(0.52, 0.2, r'$\mu$-type distortion', fontsize=20)
    plt.figtext(0.78, 0.2, 'Black body', fontsize=20)
    ax_ = ax.twinx()
    ax.set_ylim([amin(1e-3), amax(wnumsilk)])
    ax_.set_ylim([2. * pi / 1e-3, amin(wlensilk)])
    ax_.set_yscale('log')
    ax_.set_ylabel(r'$\lambda$ [Mpc]')
    ax.set_title('Silk dissipation scale', fontsize=20)

    plt.savefig(pathplot + 'wlensilk.png')
    fig.subplots_adjust(top=0.9)
    plt.close()


def plot_resi_post(gdat, postflux):

    fig = plt.figure(figsize=(16, 12))
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 3]) 
    ax = [] 
    ax.append(plt.subplot(gs[0]))
    ax.append(plt.subplot(gs[1], sharex=ax[0]))

    ax[1].set_title('')
    ax[1].errorbar(freqexpr * 1e-9, gdat.dataflux * 1e-6, ls='none', yerr=gdat.datafluxstdv*1e-6, xerr=gdat.datafluxstdv*1e-9, label=gdat.datalabl, marker='o', markersize=5, color='k')
    
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxtotl * 1e-6, lcol='salmon', alpha=0.5, dcol='red', mcol='black')
    if gdat.inclcmbrmono:
        tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxcmbr * 1e-6, lcol='lightblue', alpha=0.5, dcol='blue', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxdustcold * 1e-6, lcol='lightgreen', alpha=0.5, dcol='green', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxdustwarm * 1e-6, lcol='lightgreen', alpha=0.5, dcol='green', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxsync * 1e-6, lcol='lightyellow', alpha=0.5, dcol='yellow', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxfree * 1e-6, lcol='lightcoral', alpha=0.5, dcol='coral', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxydis * 1e-6, lcol='lightyellow', alpha=0.5, dcol='yellow', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], gdat.freqmodl * 1e-9, listfluxdeca * 1e-6, lcol='lightcyan', alpha=0.5, dcol='cyan', mcol='black')
    
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel(r'$\nu$ [GHz]')
    ax[1].set_ylabel(r'$I_\nu$ [MJy/sr]')
    ax[1].set_ylim([1e-7, 1e4])
    ax[1].legend(loc=9, ncol=4)

    tdpy.mcmc.plot_braz(ax[0], freqexpr * 1e-9, listresiflux, lcol='lightgrey',    alpha=0.5, dcol='darkgrey', mcol='black')
    ax[0].set_ylabel(r'$I_\nu^{res}$ [Jy/sr]')
    ax[0].axhline(0., ls='--', color='black', alpha=0.1)

    if path == None:
        plt.show()
    else:
        plt.savefig(path)
        plt.close(fig)


def plot_dataflux(gdat):
        
    figr = plt.figure(figsize=(16, 12))
    
    xlim = [gdat.minmfreqmodl * 1e-9, gdat.maxmfreqmodl * 1e-9]

    if gdat.datatype == 'expr':
        axis = [plt.gca()]
    if gdat.datatype == 'mock':
        axisgridtemp = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 3]) 
        axis = [] 
        axis.append(plt.subplot(axisgridtemp[0]))
        axis.append(plt.subplot(axisgridtemp[1], sharex=axis[0]))

    yerr = gdat.datafluxstdv * 1e-6
    xerr = gdat.datafluxstdv * 1e-9
    axis[0].errorbar(gdat.freqexpr * 1e-9, gdat.dataflux * 1e-6, ls='none', xerr=xerr, yerr=yerr, label=gdat.datalabl, marker='o', markersize=5, color='k')
    
    if gdat.datatype == 'mock':
        axis[0].plot(gdat.freqmodl * 1e-9, modlfluxtotl * 1e-6, label='Total model', color='r')
        if gdat.inclcmbrmono:
            axis[0].plot(gdat.freqmodl * 1e-9, modlfluxcmbr * 1e-6, label='CMB', color='b')
        axis[0].plot(gdat.freqmodl * 1e-9, modlfluxdustcold * 1e-6, label='Cold Dust', color='g', ls='--')
        axis[0].plot(gdat.freqmodl * 1e-9, modlfluxdustwarm * 1e-6, label='Warm Dust', color='g', ls='-.')
        axis[0].plot(gdat.freqmodl * 1e-9, modlfluxsync * 1e-6, label='Synchrotron', color='y')
        axis[0].plot(gdat.freqmodl * 1e-9, modlfluxfree * 1e-6, label='Brem', color='coral')
        axis[0].plot(gdat.freqmodl * 1e-9, abs(modlfluxydis) * 1e-6, label='Reionization', color='m')
        axis[0].plot(gdat.freqmodl * 1e-9, abs(modlfluxdeca) * 1e-6, label='Particle decay', color='cyan')
    
    axis[0].set_title('')
    axis[0].set_xscale('log')
    axis[0].set_yscale('log')
    axis[0].set_xlabel(r'$\nu$ [GHz]')
    axis[0].set_ylabel(r'$I_\nu$ [MJy/sr]')
    axis[0].set_ylim([1e-7, 1e4])
    axis[0].legend(loc=9, ncol=4)
    axis[0].set_xlim(xlim)

    if gdat.datatype == 'mock':
        axis[1].errorbar(freqexpr * 1e-9, resiflux, yerr=gdat.datafluxstdv, xerr=freqexprstdv*1e-9, marker='o', lw=1, ls='none', markersize=5, color='k')
        axis[1].set_ylabel(r'$I_\nu^{res}$ [Jy/sr]')
        axis[1].axhline(0., ls='--', color='black', alpha=0.1)
        axis[1].set_xlim(xlim)


    plt.savefig(path)
    plt.close(fig)


def retr_ydis_trac():
    
    path = pathplot + 'electron_photonspectrum_results.fits'
    data, hdr = pf.getdata(path, 1, header=True)

    redsoutp = data['OUTPUT_REDSHIFT'].squeeze()
    enerpart = 10**data['LOG10ENERGY'].squeeze()
    redsinpt = data['INPUT_REDSHIFT'].squeeze()
    enerphot = data['PHOTON_ENERGY'].squeeze()
    freqphot = enerphot / gdat.plnkcons

    nredsoutp = redsoutp.size
    nenerpart = enerpart.size
    nredsinpt = redsinpt.size
    nfreqphot = freqphot.shape[0]

    njredsoutp = 5
    njenerpart = 1
    njredsinpt = 5
    jredsoutp = [k * nredsoutp / njredsoutp for k in range(njredsoutp)]
    jenerpart = [k * nenerpart / njenerpart for k in range(njenerpart)]
    jredsinpt = [k * nredsinpt / njredsinpt for k in range(njredsinpt)]

    specchek = 8. * pi * enerphot[:,:,None]**2 / (gdat.velolght * gdat.plnkcons)**3 / (exp(enerphot[:,:,None] / tempcmbrnunc / gdat.boltcons / (1. + redsinpt[None,None,:])) - 1.) # [1/cm^3/eV]

    diffydisdiffreds = data['YDISTORTION'].squeeze()
    ydis = data['CUMULATIVE_YDISTORTION'].squeeze()
    spec = data['CMB_SPECTRUM'].squeeze()
    nredsout = 63
    nredsinp = 63
    nenerpart = 40
    nfreqphot = 500

    for a in range(njenerpart):
        for f in range(njredsinpt):
            labl = '$E_{inj} = %.3g, z_{inp} = %.3g$' % (enerpart[jenerpart[a]], redsinpt[jredsinpt[f]])
            plt.plot(freqphot[:,jenerpart[a]] * 1e-9, spec[:,jenerpart[a],jredsinpt[f]] * enerphot[:,jenerpart[a]], label=labl)
            plt.plot(freqphot[:,jenerpart[a]] * 1e-9, specchek[:,jenerpart[a],jredsinpt[f]] * enerphot[:,jenerpart[a]], '--')
        plt.xlabel(r'$\nu_\gamma$ [GHz]')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc=10, bbox_to_anchor=[1.3,0.5])
        plt.show()

    plt.figure(figsize=(12,2))
    plt.text(0.5, 0.5,'Differential y-distortion', ha='center', va='center')
    plt.axis('off')
    plt.show()

    for a in range(njenerpart):
        for c in range(njredsoutp):
            for f in range(njredsinpt):
                labl = '$z_{out} = %.3g, E_{inj} = %.3g, z_{inp} = %.3g$' %     (redsoutp[jredsoutp[c]], enerpart[jenerpart[a]], redsinpt[jredsinpt[f]])
                plt.plot(freqphot[:,jenerpart[a]] * 1e-9, diffydisdiffreds[:,jredsinpt[f],jenerpart[a],jredsoutp[c]], label=labl)
            plt.xlabel(r'$\nu_\gamma$ [GHz]')
            plt.xscale('log')
            plt.legend(loc=10, bbox_to_anchor=[1.3,0.5])
            plt.show()

    plt.figure(figsize=(12,2))
    plt.text(0.5, 0.5,'Cumulative y-distortion', ha='center', va='center')
    plt.axis('off')
    plt.show()

    for a in range(njenerpart):
        for c in range(njredsoutp):
            for f in range(njredsinpt):
                labl = '$z_{out} = %.3g, E_{inj} = %.3g, z_{inp} = %.3g$' % (redsoutp[jredsoutp[c]], enerpart[jenerpart[a]], redsinpt[jredsinpt[f]])
                plt.plot(freqphot[:,jenerpart[a]] * 1e-9, ydis[:,jredsinpt[f],jenerpart[a],jredsoutp[c]], label=labl)
            plt.xlabel(r'$\nu_\gamma$ [GHz]')
            plt.xscale('log')
            plt.legend(loc=10, bbox_to_anchor=[1.3,0.5])
            plt.show()

    fig, ax = plt.subplots(figsize=(10,10));
    line, = ax.plot(freqphot[:,0], spec[:,0,c])
    plt.close()
    ax.set_xscale('log')
    ax.set_yscale('log')


def plot_cmbrtemppsec():
    
    nside = 256
    lghp, bghp, nside, npixl, apix = tdpy.util.retr_heal(nside) 

    fname = os.environ["CMBR_DIST_DATA_PATH"] + '/dat/COM_PowerSpect_CMB-TT-loL-full_R2.01.txt'
    plnkdata = loadtxt(fname)
    lmodloww = plnkdata[:, 0]
    clhploww = plnkdata[:, 1]
    fname = os.environ["CMBR_DIST_DATA_PATH"] + '/dat/COM_PowerSpect_CMB-TT-hiL-full_R2.01.txt'
    plnkdata = loadtxt(fname)
    lmodhigh = plnkdata[:, 0]
    clhphigh = plnkdata[:, 1]

    flux = hp.synfast(clhphigh, nside)
    fluxcart = tdpy.util.retr_cart(flux)
    fig, ax = plt.subplots()
    imag = ax.imshow(fluxcart, origin='lower', cmap='Reds')
    plt.colorbar(imag, ax=axis, fraction=0.05)

    flux = hp.synfast(clhploww, nside)
    fluxcart = tdpy.util.retr_cart(flux)
    fig, ax = plt.subplots()
    imag = ax.imshow(fluxcart, origin='lower', cmap='Reds')
    plt.colorbar(imag, ax=axis, fraction=0.05)

    fig, ax = plt.subplots()
    ax.plot(lmodhigh, clhphigh)
    ax.plot(lmodloww, clhploww)               


def retr_flux(gdat, thisfreq, thispara):

    tempcmbr = thispara[0]
    dustodep = thispara[1]
    dustemisrati = thispara[2]
    dustpowrfrac = thispara[3]
    dustwarmindx = thispara[4]
    dustwarmtemp = thispara[5]
    dustcoldindx = thispara[6]
    syncnorm = thispara[7]
    syncindx = thispara[8]
    emmefree = thispara[9]
    tempfree = thispara[10]
    ydisampl = thispara[11]
    ampldeca = thispara[12]
    timedeca = thispara[13]

    if gdat.verbtype > 1:
        print 'tempcmbr: ', tempcmbr
        print 'dustodep: ', dustodep
        print 'dustemisrati: ', dustemisrati
        print 'dustpowrfrac: ', dustpowrfrac
        print 'dustwarmindx: ', dustwarmindx
        print 'dustwarmtemp: ', dustwarmtemp
        print 'dustcoldindx: ', dustcoldindx
        print 'syncnorm: ', syncnorm
        print 'syncindx: ', syncindx
        print 'emmefree: ', emmefree
        print 'tempfree: ', tempfree
        print 'ydisampl: ', ydisampl
        print 'ampldeca: ', ampldeca
        print 'timedeca: ', timedeca
        
    # CMB
    fluxcmbr = retr_fluxcmbr(gdat, thisfreq, tempcmbr)

    # dust
    fluxdust, fluxdustcold, fluxdustwarm = retr_fluxdust(gdat, thisfreq, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx)
   
    # synchrotron
    fluxsync = retr_fluxsync(thisfreq, syncnorm, syncindx)
    
    # Bremsstrahlung
    fluxfree = retr_fluxfree(gdat, thisfreq, emmefree, tempfree)
        
    fluxydisgren, fluxodisgren = retr_fluxgren(gdat, thisfreq, tempcmbr, disttype='ydisodis')
    
    # free y-distortion
    fluxydis = ydisampl * fluxydisgren
    
    # Particle decay
    fluxdeca = retr_fluxdeca(gdat, fluxodisgren, ampldeca, timedeca)
    
    fluxtotl = fluxdust + fluxsync + fluxfree + fluxydis + fluxdeca
    if gdat.inclcmbrmono:
        fluxtotl += fluxcmbr

    return fluxtotl, fluxcmbr, fluxdustcold, fluxdustwarm, fluxsync, fluxfree, fluxydis, fluxdeca


def retr_fluxdust(gdat, thisfreq, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx):

    factwarm = (sp.special.zetac(4. + dustwarmindx) + 1) * sp.special.gamma(4. + dustwarmindx)
    factcold = (sp.special.zetac(4. + dustcoldindx) + 1) * sp.special.gamma(4. + dustcoldindx)
    dustcoldtemp = (factwarm / factcold / dustemisrati * (gdat.plnkcons * gdat.freqpivt / gdat.boltcons)**(dustcoldindx - dustwarmindx) * \
            dustwarmtemp**(4. + dustwarmindx))**(1. / (4. + dustcoldindx))
    
    fluxdustcoldfact = dustpowrfrac * dustemisrati * (thisfreq / 3e12)**dustcoldindx
    fluxdustwarmfact = (1. - dustpowrfrac) * (thisfreq / 3e12)**dustwarmindx
        
    fluxdustcold = dustodep * fluxdustcoldfact * retr_plnkfunc(gdat, thisfreq, dustcoldtemp) * 1e26 / (fluxdustcoldfact + fluxdustwarmfact)
        
    fluxdustwarm = dustodep * fluxdustwarmfact * retr_plnkfunc(gdat, thisfreq, dustwarmtemp) * 1e26 / (fluxdustcoldfact + fluxdustwarmfact)
        
    fluxdust = fluxdustcold + fluxdustwarm
    
    return fluxdust, fluxdustcold, fluxdustwarm


def retr_llik(sampvarb, gdat):
    
    modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, modlfluxfree, modlfluxydis, modlfluxdeca = retr_flux(gdat, gdat.freqmodl, sampvarb)
       
    modlfluxintp = interp1d(gdat.freqmodl, modlfluxtotl)(freqexpr)
    resiflux = gdat.dataflux - modlfluxintp

    llik = sum(-log(sqrt(2. * pi) * gdat.datafluxstdv) - 0.5 * (modlfluxintp - gdat.dataflux) / gdat.datafluxstdv**2 * (modlfluxintp - gdat.dataflux))

    if savepost:
        sampcalc = modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, modlfluxfree, modlfluxydis, modlfluxdeca, resiflux
    else:
        sampcalc = []

    return llik, sampcalc


def retr_mocksampvarb():
    
    numbpara = 14
    
    mocksampvarb = empty(numbpara)
    mocksampvarb[0] = 2.725
    mocksampvarb[1] = 1e-3
    mocksampvarb[2] = 1e1
    mocksampvarb[3] = 0.025
    mocksampvarb[4] = 2.75
    mocksampvarb[5] = 15.
    mocksampvarb[6] = 1.75
    mocksampvarb[7] = 1e5
    mocksampvarb[8] = -0.75
    mocksampvarb[9] = 1e2
    mocksampvarb[10] = 1e2
    mocksampvarb[11] = 1e-7
    mocksampvarb[12] = 1e-5
    mocksampvarb[13] = 1e10
    
    return mocksampvarb


def retr_datapara():
    
    numbpara = 14
   
    datapara = tdpy.util.gdatstrt()

    datapara.indx = dict()
    datapara.minm = zeros(numbpara)
    datapara.maxm = zeros(numbpara)
    datapara.name = empty(numbpara, dtype=object)
    datapara.scal = empty(numbpara, dtype=object)
    datapara.labl = empty(numbpara, dtype=object)
    datapara.unit = empty(numbpara, dtype=object)
    datapara.vari = zeros(numbpara)
    
    datapara.indx['tempcmbr'] = 0
    datapara.name[0] = 'tempcmbr'
    datapara.minm[0] = 2.72
    datapara.maxm[0] = 2.73
    datapara.scal[0] = 'self'
    datapara.labl[0] = '$T_{cmb}$'
    datapara.unit[0] = '[K]'
    datapara.vari[0] = 1e-6

    datapara.indx['dustodep'] = 1
    datapara.name[1] = 'dustodep'
    datapara.minm[1] = 1e-4
    datapara.maxm[1] = 1e-2
    datapara.scal[1] = 'logt'
    datapara.labl[1] = r'$\tau_{dust}$'
    datapara.unit[1] = ''
    datapara.vari[1] = 1e-8

    datapara.indx['dustemisrati'] = 2
    datapara.name[2] = 'dustemisrati'
    datapara.minm[2] = 1e0
    datapara.maxm[2] = 1e2
    datapara.scal[2] = 'logt'
    datapara.labl[2] = '$q_1/q_2$'
    datapara.unit[2] = ''
    datapara.vari[2] = 1e-6
    
    datapara.indx['dustpowrfrac'] = 3
    datapara.name[3] = 'dustpowrfrac'
    datapara.minm[3] = 0.
    datapara.maxm[3] = 0.05
    datapara.scal[3] = 'self'
    datapara.labl[3] = '$f_1$'
    datapara.unit[3] = ''
    datapara.vari[3] = 4e-6
    
    datapara.indx['dustwarmindx'] = 4
    datapara.name[4] = 'dustwarmindx'
    datapara.minm[4] = 2.5
    datapara.maxm[4] = 3.
    datapara.scal[4] = 'self'
    datapara.labl[4] = r'$\beta_2$'
    datapara.unit[4] = ''
    datapara.vari[4] = 1e-6
    
    datapara.indx['dustwarmtemp'] = 5
    datapara.name[5] = 'dustwarmtemp'
    datapara.minm[5] = 10.
    datapara.maxm[5] = 20.
    datapara.scal[5] = 'self'
    datapara.labl[5] = '$T_2$'
    datapara.unit[5] = '[K]'
    datapara.vari[5] = 5e-8
    
    datapara.indx['dustcoldindx'] = 6
    datapara.name[6] = 'dustcoldindx'
    datapara.minm[6] = 1.
    datapara.maxm[6] = 2.5
    datapara.scal[6] = 'self'
    datapara.labl[6] = r'$\beta_1$'
    datapara.unit[6] = ''
    datapara.vari[6] = 8e-6
    
    datapara.indx['syncnorm'] = 7
    datapara.name[7] = 'syncnorm'
    datapara.minm[7] = 1e3
    datapara.maxm[7] = 1e7
    datapara.scal[7] = 'logt'
    datapara.labl[7] = '$A_{sync}$'
    datapara.unit[7] = ''
    datapara.vari[7] = 2e-6
    
    datapara.indx['syncindx'] = 8
    datapara.name[8] = 'syncindx'
    datapara.minm[8] = -1.5
    datapara.maxm[8] = 0.
    datapara.scal[8] = 'self'
    datapara.labl[8] = r'$\alpha_{sync}$'
    datapara.unit[8] = ''
    datapara.vari[8] = 8e-6
    
    datapara.indx['emmefree'] = 9
    datapara.name[9] = 'emmefree'
    datapara.minm[9] = 1e0
    datapara.maxm[9] = 1e4
    datapara.scal[9] = 'logt'
    datapara.labl[9] = 'EM'
    datapara.unit[9] = '[pc/cm$^6$]'
    datapara.vari[9] = 5e-5
    
    datapara.indx['tempfree'] = 10
    datapara.name[10] = 'tempfree'
    datapara.minm[10] = 1e1
    datapara.maxm[10] = 1e3
    datapara.scal[10] = 'logt'
    datapara.labl[10] = r'$T_e$'
    datapara.unit[10] = '[K]'
    datapara.vari[10] = 1e-3
    
    datapara.indx['ydisampl'] = 11
    datapara.name[11] = 'ydisampl'
    datapara.minm[11] = 1e-9
    datapara.maxm[11] = 1e-5
    datapara.scal[11] = 'logt'
    datapara.labl[11] = '$y_{ri}$'
    datapara.unit[11] = ''
    datapara.vari[11] = 1e-3
    
    datapara.indx['ampldeca'] = 12
    datapara.name[12] = 'ampldeca'
    datapara.minm[12] = 1e-7
    datapara.maxm[12] = 1e-3
    datapara.scal[12] = 'logt'
    datapara.labl[12] = '$f_X$'
    datapara.unit[12] = '[eV]'
    datapara.vari[12] = 2e-1
    
    datapara.indx['timedeca'] = 13
    datapara.name[13] = 'timedeca'
    datapara.minm[13] = 1e8
    datapara.maxm[13] = 1e12
    datapara.scal[13] = 'logt'
    datapara.labl[13] = r'$\tau_X$'
    datapara.unit[13] = '[s]'
    datapara.vari[13] = 2e-1
    
    datapara.strg = datapara.labl + ' ' + datapara.unit
    
    return datapara


def retr_egbl(thisfreq):
    
    path = os.environ['CMBR_DIST_DATA_PATH'] + '/egbl.csv'
    egbldata = loadtxt(path)
    wlenegbl = egbldata[:, 0] * 1e-6 # [m]
    freqegbl = flipud(gdat.velolght / wlenegbl) # [Hz]
    
    minmfreqegbl = amin(freqegbl)
    maxmfreqegbl = amax(freqegbl)
    
    fluxegbltemp = flipud(egbldata[:, 1]) * 1e26 / freqegbl # [Jy/sr]
    
    fluxegbl = zeros_like(thisfreq)
    jthisfreq = where((thisfreq < maxmfreqegbl) & (minmfreqegbl < thisfreq))[0]
    fluxegbl[jthisfreq] = interp1d(freqegbl, fluxegbltemp)(thisfreq[jthisfreq])

    return fluxegbl
        

def init( \
         datatype='mock', \
         exprtype='pixi', \
         datalabl='PIXIE', \
         numbswep=100000, \
         numbburn=0, \
         factthin=1, \
         exprflux=None, \
         exprfluxstdv=None, \
         freqexpr=None, \
         numbfreqexpr=None, \
         freqexprstdv=None, \
         minmfreqexpr=None, \
         maxmfreqexpr=None, \
         exprfluxstdvinst=None, \
         exprfluxstdvfrac=None, \
         inclcmbrmono=True, \
         plotperd=10000, \
         verbtype=1, \
         optiprop=False, \
         makeplot=True, \
        ):
    
    gdat = tdpy.util.gdatstrt()
   
    gdat.datatype = datatype
    gdat.exprtype = exprtype
    gdat.datalabl = datalabl
    
    gdat.verbtype = verbtype
    
    gdat.numbfreqexpr = numbfreqexpr
    gdat.freqexprstdv = freqexprstdv
    gdat.minmfreqexpr = minmfreqexpr
    gdat.maxmfreqexpr = maxmfreqexpr
    gdat.exprfluxstdvinst = exprfluxstdvinst
    gdat.exprfluxstdvfrac = exprfluxstdvfrac

    gdat.inclcmbrmono = inclcmbrmono

    rtag = retr_rtag(gdat)

    pathbase = os.environ["CMBR_DIST_DATA_PATH"]
    pathplot = pathbase + '/png/' + rtag + '/'
    cmnd = 'mkdir -p ' + pathplot
    os.system(cmnd)
    
    if freqexpr == None:
        gdat.freqexpr = logspace(log10(gdat.minmfreqexpr), log10(gdat.maxmfreqexpr), gdat.numbfreqexpr) # [Hz]
    
    numbfreqmodl = 1000
    gdat.minmfreqmodl = 1e9
    gdat.maxmfreqmodl = 1e13
    gdat.freqmodl = logspace(log10(gdat.minmfreqmodl), log10(gdat.maxmfreqmodl), numbfreqmodl) # [Hz]
    
    # physical constants
    gdat.velolght = 3e8 # [m/s]
    mp2m = 3.1e22 # [Mpc/m]
    yr2s = 364. * 3600. * 24. # [year/s]
    gdat.plnkcons = 6.63e-34 # [J s]
    gdat.boltcons = 1.38e-23 # [J/K]
    gdat.freqpivt = 1e12 # [Hz]
    gdat.tempcmbrconc = 2.725 # Planck concordance model temperature of the CMB today [K]
    gdat.boltconsnatu = 8.6173e-5 # Boltzmann constant [eV/K]
    massprot = 9.38e8 # [eV]
    
    # distortion related constants
    alph = 1.401
    redsdism = 2e6

    # temp
    savepost = False
    
    # redshift axis
    numbreds = 100
    minmreds = 1e2
    maxmreds = 1e10
    gdat.meanreds = logspace(log10(minmreds), log10(maxmreds), numbreds)
    numbredsplot = 10
    indxredsplot = [k * numbreds / numbredsplot for k in range(numbredsplot)]
    diffreds = gdat.meanreds[1:] - gdat.meanreds[:-1]
    gdat.meanredsextd = logspace(2., 10., numbreds)
    
    # time
    gdat.timehubb = 4.5e17 * (0.27 * (1. + gdat.meanreds)**3 + 9.2e-5 * (1. + gdat.meanreds)**4.)**(-0.5)
    gdat.time = zeros(numbreds)
    for c in range(numbreds-1):
        gdat.time[c] = trapz(gdat.timehubb[c:] / (1. + gdat.meanreds[c:]), gdat.meanreds[c:])
        
    thertempantntemp = (1.76e-11 * gdat.freqmodl)**2 * exp(1.76e-11 * gdat.freqmodl) / (exp(1.76e-11 * gdat.freqmodl) - 1.)**2
    antntempflux = 0.0307 * (gdat.freqmodl / 1e9)**2

    # scaled frequency axis
    sfrqconc = gdat.plnkcons * gdat.freqmodl / gdat.boltcons / gdat.tempcmbrconc
    
    # cosmological constants
    edencrit = 4e9 # [eV/m^3]
    hubb = 0.68
    omegbmat = 0.049
    omegdmat = 0.26
    omegmatt = omegbmat + omegdmat
    omegradi = 4.8e-5
    omegdene = 0.69

    # cosmological setup
    edenbmat = omegbmat * edencrit * (1. + gdat.meanreds)**3
    edendmat = omegdmat * edencrit * (1. + gdat.meanreds)**3
    edenmatt = omegmatt * edencrit * (1. + gdat.meanreds)**3
    gdat.edenradi = omegradi * edencrit * (1. + gdat.meanreds)**4
    edendene = omegdene * edencrit
    gdat.ndenbmat = edenbmat / massprot
    
    # distortion visibility function
    gdat.vifm = 1. - exp(-((1. + gdat.meanreds) / 5.8e4)**1.88)
    gdat.vify = 1. / (1. + ((1. + gdat.meanreds) / 6e4)**2.58)
    gdat.vift = exp(-(gdat.meanreds / redsdism)**2.5)
    
    fluxcmbrconc = retr_fluxcmbr(gdat, gdat.freqmodl, gdat.tempcmbrconc)
    fluxtdisgrenconc, fluxydisgrenconc, fluxmdisgrenconc, fluxfdisgrenconc, fluxodisgrenconc = retr_fluxgren(gdat, gdat.freqmodl, gdat.tempcmbrconc, disttype='full')
        
    #if makeplot:
        #plot_pure_dist()
        #plot_grnf()
        #plot_cros_plnk()

    numbbins = 20
    
    mocksampvarb = retr_mocksampvarb()
    datapara = retr_datapara()
    numbpara = len(datapara.name)
    indxpara = arange(numbpara)
    
    thissampvarb = copy(mocksampvarb)
    thissamp = tdpy.mcmc.cdfn_samp(thissampvarb, datapara)

    if gdat.verbtype > 1:
        print 'thissampvarb'
        print thissampvarb
        print 'thissamp'
        print thissamp

    # get expected instrumental uncertainty for PIXIE
    path = os.environ["CMBR_DIST_DATA_PATH"] + '/pixifluxstdv.csv'
    gdat.exprfluxstdvinstfreq = loadtxt(path)
    gdat.exprfluxstdvinstfreq = interp1d(gdat.exprfluxstdvinstfreq[:, 0] * 1e9, gdat.exprfluxstdvinstfreq[:, 1] * 1e26)(freqexpr)
    
    if gdat.datatype == 'mock':
         
        mockfluxtotl, mockfluxcmbr, mockfluxdustcold, mockfluxdustwarm, mockfluxsync, mockfluxfree, mockfluxydis, mockfluxdeca = retr_flux(gdat, gdat.freqmodl, mocksampvarb)
        mockfluxintp = interp1d(gdat.freqmodl, mockfluxtotl)(freqexpr)
                
        if gdat.exprtype == 'pixi':
            gdat.datafluxstdv = mockfluxintp * gdat.exprfluxstdvfrac + gdat.exprfluxstdvinstfreq * gdat.exprfluxstdvinst
        if gdat.exprtype == 'plnk':
            gdat.datafluxstdv = mockfluxintp * gdat.exprfluxstdvfrac + gdat.exprfluxstdvinst
            
        gdat.dataflux = mockfluxintp + gdat.datafluxstdv * randn(gdat.datafluxstdv.size)
        
        # temp
        #fluxegbl = retr_egbl(freqexpr)
        #gdat.dataflux += fluxegbl
        
        gdat.resiflux = gdat.dataflux - mockfluxintp
        gdat.path = pathplot + 'mockresi.png'
        # temp
        plot_dataflux(gdat)
    
    else:
        
        gdat.datafluxstdv = exprfluxstdv
        gdat.dataflux = exprflux
        
    #plot_llik()
    
    gdat.numbfreqexpr = freqexpr.size

    # sampler setup
    numbsamp = tdpy.mcmc.retr_numbsamp(numbswep, numbburn, factthin)
    nproc = 1

    swepcntr = 0
    
    listflux = zeros((numbsamp, gdat.numbfreqexpr))
    
    listsampunitfull = []
    
    numbproc = 1
    sampbund = tdpy.mcmc.init(numbproc, numbswep, retr_llik, datapara, numbburn=numbburn, gdatextr=gdat, factpropeffi=3., \
                                    factthin=factthin, optiprop=optiprop, verbtype=gdat.verbtype, pathbase=pathbase, rtag=rtag)
    listsampvarb = sampbund[0]
    listsamp = sampbund[1]
    listsampcalc = sampbund[2]
    listllik = sampbund[3]
    listaccp = sampbund[4]
    listjsampvari = sampbund[5]

    statpara = zeros((numbpara, 3))
    statpara[:, 0] = percentile(listsampvarb, 10., axis=0)
    statpara[:, 1] = percentile(listsampvarb, 50., axis=0)
    statpara[:, 2] = percentile(listsampvarb, 90., axis=0)

    # smooth the posterior using the KDE
    #meshtimedeca, meshampl = meshgrid(binstimedeca, binsampldeca)
    #grid = vstack([meshtimedeca.ravel(), meshampl.ravel()])
    #kernel = sp.stats.gaussian_kde(samp.T, 0.1)
    #sampsmth = reshape(kernel(grid).T, meshtimedeca.shape)

    if makeplot:

        indxparaself = where(scalpara == 'self')[0]
        indxparalogt = where(scalpara == 'logt')[0]

        listsampvarbtran = empty_like(listsampvarb)
        strgparatran = empty(numbpara, dtype=object)
        mocksampvarbtran = zeros_like(mocksampvarb)
        scalparatran = empty(numbpara, dtype=object)
        scalparatran[:] = 'self'

        listsampvarbtran[:, indxparaself] = listsampvarb[:, indxparaself] / mocksampvarb[None, indxparaself] - 1.
        listsampvarbtran[:, indxparalogt] = log10(listsampvarb[:, indxparalogt] / mocksampvarb[None, indxparalogt])
        listsampvarbtran[:, :numbpara-3] *= 1e6

        strgparatran[:] = ''
        strgparatran[:numbpara-3] += r'$10^6 \times$ '
        strgparatran[indxparaself] += '(' + lablpara[indxparaself] + ' / ' + lablpara[indxparaself] + r'$^{mock}$ - 1)'
        strgparatran[indxparalogt] += 'log(' + lablpara[indxparalogt] + ' / ' + lablpara[indxparalogt] + r'$^{mock}$)' 
        
        if savepost:
            listfluxtotl = empty((numbsamp, numbfreqmodl))
            listfluxcmbr = empty((numbsamp, numbfreqmodl))
            listfluxdustcold = empty((numbsamp, numbfreqmodl))
            listfluxdustwarm = empty((numbsamp, numbfreqmodl))
            listfluxsync = empty((numbsamp, numbfreqmodl))
            listfluxfree = empty((numbsamp, numbfreqmodl))
            listfluxydis = empty((numbsamp, numbfreqmodl))
            listfluxdeca = empty((numbsamp, numbfreqmodl))
            listresiflux = empty((numbsamp, gdat.numbfreqexpr))
            for k in range(numbsamp):
                listfluxtotl[k, :] = listsampcalc[0][k]
                listfluxcmbr[k, :] = listsampcalc[1][k]
                listfluxdustcold[k, :] = listsampcalc[2][k]
                listfluxdustwarm[k, :] = listsampcalc[3][k]
                listfluxsync[k, :] = listsampcalc[4][k]
                listfluxfree[k, :] = listsampcalc[5][k]
                listfluxydis[k, :] = listsampcalc[6][k]
                listfluxdeca[k, :] = listsampcalc[7][k]
                listresiflux[k, :] = listsampcalc[8][k]

            path = pathplot + 'postresi.png'
            plot_resiflux(gdat, postflux=postflux) 
            #postfluxtotl, postfluxcmbr, postfluxdustwarm, postfluxdustcold, postfluxsync, postfluxfree, postfluxydis, postfluxdeca)

    return statpara
    

def intr_fluxdust():
    
    numbbins = 4
    
    logtminmdustodep = log10(minmpara[indxpara['dustodep']])
    logtmaxmdustodep = log10(maxmpara[indxpara['dustodep']])
    ndustodep = (logtmaxmdustodep - logtminmdustodep) / numbbins
    
    print 'logtminmdustodep'
    print logtminmdustodep
    print 'logtmaxmdustodep'
    print logtmaxmdustodep
    print 'ndustodep'
    print ndustodep
    
    
    logtminmdustemisrati = log10(minmpara[indxpara['dustemisrati']])
    logtmaxmdustemisrati = log10(maxmpara[indxpara['dustemisrati']])
    ndustemisrati = (logtmaxmdustemisrati - logtminmdustemisrati) / numbbins
    
    minmdustpowrfrac = minmpara[indxpara['dustpowrfrac']]
    maxmdustpowrfrac = maxmpara[indxpara['dustpowrfrac']]
    ndustpowrfrac = (maxmdustpowrfrac - minmdustpowrfrac) / numbbins
    
    minmdustwarmindx = minmpara[indxpara['dustwarmindx']]
    maxmdustwarmindx = maxmpara[indxpara['dustwarmindx']]
    ndustwarmindx = (maxmdustwarmindx - minmdustwarmindx) / numbbins

    minmdustwarmtemp = minmpara[indxpara['dustwarmtemp']]
    maxmdustwarmtemp = maxmpara[indxpara['dustwarmtemp']]
    ndustwarmtemp = (maxmdustwarmtemp - minmdustwarmtemp) / numbbins
    
    minmdustcoldindx = minmpara[indxpara['dustcoldindx']]
    maxmdustcoldindx = maxmpara[indxpara['dustcoldindx']]
    ndustcoldindx = (maxmdustcoldindx - minmdustcoldindx) / numbbins


def plot_fluxdust_wrap(   # logtdustodep, logtdustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, \
  dustcoldindx):
    
    logtdustodep = -4.
    logtdustemisrati  = 1.
    dustpowrfrac = 0.05
    dustwarmindx = 2.5
    dustwarmtemp = 20.
    plot_fluxdust(10**logtdustodep, 10**logtdustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx)
    
    
def plot_fluxdust(gdat, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx):
    
    fluxdust, fluxdustcold, fluxdustwarm = retr_fluxdust(gdat, gdat.freqmodl, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.loglog(gdat.freqmodl * 1e-9, 1e-6 * fluxdust, label='Total')
    ax.loglog(gdat.freqmodl * 1e-9, 1e-6 * fluxdustcold, label='Cold')
    ax.loglog(gdat.freqmodl * 1e-9, 1e-6 * fluxdustwarm, label='Warm')
    ax.set_title('Two-component thermal dust SED')
    ax.set_xlabel(r'$\nu$ [GHz]')
    ax.set_ylabel(r'$I_\nu$ [MJy/sr]')
    ax.legend(loc=2)
    ax.set_ylim([1e-6, 1e5])
    plt.show()


def retr_plnkflux():
    
    plnkflux = loadtxt(os.environ["CMBR_DIST_DATA_PATH"] + '/plnkflux.dat')
    plnkfluxstdv = loadtxt(os.environ["CMBR_DIST_DATA_PATH"] + '/plnkfluxstdv.dat')

    return plnkflux, plnkfluxstdv


def retr_plnkfreq():
    
    freqexpr = array([3e10, 4.4e10, 7e10, 1e11, 1.43e11, 2.17e11, 3.53e11, 5.45e11, 8.57e11]) # [Hz]
    gdat.numbfreqexpr = freqexpr.size
    freqexprstdv = empty(gdat.numbfreqexpr)
    freqexprstdv[:3] = 0.2
    freqexprstdv[3:] = 0.33
    
    gdat.exprfluxstdvinst = 5e3 # [Jy/sr]
    gdat.exprfluxstdvfrac = 1e-3
    
    return freqexpr, freqexprstdv, gdat.exprfluxstdvinst, gdat.exprfluxstdvfrac
    

def writ_plnk():

    nfreqplnk = 9
    
    freqexpr, freqexprstdv, gdat.exprfluxstdvinst, gdat.exprfluxstdvfrac = retr_plnkfreq()

    nside = 256
    npixl = nside**2 * 12
    flux = zeros((npixl, nfreqplnk))
    fluxstdv = zeros((npixl, nfreqplnk))
    mpixl = []

    
    #ipixlring = hp.pixelfunc.nest2ring(nside, arange(npixl))
    
    for k in range(nfreqplnk):
        
        print 'Processing Planck Map at %d GHz...' % (freqexpr[k] / 1e9)
        
        if k >= 3:
            path = os.environ["CMBR_DIST_DATA_PATH"] + '/HFI_SkyMap_%03d_2048_R2.02_nominal.fits' % (freqexpr[k] / 1e9)
        else:
            path = os.environ["CMBR_DIST_DATA_PATH"] + '/LFI_SkyMap_%03d_1024_R2.01_full.fits' % (freqexpr[k] / 1e9)
            
        fluxtemp = pf.getdata(path, 1)

        if k < 7:
            frac = 1e6 * interp1d(gdat.freqmodl, thertempantntemp * antntempflux)(freqexpr[k]) * 1e6
        else:
            frac = 1e6
    
        flux[:, k] = hp.pixelfunc.ud_grade(fluxtemp['I_Stokes'], nside, order_in='NESTED', order_out='RING') * frac
        fluxstdv[:, k] = sqrt(hp.pixelfunc.ud_grade(fluxtemp['II_cov'], nside, order_in='NESTED', order_out='RING')) * frac
        
    path = os.environ["CMBR_DIST_DATA_PATH"] + '/plnkflux.dat'
    savetxt(path, flux)
    
    path = os.environ["CMBR_DIST_DATA_PATH"] + '/plnkfluxstdv.dat'
    savetxt(path, fluxstdv)



def plot_plnktran():
    
    freq, tran = retr_plnktran()

    fig, ax = plt.subplots(figsize=(10, 8))
    for k in range(9):
        ax.loglog(freq[k] * 1e-9, tran[k])
    ax.set_ylim([1e-13, 1e-9])
    ax.set_xlim([1e1, 1e3])
    ax.set_ylabel(r'$T(\nu)$')
    ax.set_xlabel(r'$\nu$ [GHz]')
    plt.savefig(os.environ["CMBR_DIST_DATA_PATH"] + '/png/plnktran.png')
    plt.close()
    
    plt.show()


def retr_plnktran():
    
    freq = []
    tran = []

    path = os.environ["CMBR_DIST_DATA_PATH"] + '/LFI_RIMO_R1.12.fits'
    for i in range(3):
        data, hdr = pf.getdata(path, i+2, header=True)
        freqtemp = 1e9 * data['WAVENUMBER']
        trantemp = data['TRANSMISSION']
        trantemp /= trapz(trantemp, freqtemp) 
        freq.append(freqtemp)
        tran.append(trantemp)

    path = os.environ["CMBR_DIST_DATA_PATH"] + '/HFI_RIMO_R1.10.fits'
    for i in range(6):
        data, hdr = pf.getdata(path, i+2, header=True)
        freqtemp = 1e2 * gdat.velolght * data['WAVENUMBER']
        trantemp = data['TRANSMISSION']
        trantemp /= trapz(trantemp, freqtemp) 
        freq.append(freqtemp)
        tran.append(trantemp)
        
    return freq, tran
    

def retr_rtag(gdat):
    
    rtag = '_%d_%03.1f_%03.1f_%03.1f_%03.1f' % (gdat.numbfreqexpr, log10(gdat.minmfreqexpr), log10(gdat.maxmfreqexpr), \
                log10(gdat.exprfluxstdvinst * 1e3), -log10(gdat.exprfluxstdvfrac))
 
    return rtag


def chek_plnk():
    
    indxpixltemp = random_integers(0, 12*256**2, size=100)
    exprflux = loadtxt(os.environ["CMBR_DIST_DATA_PATH"] + '/plnkflux.dat')[indxpixltemp, :]
    exprfluxstdv = loadtxt(os.environ["CMBR_DIST_DATA_PATH"] + '/plnkfluxstdv.dat')[indxpixltemp, :]

    freqexpr, freqexprstdv, gdat.exprfluxstdvinst, gdat.exprfluxstdvfrac = retr_plnkfreq()

    freqtran, tran = retr_plnktran()
    gdat.numbfreqexpr = freqexpr.size

    #bolo = empty(gdat.numbfreqexpr)
    #for i in range(gdat.numbfreqexpr):
    #    bolo[i] = trapz(tran[i][1:] * exprflux[0, i] * freqexpr[i] / freqtran[i][1:], freqtran[i][1:])
    #    plt.loglog(freqtran[i] * 1e-9, tran[i])  
    #plt.show()
    #plt.loglog(freqexpr * 1e-9, exprfluxstdv)
    #plt.show()

    for k in range(10):
        plt.loglog(freqexpr * 1e-9, exprflux[k, :])
    plt.show()

    plt.loglog(freqexpr * 1e-9, bolo)
    plt.show()
    plt.loglog(freqexpr * 1e-9, 100. * (bolo - exprflux) / exprflux)
    plt.show()


def plot_cros_plnk():

    plnkflux, plnkfluxstdv = retr_plnkflux()
    
    nfreqplnk = plnkflux.shape[1]
    
    strgpara = []
    for k in range(nfreqplnk):
        strgpara.append(r'$\mathcal{I}_{%d}$' % (freqexpr[k] / 1e9))
    scalpara = ['self'] * nfreqplnk
    path = os.environ["CMBR_DIST_DATA_PATH"] + '/png/plnkcros'
    tdpy.mcmc.plot_grid(plnkflux * 1e-6, strgpara, path=path, scalpara=scalpara)
    

def cnfg_arca():

    # ARCADE frequency axis
    cnfg['freqarca'] = array([3., 5., 8., 10., 30., 90.]) * 1e9
    cnfg['fluxarca'] = 24.1 * (freqarca / 0.31e9)**(-2.6)


def cnfg_plnk_mock():

    freqexpr, freqexprstdv, gdat.exprfluxstdvinst, gdat.exprfluxstdvfrac = retr_plnkfreq()
    gdat.numbfreqexpr = freqexpr.size
    
    statpara = init( \
                    freqexpr=freqexpr, \
                    inclcmbrmono=False, \
                    freqexprstdv=freqexprstdv, \
                    exprfluxstdvinst=gdat.exprfluxstdvinst, \
                    exprfluxstdvfrac=gdat.exprfluxstdvfrac \
                   )
    

def cnfg_plnk_expr():
    
    exprflux, exprfluxstdv = retr_plnkflux()
    
    # temp
    exprflux = exprflux[0, :]
    exprfluxstdv = exprfluxstdv[0, :]
    
    freqexpr, freqexprstdv, gdat.exprfluxstdvinst, gdat.exprfluxstdvfrac = retr_plnkfreq()
    gdat.numbfreqexpr = freqexpr.size

    cnfg = init( \
                numbswep=100, \
                datatype='inpt', \
                inclcmbrmono=False, \
                freqexpr=freqexpr, \
                freqexprstdv=freqexprstdv, \
                exprfluxstdvinst=gdat.exprfluxstdvinst, \
                exprfluxstdvfrac=gdat.exprfluxstdvfrac, \
                exprflux=exprflux, \
                exprfluxstdv=exprfluxstdv \
               )
    

def cnfg_pixi_mock():
    
    minmfreqexpr = 3e10 # [Hz]
    maxmfreqexpr = 6e12 # [Hz]
    numbfreqexpr = 400
    
    freqexprstdv = ones(numbfreqexpr) * 0.01
    
    exprfluxstdvinst = 5e0 # [Jy/sr]
    exprfluxstdvfrac = 1e-20
    
    cnfg = init( \
                numbswep=5000, \
                verbtype=1, \
                inclcmbrmono=True, \
                numbfreqexpr=numbfreqexpr, \
                minmfreqexpr=minmfreqexpr, \
                maxmfreqexpr=maxmfreqexpr, \
                freqexprstdv=freqexprstdv, \
                exprfluxstdvinst=exprfluxstdvinst, \
                exprfluxstdvfrac=exprfluxstdvfrac \
               )
    

def cnfg_pixi_mock_stdv():
    
    numbswep = 10000
        
    datapara = retr_datapara()
    numbpara = len(lablpara)
    
    strgpertpara = [r'$\nu_{min}$', r'$\nu_{max}$', r'$N_\nu$', r'$\sigma$', r'$\sigma_f$']
    numbpertpara = 4
    numbpert = 1
    arryminmfreqexpr = logspace(log10(3e9), log10(3e11), numbpert)
    arrymaxmfreqexpr = logspace(log10(1.5e12), log10(6e12), numbpert)
    arrynumbfreqexpr = linspace(100, 700, numbpert)
    arryexprfluxstdvinst = logspace(log10(5e0), log10(5e4), numbpert)
    
    arryexprfluxstdvfrac = logspace(-10., -6., numbpert)
    
    statparagrid = zeros((numbpertpara, numbpert, numbpara, 3))
    for k in range(numbpertpara):
        
        minmfreqexpr = 3e10 # [Hz]
        maxmfreqexpr = 3e12 # [Hz]
        numbfreqexpr = 400
        exprfluxstdvinst = 5e0 # [Jy/sr]
        exprfluxstdvfrac = 1e-10
        
        for l in range(numbpert):
            
            if k == 0:
                gdat.minmfreqexpr = arrygdat.minmfreqexpr[l]
            if k == 1:
                gdat.maxmfreqexpr = arrygdat.maxmfreqexpr[l]
            if k == 2:
                gdat.numbfreqexpr = arrygdat.numbfreqexpr[l]
            if k == 3:
                gdat.exprfluxstdvinst = arrygdat.exprfluxstdvinst[l]
            if k == 4:
                gdat.exprfluxstdvfrac = arrygdat.exprfluxstdvfrac[l]

            freqexprstdv = ones(gdat.numbfreqexpr) * 0.01
            
            statparagrid[k, l, :, :] = init( \
                                            numbswep=numbswep, \
                                            verbtype=1, \
                                            inclcmbrmono=True, \
                                            numbfreqexpr=numbfreqexpr, \
                                            minmfreqexpr=minmfreqexpr, \
                                            maxmfreqexpr=maxmfreqexpr, \
                                            freqexprstdv=freqexprstdv, \
                                            exprfluxstdvinst=exprfluxstdvinst, \
                                            exprfluxstdvfrac=exprfluxstdvfrac \
                                           )

    for k in range(numbpertpara):
        
        path = os.environ["CMBR_DIST_DATA_PATH"] + '/png/stdv%d.png' % k
        fig, axgr = plt.subplots(numbpara / 2, 2, figsize=(14, numbpara * 3))

        if k == 0:
            xaxi = arrygdat.minmfreqexpr
        if k == 1:
            xaxi = arrygdat.maxmfreqexpr
        if k == 2:
            xaxi = arrygdat.numbfreqexpr
        if k == 3:
            xaxi = arrygdat.exprfluxstdvinst
        if k == 4:
            xaxi = arrygdat.exprfluxstdvfrac
            
        for a, axrw in enumerate(axgr):
            for b, ax in enumerate(axrw):

                m = 2 * a + b
                
                ax.plot(xaxi, statparagrid[k, :, m, 0])
                ax.plot(xaxi, statparagrid[k, :, m, 1])
                ax.plot(xaxi, statparagrid[k, :, m, 2])
                ax.set_xlabel(strgpertpara[k])
                ax.set_ylabel(lablpara[m] + ' ' + unitpara[m])
                if scalpara[m] == 'logt':
                    ax.set_yscale('log')
                if k != 2:
                    ax.set_xscale('log')
                    
                ax.set_xlim([amin(xaxi), amax(xaxi)])
                ax.set_ylim([amin(statparagrid[k, :, m, :]), amax(statparagrid[k, :, m, :])])

        plt.savefig(path)
        plt.close(fig)


def plot_plnkmaps():
    
    exprflux, exprfluxstdv = retr_plnkflux()
    nside = exprflux.shape[0]
    npixl = 12 * nside**2
    for k in range(9):
        hp.visufunc.mollview(exprflux[:, k])
        plt.show()
        hmapcart = tdpy.util.retr_cart(exprflux[:, k])
        plt.imshow(hmapcart)


if __name__ == '__main__':   
    
    #cnfg_pixi_mock_stdv('tdpy')
    #cnfg_pixi_mock_stdv('emce')
    cnfg_pixi_mock()
    #cnfg_plnk_expr()
    #intr_fluxdust()
    #writ_plnk()
    
    if os.uname()[1] == 'fink1.rc.fas.harvard.edu':
        cmnd = 'cp -r ' + os.environ["CMBR_DIST_DATA_PATH"] + '/png/* /n/pan/www/tansu/png/cmbr_dist/'
        os.system(cmnd)

