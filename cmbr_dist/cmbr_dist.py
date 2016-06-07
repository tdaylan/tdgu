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
    ax[0].plot(reds, vify, label=r'$\mathcal{J}_y(z_i) = \left(1 + \left(\frac{1+z_i}{6\times10^4}\right)\right)^{-1}$')
    ax[0].plot(reds, vift, label=r'$\mathcal{J_B}(z_i) = \exp\left(-\frac{z_i}{2\times10^6}\right)^{2.5}$')
    ax[0].plot(reds, vifm, label=r'$\mathcal{J}_\mu(z_i) = 1 - \exp\left(-\frac{1+z_i}{6\times10^4}\right)^{1.88}$')
    ax[0].plot(reds, vifm * vift, label=r'$\mathcal{J}_B(z_i)\mathcal{J}_\mu(z_i)$')
    ax[0].set_xscale('log')
    ax[1].set_xlabel('$z$')
    ax[0].set_ylabel(r'$\mathcal{J}(z)$')
    ax[0].legend(fontsize=16, loc='center left', bbox_to_anchor=[0.,0.5])
    ax[1].plot(reds, vifm * vift + vify - vift, 'y')
    ax[1].set_xscale('log')
    ax[1].set_ylabel('$\mathcal{J}_B\mathcal{J}_\mu + \mathcal{J}_y - \mathcal{J}_B$')
    plt.savefig(pathplot + 'visifunc.png')
    plt.close()

    with sns.color_palette("Blues", njreds):
        fig, ax = plt.subplots(2, 1, sharex='all')
        for c in range(njreds):
            ax[0].plot(freqmodl * 1e-9, 1e-6 * fluxfdisgrenconc[:,jreds[c]], label='$z$ = %3.1g' % reds[jreds[c]])
        ax[0].set_xscale('log')
        ax[0].set_title('Full distortion')
        for c in range(njreds):
            ax[1].plot(freqmodl * 1e-9, 1e-6 * fluxodisgrenconc[:,jreds[c]], label='$z$ = %3.1g' % reds[jreds[c]])
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
    
    timedcom = 1e13 * (reds / 9e3)**(-5.) 
    timebrem = 1e2 * (reds / 1e7)**(-2.5)

    timecomp = 1e-10 * (reds / 6e5)**(-4.)
    timecoul = 1e-3 * (reds / 1e3)**(-1.5)

    timebose = 4e-7 * (reds / 1e7)**(-4.)

    timether0 = 5e12 * (reds / 1e3)**(-4.)
    timether1 = 5e-3 * (reds / 1e7)**(-5)
    timether = minimum(timether0, timether1)

    fig, ax = plt.subplots()
    fig.suptitle('Time scales relevant to Comptonization in the early Universe', fontsize=20)

    #ax.loglog(reds, timebrem, label='Bremsstrahlung')
    #ax.loglog(reds, timedcom, label='Double Compton')
    #ax.loglog(reds, timecomp, label='Compton')
    ax.loglog(reds, timecoul, label='Coulomb')
    ax.loglog(reds, timebose, label='Compton')
    ax.loglog(reds, timether, label='Thermalization')
    ax.loglog(reds, timehubb / yr2s, color='black', label='Hubble')

    line0 = ax.axvline(5e4, ls='--', color='grey')
    line1 = ax.axvline(2e6, ls='-.', color='grey')

    leg0 = plt.legend([line0, line1], ['Compton Surface', 'Blackbody Surface'], loc='center', bbox_to_anchor=[0.2, 0.1])
    leg1 = ax.legend(loc='center', bbox_to_anchor=[0.65, 0.83])
    ax.add_artist(leg0)

    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$\tau$ [yr]')

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amax(time[0:-1]), amin(time[0:-1])])
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
    ax.plot(freqmodl * 1e-9, fluxcmbrconc * 1e-6, color='black', label=r'$I^B_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{1}{e^x-1}$')
    ax.plot(freqmodl * 1e-9, fluxtdisgrenconc * 1e-6, label=r'$I^T_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{xe^x}{(e^x-1)^2}$')
    ax.plot(freqmodl * 1e-9, fluxydisgrenconc * 1e-6, label=r'$I^Y_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{xe^x}{(e^x-1)^2}(x\coth (x/2) - 4)$')
    ax.plot(freqmodl * 1e-9, fluxmdisgrenconc * 1e-6, label=r'$I^M_\nu(\nu) = \frac{2h\nu^3}{c^2}\frac{e^x}{(e^x-1)^2}\left(x/2.2-1\right)$')
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
    
    heatstrg = ['Silk damping',                 r's-wave annihilation, $<\sigma v> = 3 \times 10^{-26}$ cm$^3$/s',                 r'p-wave annihilation (S-enhanced), $<\sigma v> = 3 \times 10^{-31}$ cm$^3$/s',                 r'p-wave annihilation, $<\sigma v> = 3 \times 10^{-36}$ cm$^3$/s',                 r'Decay, $\tau = 1$ year']
    heatrtag = ['silk', 'swav', 'pwrl', 'pwnr', 'deca']

    csecdmatswav = 3e-20
    csecdmatpwrl = 3e-25
    csecdmatpwnr = 3e-30

    massdmat = 1e11 # [eV]
    timedeca = 1. # [yr]
    ratedeca = 1. / timedeca
    redsdeca = interp1d(time[::-1], reds[::-1])(timedeca)
    ndendmat = edendmat / massdmat

    ninjetype = 5
    heat = zeros((nreds, ninjetype))
    heat[:, 0] = 0.1 * edendmat / massdmat * ratedeca
    heat[:, 1] = ndendmat**2 * csecdmatswav
    heat[:, 2] = ndendmat**2 * csecdmatpwrl * (1. + reds)
    heat[:, 3] = ndendmat**2 * csecdmatpwnr * (1. + reds)**2
    heat[:, 4] = 0.05 * ndendmat * ratedeca * exp(-(redsdeca / reds)**2)

    for k in range(ninjetype):
        heat[:, k] *= timehubb / (1. + reds) / edenradi

    fig, ax = plt.subplots()
    for k in range(1, ninjetype):
        ax.loglog(reds, reds * heat[:, k], label=heatstrg[k])
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$d\kappa/d\ln z$')
    ax.set_ylim([1e-12, 1e-6])
    ax.legend(loc=2)

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amax(time[0:-1]), amin(time[0:-1])])
    ax_.set_xlabel('$t$ [yr]')

    plt.savefig(pathplot + 'heatrate.png')
    plt.close()
    
    edendmatextd = omegdmat * edencrit * (1. + redsextd)**3
    edenradiextd = omegradi * edencrit * (1. + redsextd)**4
    ndendmatextd = edendmatextd / massdmat
    timehubbextd = 4.5e17 * (omegmatt * (1. + redsextd)**3 + omegradi * (1. + redsextd)**4.)**(-0.5)
    timeextd = zeros_like(redsextd)
    for c in range(nreds-1):
        timeextd[c] = trapz(timehubbextd[c:] / (1. + redsextd[c:]), redsextd[c:])

    heatextd = zeros((nreds, ninjetype))
    heatextd[:, 1] = ndendmatextd**2 * csecdmatswav
    heatextd[:, 2] = ndendmatextd**2 * csecdmatpwrl * (1. + redsextd)
    heatextd[:, 3] = ndendmatextd**2 * csecdmatpwnr * (1. + redsextd)**2
    for k in range(1, 4):
        heatextd[:, k] *= timehubbextd / (1. + redsextd) / edenradiextd

    fig, ax = plt.subplots()
    for k in range(1, 4):
        ax.loglog(redsextd, redsextd * heatextd[:, k], label=heatstrg[k])
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$d\kappa/d\ln z$')
    #ax.set_ylim([1e-12, 1e-1])
    ax.legend(loc=2)

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amax(timeextd[0:-1]), amin(timeextd[0:-1])])
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
    
    dist = zeros((nreds, ndeca))
    heat = zeros((nreds, ndeca))
    redsdeca = zeros(ndeca)
    for k in range(ndeca):
        disttemp, heattemp, redsdecatemp = retr_fluxdeca(fluxodisgrenconc, ampldeca, timedecalist[k])
        heatdeca[:, k] = heattemp
        heatdeca[:, k] = dist
        redsdeca[k] = redsdecatemp


    fig, ax = plt.subplots()
    for k in range(njdeca):
        ax.loglog(reds, reds * heatdeca[:, jdeca[k]], label=r'$\tau = %d$ yr' % timedecalist[jdeca[k]])
        ax.set_xlabel(r'$z$')
        ax.set_ylabel(r'$d\kappa/d\ln z$')
        ax.set_ylim([1e-10, 1e-6])
        ax.legend(loc=2, ncol=4)

        ax_ = ax.twiny()
        ax_.set_xscale('log')
        ax_.set_xlim([amax(time[0:-1]), amin(time[0:-1])])
        ax_.set_xlabel('$t$ [yr]')

    plt.savefig(pathplot + 'heatdeca.png')
    plt.close()
    

    timedecaplot = logspace(-4., 6, 10)
    with sns.color_palette("Blues", njreds): 
        fig, ax = plt.subplots()    
        for k in range(timedecaplot.size):
            ax.loglog(freqmodl * 1e-9, abs(retr_fluxdeca(fluxodisgren, 1., timedecaplot[k])), label='$t = %.3g$' % timedecaplot[k])
        ax.set_xscale('log')
        ax.set_xlabel(r'$\nu$ [GHz]')
        ax.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')
        ax.legend(loc=4)
        plt.close()


def plot_sampdist():
    diffdistdiffreds = heat[None, :, :] * fluxodisgren[:, :, None]
    dist = zeros((nfreq, ninjetype))
    for k in range(ninjetype):
        dist[:, k] = trapz(diffdistdiffreds[:, :, k], reds, axis=1)

    fig, ax = plt.subplots()
    for k in range(1, ninjetype):
        ax.plot(freqmodl * 1e-9, dist[:, k], label=heatstrg[k])
    ax.set_xscale('log')
    ax.set_xlabel(r'$\nu$ [GHz]')
    ax.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')
    ax.legend(loc=2)
    ax.axhline(5., ls='--', color='grey')
    ax.axhline(-5., ls='--', color='grey')
    ax.fill_between(freqmodl * 1e-9, ones_like(freqmodl) * 5., ones_like(freqmodl) * -5., color='grey')
    ax.text(2, 10, r'PIXIE 1-$\sigma$ sensitivity', color='grey', fontsize=20)
    plt.savefig(pathplot + 'totldist.png')
    plt.close()
    
    
    diffdistdifflred = heat[None, :, :] * fluxodisgren[:, :, None] * reds[None, :, None]
    for k in range(1, 2):
        ylim = [amin(diffdistdifflred[:, :, k]), amax(diffdistdifflred[:, :, k])]
        for c in range(njreds):
            fig, ax = plt.subplots()
            ax.plot(freqmodl * 1e-9, diffdistdifflred[:, nreds-1-jreds[c], k])
            text = plt.figtext(0.2, 0.8, '$z_i$ = %.3g' % reds[nreds-jreds[c]-1], fontsize=20)

            ax.set_xscale('log')
            ax.set_xlabel(r'$\nu$ [GHz]')
            ax.set_ylabel(r'$d\Delta I_\nu/ d\ln z$ [Jy/sr]')
            ax.set_ylim(ylim)
            ax.set_title(heatstrg[k])

            plt.legend(loc=2, ncol=2)
            plt.savefig(pathplot + 'diffdistdifflred_' + heatrtag[k] + '_%d.png' % jreds[c])
            plt.close(fig)
        

def retr_fluxdeca(fluxodisgren, ampldeca, timedeca):

    ratedeca = 1. / timedeca
    redsdeca = interp1d(time, reds)(timedeca)
    
    heatdeca = ampldeca * ndenbmat * ratedeca * exp(-timedeca / time) * timehubb / (1. + reds) / edenradi

    difffluxdecadiffreds = heatdeca[None, :] * fluxodisgren
    fluxdeca = trapz(difffluxdecadiffreds, reds, axis=1)

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
    imag = ax.imshow(lliktopo, origin='lower', interpolation='none', cmap='Reds',               extent=[minmtimedeca, maxmtimedeca, minmampldeca, maxmampldeca])
    
    sp.stats.chi2.ppf(1 - 0.05, 1) / 2.
    
    levl = zeros(2)
    levl[0] = array([amax(lliktopo) - sp.stats.chi2.ppf(1 - 0.37, 1) / 2.])
    levl[1] = array([amax(lliktopo) - sp.stats.chi2.ppf(1 - 0.05, 1) / 2.])
    cont = ax.contour(binstimedeca, binsampldeca, lliktopo, origin='lower', color='b', levels=levl)
    
    ax.set_ylabel("$f_X$ [eV]")
    ax.set_xlabel(r"$\tau_X$ [year]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    #plt.colorbar(imag, ax=ax)
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


def retr_fluxcmbr(thisfreq, tempcmbr):
    
    sfrq = plnkcons * thisfreq / boltcons / tempcmbr
    
    occpplnk = retr_occp(sfrq)

    fluxcmbr = 2. * plnkcons * thisfreq**3 / velolght**2 * occpplnk * 1e26
    
    return fluxcmbr


def retr_fluxgren(thisfreq, tempcmbr, disttype):
    
    sfrq = plnkcons * thisfreq / boltcons / tempcmbr

    occpplnk, occptdis, occpydis, occpmdis = retr_occp(sfrq, dist=True)
    
    fluxtdisgren = 0.25 * 2. * plnkcons * thisfreq**3 / velolght**2 * occptdis * 1e26
    fluxydisgren = 0.25 * 2. * plnkcons * thisfreq**3 / velolght**2 * occpydis * 1e26
    fluxmdisgren = 1.41 * 2. * plnkcons * thisfreq**3 / velolght**2 * occpmdis * 1e26
    
    fluxfdisgren = vifm[None, :] * vift[None, :] * fluxmdisgren[:, None] +         vify[None, :] * fluxydisgren[:, None] + (1. - vift[None, :]) * fluxtdisgren[:, None]
    fluxodisgren = vifm[None, :] * vift[None, :] * fluxmdisgren[:, None] + vify[None, :] * fluxydisgren[:, None]
    
    if disttype == 'full':
        return fluxtdisgren, fluxydisgren, fluxmdisgren, fluxfdisgren, fluxodisgren
    elif disttype == 'ydisodis':
        return fluxydisgren, fluxodisgren


def retr_fluxsync(thisfreq, syncnorm, syncindx):
    
    fluxsync = syncnorm * (thisfreq / 1e9)**syncindx
    
    return fluxsync
  
    
def retr_fluxfree(thisfreq, emmefree, tempfree):
    
    thisplnkfunc = retr_plnkfunc(thisfreq, tempfree) * 1e26
    
    #print 'thisplnkfunc'
    #print thisplnkfunc
    #print 'emmefree / thisfreq**2.1 / tempfree**1.5'
    #print emmefree / thisfreq**2.1 / tempfree**1.5
    
    gaunfact = log(4.955e2 * (thisfreq / 1e9)**(-1.)) + 1.5 * log(tempfree)
    odepfree = 3.014e-2 * (tempfree)**(-1.5) * (thisfreq / 1e9)**(-2.) * emmefree * gaunfact
    fluxfree = (1. - exp(-odepfree)) * thisplnkfunc
    
    return fluxfree


def retr_plnkfunc(thisfreq, temp):
    
    thissfrq = plnkcons * thisfreq / boltcons / temp
    
    plnk = 2. * plnkcons * thisfreq**3 / velolght**2 / (exp(thissfrq) - 1.)
    
    return plnk


def plot_temp():
    
    #def difftempdmatdiffreds(tempdmat, thisreds):
    #    thishubb = interp1d(reds, hubb)(thisreds)
    #    difftempdmatdiffreds = -2. * thishubb * tempdmat# + gamm * b * (tempbmat - tempdmat)
    #    return difftempdmatdiffreds
    #inittempdmat = array([3000.])
    #tempdmat = odeint(difftempdmatdiffreds, inittempdmat, reds)

    tempmatt = tempcmbrconc * (1. + reds) / (1. + 119.  / (1. + reds) / (1. + ((1. + reds) / 115.)**1.5))
    tempdmatcold = tempcmbrconc * (1. + reds) / (1. + 1e9  / (1. + reds) / (1. + ((1. + reds) / 1e9)**2.5))
    tempdmatwarm = tempcmbrconc * (1. + reds) / (1. + 5e6  / (1. + reds) / (1. + ((1. + reds) / 5e6)**2.5))
    tempcmbr = tempcmbrconc * (1. + reds)

    fig, ax = plt.subplots()
    ax.loglog(reds, tempcmbr, label='CMB')
    ax.loglog(reds, tempdmatwarm, label='WDM')
    ax.loglog(reds, tempdmatcold, label='CDM')
    ax.loglog(reds, tempmatt, label='Baryonic matter')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'T(z) [K]')

    enerradilate = tempcmbrconc * (1. + reds) * boltconsnatu

    ax_ = ax.twiny()
    ax_.set_xscale('log')
    ax_.set_xlim([amin(enerradilate), amax(enerradilate)])
    ax_.set_xlabel(r'$E_\gamma$ [eV]')

    ax.legend(loc=2)
    plt.savefig(pathplot + 'tempevol.png')
    plt.close()
    

def plot_silkscal():
    
    wnumsilk = (5.92e10)**(-0.5) * (1. + reds)**1.5
    wlensilk = 2. * pi / wnumsilk

    wlenahor = (1. + reds) * time * yr2s / mp2m * velolght / sqrt(3. * (1. + edenbmat / edenradi))
    wnumahor = 2. * pi / wlenahor

    wlenphor = (1. + reds) * time * yr2s / mp2m * velolght
    wnumphor = 2. * pi / wlenphor

    fig, ax = plt.subplots()
    ax.loglog(reds, wnumsilk, label='Dissipation scale')
    ax.set_ylabel('$k$ [Mpc$^{-1}$]')
    ax.set_xlabel('$z$')
    ax.axvline(5e4, ls='--', color='black')
    ax.axvline(2e6, ls='--', color='black')
    ax.axhline(interp1d(reds, wnumsilk)(5e4), ls='--', color='black')
    ax.axhline(interp1d(reds, wnumsilk)(2e6), ls='--', color='black')

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


def plot_resi_post(freqexpr, dataflux, listresiflux,              freqmodl, listfluxtotl, listfluxcmbr, listfluxdustwarm, listfluxdustcold,              listfluxsync, listfluxfree, listfluxydis, listfluxdeca, path):

    fig = plt.figure(figsize=(16, 12))
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 3]) 
    ax = [] 
    ax.append(plt.subplot(gs[0]))
    ax.append(plt.subplot(gs[1], sharex=ax[0]))

    ax[1].set_title('')
    ax[1].errorbar(freqexpr * 1e-9, dataflux * 1e-6, ls='none',                    yerr=datafluxstdv*1e-6, xerr=datafluxstdv*1e-9,                    #label=datalabl, \
                   marker='o', markersize=5, color='k')
    
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxtotl * 1e-6, lcol='salmon',    alpha=0.5, dcol='red', mcol='black')
    if inclcmbrmono:
        tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxcmbr * 1e-6, lcol='lightblue',        alpha=0.5, dcol='blue', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxdustcold * 1e-6, lcol='lightgreen',    alpha=0.5, dcol='green', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxdustwarm * 1e-6, lcol='lightgreen',    alpha=0.5, dcol='green', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxsync * 1e-6, lcol='lightyellow',    alpha=0.5, dcol='yellow', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxfree * 1e-6, lcol='lightcoral',    alpha=0.5, dcol='coral', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxydis * 1e-6, lcol='lightyellow',    alpha=0.5, dcol='yellow', mcol='black')
    tdpy.mcmc.plot_braz(ax[1], freqmodl * 1e-9, listfluxdeca * 1e-6, lcol='lightcyan',    alpha=0.5, dcol='cyan', mcol='black')
    
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


def plot_resi(freqexpr, freqexprstdv, dataflux, datafluxstdv, resiflux, datalabl, freqmodl,
              modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, \
              modlfluxfree, modlfluxydis, modlfluxdeca, path=None):
        
    fig = plt.figure(figsize=(16, 12))
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 3]) 
    ax = [] 
    ax.append(plt.subplot(gs[0]))
    ax.append(plt.subplot(gs[1], sharex=ax[0]))

    ax[1].set_title('')
    ax[1].errorbar(freqexpr * 1e-9, dataflux * 1e-6, ls='none',                    yerr=datafluxstdv*1e-6, xerr=datafluxstdv*1e-9,                    #label=datalabl, \
                   marker='o', markersize=5, color='k')
    
    ax[1].plot(freqmodl * 1e-9, modlfluxtotl * 1e-6, label='Total model', color='r')
    if inclcmbrmono:
        ax[1].plot(freqmodl * 1e-9, modlfluxcmbr * 1e-6, label='CMB', color='b')
    ax[1].plot(freqmodl * 1e-9, modlfluxdustcold * 1e-6, label='Cold Dust', color='g', ls='--')
    ax[1].plot(freqmodl * 1e-9, modlfluxdustwarm * 1e-6, label='Warm Dust', color='g', ls='-.')
    ax[1].plot(freqmodl * 1e-9, modlfluxsync * 1e-6, label='Synchrotron', color='y')
    ax[1].plot(freqmodl * 1e-9, modlfluxfree * 1e-6, label='Brem', color='coral')
    ax[1].plot(freqmodl * 1e-9, abs(modlfluxydis) * 1e-6, label='Reionization', color='m')
    ax[1].plot(freqmodl * 1e-9, abs(modlfluxdeca) * 1e-6, label='Particle decay', color='cyan')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel(r'$\nu$ [GHz]')
    ax[1].set_ylabel(r'$I_\nu$ [MJy/sr]')
    ax[1].set_ylim([1e-7, 1e4])
    ax[1].legend(loc=9, ncol=4)

    # temp
    ax[0].set_xlim([minmfreqmodl * 1e-9, maxmfreqmodl * 1e-9])
    ax[1].set_xlim([minmfreqmodl * 1e-9, maxmfreqmodl * 1e-9])

    ax[0].errorbar(freqexpr * 1e-9, resiflux, yerr=datafluxstdv, xerr=freqexprstdv*1e-9,                marker='o', lw=1, ls='none', markersize=5, color='k')
    ax[0].set_ylabel(r'$I_\nu^{res}$ [Jy/sr]')
    ax[0].axhline(0., ls='--', color='black', alpha=0.1)
    
    if path == None:
        plt.show()
    else:
        plt.savefig(path)
        plt.close(fig)


def retr_ydis_trac():
    
    path = pathplot + 'electron_photonspectrum_results.fits'
    data, hdr = pf.getdata(path, 1, header=True)

    redsoutp = data['OUTPUT_REDSHIFT'].squeeze()
    enerpart = 10**data['LOG10ENERGY'].squeeze()
    redsinpt = data['INPUT_REDSHIFT'].squeeze()
    enerphot = data['PHOTON_ENERGY'].squeeze()
    freqphot = enerphot / plnkcons

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

    specchek = 8. * pi * enerphot[:,:,None]**2 / (velolght * plnkcons)**3 / (exp(enerphot[:,:,None] / tempcmbrnunc / boltcons / (1. + redsinpt[None,None,:])) - 1.) # [1/cm^3/eV]

    diffydisdiffreds = data['YDISTORTION'].squeeze()
    ydis = data['CUMULATIVE_YDISTORTION'].squeeze()
    spec = data['CMB_SPECTRUM'].squeeze()
    nredsout = 63
    nredsinp = 63
    nenerpart = 40
    nfreqphot = 500

    for a in range(njenerpart):
        for f in range(njredsinpt):
            plt.plot(freqphot[:,jenerpart[a]] * 1e-9, spec[:,jenerpart[a],jredsinpt[f]] * enerphot[:,jenerpart[a]], label='$E_{inj} = %.3g, z_{inp} = %.3g$' % (enerpart[jenerpart[a]], redsinpt[jredsinpt[f]]))
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
                plt.plot(freqphot[:,jenerpart[a]] * 1e-9, diffydisdiffreds[:,jredsinpt[f],jenerpart[a],jredsoutp[c]],     label='$z_{out} = %.3g, E_{inj} = %.3g, z_{inp} = %.3g$' %     (redsoutp[jredsoutp[c]], enerpart[jenerpart[a]], redsinpt[jredsinpt[f]]))
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
                plt.plot(freqphot[:,jenerpart[a]] * 1e-9, ydis[:,jredsinpt[f],jenerpart[a],jredsoutp[c]],     label='$z_{out} = %.3g, E_{inj} = %.3g, z_{inp} = %.3g$' %     (redsoutp[jredsoutp[c]], enerpart[jenerpart[a]], redsinpt[jredsinpt[f]]))
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
    plt.colorbar(imag, ax=ax, fraction=0.05)

    flux = hp.synfast(clhploww, nside)
    fluxcart = tdpy.util.retr_cart(flux)
    fig, ax = plt.subplots()
    imag = ax.imshow(fluxcart, origin='lower', cmap='Reds')
    plt.colorbar(imag, ax=ax, fraction=0.05)

    fig, ax = plt.subplots()
    ax.plot(lmodhigh, clhphigh)
    ax.plot(lmodloww, clhploww)               


def retr_flux(thisfreq, thispara):

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

    if verbtype > 1:
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
    fluxcmbr = retr_fluxcmbr(thisfreq, tempcmbr)

    # dust
    fluxdust, fluxdustcold, fluxdustwarm =         retr_fluxdust(thisfreq, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx)
    
    # synchrotron
    fluxsync = retr_fluxsync(thisfreq, syncnorm, syncindx)
    
    # Bremsstrahlung
    fluxfree = retr_fluxfree(thisfreq, emmefree, tempfree)
        
    fluxydisgren, fluxodisgren = retr_fluxgren(thisfreq, tempcmbr, disttype='ydisodis')
    
    # free y-distortion
    fluxydis = ydisampl * fluxydisgren
    
    # Particle decay
    fluxdeca = retr_fluxdeca(fluxodisgren, ampldeca, timedeca)
    
    fluxtotl = fluxdust + fluxsync + fluxfree + fluxydis + fluxdeca
    if inclcmbrmono:
        fluxtotl += fluxcmbr

    return fluxtotl, fluxcmbr, fluxdustcold, fluxdustwarm, fluxsync, fluxfree, fluxydis, fluxdeca


def retr_fluxdust(thisfreq, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx):

    factwarm = (sp.special.zetac(4. + dustwarmindx) + 1) * sp.special.gamma(4. + dustwarmindx)
    factcold = (sp.special.zetac(4. + dustcoldindx) + 1) * sp.special.gamma(4. + dustcoldindx)
    dustcoldtemp = (factwarm / factcold / dustemisrati *                     (plnkcons * freqpivt / boltcons)**(dustcoldindx - dustwarmindx) *                     dustwarmtemp**(4. + dustwarmindx))**(1. / (4. + dustcoldindx))
    
    fluxdustcoldfact = dustpowrfrac * dustemisrati * (freqmodl / 3e12)**dustcoldindx
    fluxdustwarmfact = (1. - dustpowrfrac) * (freqmodl / 3e12)**dustwarmindx
        
    fluxdustcold = dustodep * fluxdustcoldfact * retr_plnkfunc(thisfreq, dustcoldtemp) * 1e26       / (fluxdustcoldfact + fluxdustwarmfact)
        
    fluxdustwarm = dustodep * fluxdustwarmfact * retr_plnkfunc(thisfreq, dustwarmtemp) * 1e26       / (fluxdustcoldfact + fluxdustwarmfact)
        
    fluxdust = fluxdustcold + fluxdustwarm
    
    return fluxdust, fluxdustcold, fluxdustwarm


def retr_llik(sampvarb, init=False):
    
    if samptype == 'emce':
        global swepcntr, thisswepcntr
        nextswepcntr = int(20. * swepcntr / (numbswep * numbwalk)) * 5
        if nextswepcntr > thisswepcntr:
            print '%3d%% completed.' % nextswepcntr
            thisswepcntr = nextswepcntr
        swepcntr += 1

    modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, modlfluxfree, modlfluxydis, modlfluxdeca = retr_flux(freqmodl, sampvarb)
       
    modlfluxintp = interp1d(freqmodl, modlfluxtotl)(freqexpr)
    resiflux = dataflux - modlfluxintp

    llik = sum(-log(sqrt(2. * pi) * datafluxstdv) - 0.5 * (modlfluxintp - dataflux) / datafluxstdv**2 * (modlfluxintp - dataflux))

    if savepost:
        sampcalc = modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, modlfluxfree, modlfluxydis, modlfluxdeca, resiflux
    else:
        sampcalc = []

    if samptype == 'emce':
        return llik
    else:
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
    
    dictpara = dict()
    minmpara = zeros(numbpara)
    maxmpara = zeros(numbpara)
    namepara = empty(numbpara, dtype=object)
    scalpara = empty(numbpara, dtype=object)
    lablpara = empty(numbpara, dtype=object)
    unitpara = empty(numbpara, dtype=object)
    varipara = zeros(numbpara)
    
    dictpara['tempcmbr'] = 0
    namepara[0] = 'tempcmbr'
    minmpara[0] = 2.72
    maxmpara[0] = 2.73
    scalpara[0] = 'self'
    lablpara[0] = '$T_{cmb}$'
    unitpara[0] = '[K]'
    varipara[0] = 1e-6

    dictpara['dustodep'] = 1
    namepara[1] = 'dustodep'
    minmpara[1] = 1e-4
    maxmpara[1] = 1e-2
    scalpara[1] = 'logt'
    lablpara[1] = r'$\tau_{dust}$'
    unitpara[1] = ''
    varipara[1] = 1e-8

    dictpara['dustemisrati'] = 2
    namepara[2] = 'dustemisrati'
    minmpara[2] = 1e0
    maxmpara[2] = 1e2
    scalpara[2] = 'logt'
    lablpara[2] = '$q_1/q_2$'
    unitpara[2] = ''
    varipara[2] = 1e-6
    
    dictpara['dustpowrfrac'] = 3
    namepara[3] = 'dustpowrfrac'
    minmpara[3] = 0.
    maxmpara[3] = 0.05
    scalpara[3] = 'self'
    lablpara[3] = '$f_1$'
    unitpara[3] = ''
    varipara[3] = 4e-6
    
    dictpara['dustwarmindx'] = 4
    namepara[4] = 'dustwarmindx'
    minmpara[4] = 2.5
    maxmpara[4] = 3.
    scalpara[4] = 'self'
    lablpara[4] = r'$\beta_2$'
    unitpara[4] = ''
    varipara[4] = 1e-6
    
    dictpara['dustwarmtemp'] = 5
    namepara[5] = 'dustwarmtemp'
    minmpara[5] = 10.
    maxmpara[5] = 20.
    scalpara[5] = 'self'
    lablpara[5] = '$T_2$'
    unitpara[5] = '[K]'
    varipara[5] = 5e-8
    
    dictpara['dustcoldindx'] = 6
    namepara[6] = 'dustcoldindx'
    minmpara[6] = 1.
    maxmpara[6] = 2.5
    scalpara[6] = 'self'
    lablpara[6] = r'$\beta_1$'
    unitpara[6] = ''
    varipara[6] = 8e-6
    
    dictpara['syncnorm'] = 7
    namepara[7] = 'syncnorm'
    minmpara[7] = 1e3
    maxmpara[7] = 1e7
    scalpara[7] = 'logt'
    lablpara[7] = '$A_{sync}$'
    unitpara[7] = ''
    varipara[7] = 2e-6
    
    dictpara['syncindx'] = 8
    namepara[8] = 'syncindx'
    minmpara[8] = -1.5
    maxmpara[8] = 0.
    scalpara[8] = 'self'
    lablpara[8] = r'$\alpha_{sync}$'
    unitpara[8] = ''
    varipara[8] = 8e-6
    
    dictpara['emmefree'] = 9
    namepara[9] = 'emmefree'
    minmpara[9] = 1e0
    maxmpara[9] = 1e4
    scalpara[9] = 'logt'
    lablpara[9] = 'EM'
    unitpara[9] = '[pc/cm$^6$]'
    varipara[9] = 5e-5
    
    dictpara['tempfree'] = 10
    namepara[10] = 'tempfree'
    minmpara[10] = 1e1
    maxmpara[10] = 1e3
    scalpara[10] = 'logt'
    lablpara[10] = r'$T_e$'
    unitpara[10] = '[K]'
    varipara[10] = 1e-3
    
    dictpara['ydisampl'] = 11
    namepara[11] = 'ydisampl'
    minmpara[11] = 1e-9
    maxmpara[11] = 1e-5
    scalpara[11] = 'logt'
    lablpara[11] = '$y_{ri}$'
    unitpara[11] = ''
    varipara[11] = 1e-3
    
    dictpara['ampldeca'] = 12
    namepara[12] = 'ampldeca'
    minmpara[12] = 1e-7
    maxmpara[12] = 1e-3
    scalpara[12] = 'logt'
    lablpara[12] = '$f_X$'
    unitpara[12] = '[eV]'
    varipara[12] = 2e-1
    
    dictpara['timedeca'] = 13
    namepara[13] = 'timedeca'
    minmpara[13] = 1e8
    maxmpara[13] = 1e12
    scalpara[13] = 'logt'
    lablpara[13] = r'$\tau_X$'
    unitpara[13] = '[s]'
    varipara[13] = 2e-1
    
    strgpara = lablpara + ' ' + unitpara
    datapara = namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara
    
    return datapara


def retr_egbl(thisfreq):
    
    path = os.environ['CMBR_DIST_DATA_PATH'] + '/egbl.csv'
    egbldata = loadtxt(path)
    wlenegbl = egbldata[:, 0] * 1e-6 # [m]
    freqegbl = flipud(velolght / wlenegbl) # [Hz]
    
    minmfreqegbl = amin(freqegbl)
    maxmfreqegbl = amax(freqegbl)
    
    fluxegbltemp = flipud(egbldata[:, 1]) * 1e26 / freqegbl # [Jy/sr]
    
    fluxegbl = zeros_like(thisfreq)
    jthisfreq = where((thisfreq < maxmfreqegbl) & (minmfreqegbl < thisfreq))[0]
    fluxegbl[jthisfreq] = interp1d(freqegbl, fluxegbltemp)(thisfreq[jthisfreq])

    return fluxegbl
        

def init(cnfg):
    
    global datatype, exprtype, datalabl, verbtype, makeplot
    global datafluxstdv, dataflux, freqexpr, freqexprstdv
    global numbfreqexpr, inclcmbrmono, numbswep, samptype
    global minmfreqexpr, maxmfreqexpr, numbfreqexpr, exprfluxstdvinst, exprfluxstdvfrac
    
    datatype = cnfg['datatype']
    exprtype = cnfg['exprtype']
    datalabl = cnfg['datalabl']
    
    numbswep = cnfg['numbswep']
    numbburn = cnfg['numbburn']
    factthin = cnfg['factthin']
    samptype = cnfg['samptype']
    
    freqexpr = cnfg['freqexpr']
    numbfreqexpr = cnfg['numbfreqexpr']
    minmfreqexpr = cnfg['minmfreqexpr']
    maxmfreqexpr = cnfg['maxmfreqexpr']

    freqexprstdv = cnfg['freqexprstdv']
    
    exprfluxstdvfrac = cnfg['exprfluxstdvfrac']
    exprfluxstdvinst = cnfg['exprfluxstdvinst']
    
    exprflux = cnfg['exprflux']
    exprfluxstdv = cnfg['exprfluxstdv']
    
    inclcmbrmono = cnfg['inclcmbrmono']
    
    plotperd = cnfg['plotperd']
    verbtype = cnfg['verbtype']
    makeplot = cnfg['makeplot']
    optiprop = cnfg['optiprop']
    
    rtag = retr_rtag()
    
    global pathplot
    pathbase = os.environ["CMBR_DIST_DATA_PATH"]
    pathplot = pathbase + '/png/' + rtag + '/'
    cmnd = 'mkdir -p ' + pathplot
    os.system(cmnd)
    
    if freqexpr == None:
        freqexpr = logspace(log10(minmfreqexpr), log10(maxmfreqexpr), numbfreqexpr) # [Hz]
    
    global freqmodl, minmfreqmodl, maxmfreqmodl
    numbfreqmodl = 1000
    minmfreqmodl = 1e9
    maxmfreqmodl = 1e13
    freqmodl = logspace(log10(minmfreqmodl), log10(maxmfreqmodl), numbfreqmodl) # [Hz]
    
    
    # physical constants
    global velolght, mp2m, yr2s, plnkcons, boltcons, tempcmbrconc, boltconsnatu,         alph, redsdism, massprot, freqpivt
    velolght = 3e8 # [m/s]
    mp2m = 3.1e22 # [Mpc/m]
    yr2s = 364. * 3600. * 24. # [year/s]
    plnkcons = 6.63e-34 # [J s]
    boltcons = 1.38e-23 # [J/K]
    freqpivt = 1e12 # [Hz]
    tempcmbrconc = 2.725 # Planck concordance model temperature of the CMB today [K]
    boltconsnatu = 8.6173e-5 # Boltzmann constant [eV/K]
    massprot = 9.38e8 # [eV]
    
    # distortion related constants
    alph = 1.401
    redsdism = 2e6

    # temp
    global savepost
    savepost = False
    
    # redshift axis
    global reds, jreds, njreds
    nreds = 100
    minmreds = 1e2
    maxmreds = 1e10
    reds = logspace(log10(minmreds), log10(maxmreds), nreds)
    njreds = 10
    jreds = [k * nreds / njreds for k in range(njreds)]
    diffreds = reds[1:] - reds[0:-1]
    redsextd = logspace(2., 10., nreds)
    
    # time
    global time, timehubb
    timehubb = 4.5e17 * (0.27 * (1. + reds)**3 + 9.2e-5 * (1. + reds)**4.)**(-0.5)
    time = zeros_like(reds)
    for c in range(nreds-1):
        time[c] = trapz(timehubb[c:] / (1. + reds[c:]), reds[c:])
        
    global thertempantntemp, antntempflux
    thertempantntemp = (1.76e-11 * freqmodl)**2 * exp(1.76e-11 * freqmodl) / (exp(1.76e-11 * freqmodl) - 1.)**2
    antntempflux = 0.0307 * (freqmodl / 1e9)**2

    # scaled frequency axis
    global sfrqconc
    sfrqconc = plnkcons * freqmodl / boltcons / tempcmbrconc
    
    # cosmological constants
    edencrit = 4e9 # [eV/m^3]
    hubb = 0.68
    omegbmat = 0.049
    omegdmat = 0.26
    omegmatt = omegbmat + omegdmat
    omegradi = 4.8e-5
    omegdene = 0.69

    # cosmological setup
    global ndenbmat, edenradi, edenbmat
    edenbmat = omegbmat * edencrit * (1. + reds)**3
    edendmat = omegdmat * edencrit * (1. + reds)**3
    edenmatt = omegmatt * edencrit * (1. + reds)**3
    edenradi = omegradi * edencrit * (1. + reds)**4
    edendene = omegdene * edencrit
    ndenbmat = edenbmat / massprot
    
    # distortion visibility function
    global vifm, vify, vift
    vifm = 1. - exp(-((1. + reds) / 5.8e4)**1.88)
    vify = 1. / (1. + ((1. + reds) / 6e4)**2.58)
    vift = exp(-(reds / redsdism)**2.5)
    
    
    global fluxcmbrconc, fluxtdisgrenconc, fluxydisgrenconc, fluxmdisgrenconc,         fluxfdisgrenconc, fluxodisgrenconc
    fluxcmbrconc = retr_fluxcmbr(freqmodl, tempcmbrconc)
    fluxtdisgrenconc, fluxydisgrenconc, fluxmdisgrenconc,         fluxfdisgrenconc, fluxodisgrenconc = retr_fluxgren(freqmodl, tempcmbrconc, disttype='full')
        
    #if makeplot:
        #plot_pure_dist()
        #plot_grnf()
        #plot_cros_plnk()

    global numbbins
    numbbins = 20
    
    global minmpara, maxmpara, scalpara, namepara, lablpara, unitpara, varipara, dictpara, indxpara, numbpara
    mocksampvarb = retr_mocksampvarb()
    datapara = retr_datapara()
    namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara = datapara
    numbpara = len(lablpara)
    indxpara = arange(numbpara)
    
    thissampvarb = copy(mocksampvarb)
    thissamp = tdpy.mcmc.cdfn_samp(thissampvarb, datapara)

    if verbtype > 1:
        print 'thissampvarb'
        print thissampvarb
        print 'thissamp'
        print thissamp

    path = os.environ["CMBR_DIST_DATA_PATH"] + '/pixifluxstdv.csv'
    exprfluxstdvinstfreq = loadtxt(path)
    
    exprfluxstdvinstfreq = interp1d(exprfluxstdvinstfreq[:, 0] * 1e9, exprfluxstdvinstfreq[:, 1] * 1e26)(freqexpr)
    exprfluxstdvinstfreq = exprfluxstdvinstfreq / exprfluxstdvinstfreq[0]
    
    if datatype == 'mock':
         
        mockfluxtotl, mockfluxcmbr, mockfluxdustcold, mockfluxdustwarm,             mockfluxsync, mockfluxfree, mockfluxydis, mockfluxdeca = retr_flux(freqmodl, mocksampvarb)
        mockfluxintp = interp1d(freqmodl, mockfluxtotl)(freqexpr)
                
        if exprtype == 'pixi':
            datafluxstdv = mockfluxintp * exprfluxstdvfrac + exprfluxstdvinstfreq * exprfluxstdvinst
        if exprtype == 'plnk':
            datafluxstdv = mockfluxintp * exprfluxstdvfrac + exprfluxstdvinst
            
        dataflux = mockfluxintp + datafluxstdv * randn(datafluxstdv.size)
        
        # temp
        #fluxegbl = retr_egbl(freqexpr)
        #dataflux += fluxegbl
        
        resiflux = dataflux - mockfluxintp
        path = pathplot + 'mockresi.png'
        # temp
        # path = None
        plot_resi(freqexpr, freqexprstdv, dataflux, datafluxstdv, resiflux, 'Mock', freqmodl, \
            mockfluxtotl, mockfluxcmbr, mockfluxdustcold, mockfluxdustwarm, mockfluxsync, mockfluxfree, mockfluxydis, mockfluxdeca, path=path)
    
    else:
        
        datafluxstdv = exprfluxstdv
        dataflux = exprflux
        
    #plot_llik()
    
    numbfreqexpr = freqexpr.size

    # sampler setup
    numbsamp = tdpy.mcmc.retr_numbsamp(numbswep, numbburn, factthin)
    nproc = 1

    global swepcntr
    swepcntr = 0
    
    global listflux
    listflux = zeros((numbsamp, numbfreqexpr))
    
    global listsampunitfull
    listsampunitfull = []
    
    numbproc = 1
    if samptype == 'emce':
        
        global thisswepcntr, numbwalk
        thisswepcntr = -1.
        numbwalk = 100
        initsamp = thissampvarb[None, :] * (1. + 1e-2 * randn(numbwalk * numbpara).reshape((numbwalk, numbpara)))
        sampler = emcee.EnsembleSampler(numbwalk, numbpara, retr_llik)
        
        sampler.run_mcmc(initsamp, numbswep)
        listsampvarb = sampler.flatchain
        
    else:
        sampbund = tdpy.mcmc.init(numbproc, numbswep, retr_llik, datapara, numbburn=numbburn, factpropeffi=3., \
                                        factthin=factthin, optiprop=optiprop, verbtype=verbtype, pathbase=pathbase, rtag=rtag)
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
        listsampvarbtran[:, 0:numbpara-3] *= 1e6

        strgparatran[:] = ''
        strgparatran[0:numbpara-3] += r'$10^6 \times$ '
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
            listresiflux = empty((numbsamp, numbfreqexpr))
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
            plot_resi_post(freqexpr, dataflux, listresiflux, freqmodl, listfluxtotl, listfluxcmbr, listfluxdustwarm, listfluxdustcold, \
                listfluxsync, listfluxfree, listfluxydis, listfluxdeca, path)
    

    return statpara
    

def intr_fluxdust():
    
    numbbins = 4
    
    logtminmdustodep = log10(minmpara[dictpara['dustodep']])
    logtmaxmdustodep = log10(maxmpara[dictpara['dustodep']])
    ndustodep = (logtmaxmdustodep - logtminmdustodep) / numbbins
    
    print 'logtminmdustodep'
    print logtminmdustodep
    print 'logtmaxmdustodep'
    print logtmaxmdustodep
    print 'ndustodep'
    print ndustodep
    
    
    logtminmdustemisrati = log10(minmpara[dictpara['dustemisrati']])
    logtmaxmdustemisrati = log10(maxmpara[dictpara['dustemisrati']])
    ndustemisrati = (logtmaxmdustemisrati - logtminmdustemisrati) / numbbins
    
    minmdustpowrfrac = minmpara[dictpara['dustpowrfrac']]
    maxmdustpowrfrac = maxmpara[dictpara['dustpowrfrac']]
    ndustpowrfrac = (maxmdustpowrfrac - minmdustpowrfrac) / numbbins
    
    minmdustwarmindx = minmpara[dictpara['dustwarmindx']]
    maxmdustwarmindx = maxmpara[dictpara['dustwarmindx']]
    ndustwarmindx = (maxmdustwarmindx - minmdustwarmindx) / numbbins

    minmdustwarmtemp = minmpara[dictpara['dustwarmtemp']]
    maxmdustwarmtemp = maxmpara[dictpara['dustwarmtemp']]
    ndustwarmtemp = (maxmdustwarmtemp - minmdustwarmtemp) / numbbins
    
    minmdustcoldindx = minmpara[dictpara['dustcoldindx']]
    maxmdustcoldindx = maxmpara[dictpara['dustcoldindx']]
    ndustcoldindx = (maxmdustcoldindx - minmdustcoldindx) / numbbins
    
    temp = interact(plot_fluxdust_wrap,                     #logtdustodep=(logtminmdustodep, logtmaxmdustodep, ndustodep), \
                    #logtdustemisrati=(logtminmdustemisrati, logtmaxmdustemisrati, ndustemisrati), \
                    #dustpowrfrac=(minmdustpowrfrac, maxmdustpowrfrac, ndustpowrfrac), \
                    #dustwarmindx=(minmdustwarmindx, maxmdustwarmindx, ndustwarmindx), \
                    #dustwarmtemp=(minmdustwarmtemp, maxmdustwarmtemp, ndustwarmtemp), \
                    dustcoldindx=(minmdustcoldindx, maxmdustcoldindx, ndustcoldindx))


def plot_fluxdust_wrap(   # logtdustodep, logtdustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, \
  dustcoldindx):
    
    logtdustodep = -4.
    logtdustemisrati  = 1.
    dustpowrfrac = 0.05
    dustwarmindx = 2.5
    dustwarmtemp = 20.
    plot_fluxdust(10**logtdustodep, 10**logtdustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx)
    
    
def plot_fluxdust(dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx):
    

    fluxdust, fluxdustcold, fluxdustwarm = retr_fluxdust(freqmodl, dustodep, 
                  dustemisrati, dustpowrfrac, \
                  dustwarmindx, dustwarmtemp, dustcoldindx)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.loglog(freqmodl * 1e-9, 1e-6 * fluxdust, label='Total')
    ax.loglog(freqmodl * 1e-9, 1e-6 * fluxdustcold, label='Cold')
    ax.loglog(freqmodl * 1e-9, 1e-6 * fluxdustwarm, label='Warm')
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
    numbfreqexpr = freqexpr.size
    freqexprstdv = empty(numbfreqexpr)
    freqexprstdv[:3] = 0.2
    freqexprstdv[3:] = 0.33
    
    exprfluxstdvinst = 5e3 # [Jy/sr]
    exprfluxstdvfrac = 1e-3
    
    return freqexpr, freqexprstdv, exprfluxstdvinst, exprfluxstdvfrac
    

def writ_plnk():

    nfreqplnk = 9
    
    freqexpr, freqexprstdv, exprfluxstdvinst, exprfluxstdvfrac = retr_plnkfreq()

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
            frac = 1e6 * interp1d(freqmodl, thertempantntemp * antntempflux)(freqexpr[k]) * 1e6
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
        freqtemp = 1e2 * velolght * data['WAVENUMBER']
        trantemp = data['TRANSMISSION']
        trantemp /= trapz(trantemp, freqtemp) 
        freq.append(freqtemp)
        tran.append(trantemp)
        
    return freq, tran
    

def retr_cnfg(samptype='emce', \
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
              inclcmbrmono=None, \
              plotperd=10000, \
              verbtype=1, \
              optiprop=False, \
              makeplot=True, \
              ):
        
    cnfg = dict()
    
    # data type and label
    cnfg['datatype'] = datatype
    cnfg['exprtype'] = exprtype
    cnfg['datalabl'] = datalabl

    # sampler setup
    cnfg['numbswep'] = numbswep
    cnfg['numbburn'] = numbburn
    cnfg['factthin'] = factthin
    cnfg['samptype'] = samptype

    # frequency axis
    cnfg['freqexpr'] = freqexpr
    cnfg['numbfreqexpr'] = numbfreqexpr
    cnfg['minmfreqexpr'] = minmfreqexpr
    cnfg['maxmfreqexpr'] = maxmfreqexpr
    
    cnfg['freqexprstdv'] = freqexprstdv
    
    cnfg['exprfluxstdvinst'] = exprfluxstdvinst
    cnfg['exprfluxstdvfrac'] = exprfluxstdvfrac

    cnfg['exprflux'] = exprflux
    cnfg['exprfluxstdv'] = exprfluxstdv
    cnfg['inclcmbrmono'] = inclcmbrmono
    
    # plot frequency
    cnfg['plotperd'] = plotperd
    
    # verbosity level
    cnfg['verbtype'] = verbtype
    cnfg['makeplot'] = makeplot
    cnfg['optiprop'] = optiprop
    
    return cnfg


def retr_rtag():
    
    rtag = samptype + '_%d_%03.1f_%03.1f_%03.1f_%03.1f' % (numbfreqexpr, log10(minmfreqexpr), log10(maxmfreqexpr), log10(exprfluxstdvinst * 1e3), -log10(exprfluxstdvfrac))
 
    return rtag


def chek_plnk():
    
    indxpixltemp = random_integers(0, 12*256**2, size=100)
    exprflux = loadtxt(os.environ["CMBR_DIST_DATA_PATH"] + '/plnkflux.dat')[indxpixltemp, :]
    exprfluxstdv = loadtxt(os.environ["CMBR_DIST_DATA_PATH"] + '/plnkfluxstdv.dat')[indxpixltemp, :]

    freqexpr, freqexprstdv, exprfluxstdvinst, exprfluxstdvfrac = retr_plnkfreq()

    freqtran, tran = retr_plnktran()
    numbfreqexpr = freqexpr.size

    #bolo = empty(numbfreqexpr)
    #for i in range(numbfreqexpr):
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

    freqexpr, freqexprstdv, exprfluxstdvinst, exprfluxstdvfrac = retr_plnkfreq()
    numbfreqexpr = freqexpr.size
    
    cnfg = retr_cnfg( \
                     freqexpr=freqexpr, \
                     inclcmbrmono=False, \
                     freqexprstdv=freqexprstdv, \
                     exprfluxstdvinst=exprfluxstdvinst, \
                     exprfluxstdvfrac=exprfluxstdvfrac \
                    )
    statpara = init(cnfg)
    

def cnfg_plnk_expr():
    
    exprflux, exprfluxstdv = retr_plnkflux()
    
    # temp
    exprflux = exprflux[0, :]
    exprfluxstdv = exprfluxstdv[0, :]
    
    freqexpr, freqexprstdv, exprfluxstdvinst, exprfluxstdvfrac = retr_plnkfreq()
    numbfreqexpr = freqexpr.size

    cnfg = retr_cnfg( \
                     numbswep=100, \
                     datatype='inpt', \
                     inclcmbrmono=False, \
                     freqexpr=freqexpr, \
                     freqexprstdv=freqexprstdv, \
                     exprfluxstdvinst=exprfluxstdvinst, \
                     exprfluxstdvfrac=exprfluxstdvfrac, \
                     exprflux=exprflux, \
                     exprfluxstdv=exprfluxstdv \
                    )
    
    statpara = init(cnfg)
    

def cnfg_pixi_mock():
    
    minmfreqexpr = 3e10 # [Hz]
    maxmfreqexpr = 6e12 # [Hz]
    
    numbfreqexpr = 400
    
    freqexprstdv = ones(numbfreqexpr) * 0.01
    
    exprfluxstdvinst = 5e0 # [Jy/sr]
    exprfluxstdvfrac = 1e-20
    
    cnfg = retr_cnfg( \
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
    
    statpara = init(cnfg)
    

def cnfg_pixi_mock_stdv(samptype):
    
    if samptype == 'emce':
        numbswep = 10
    else:
        numbswep = 10000
        
    datapara = retr_datapara()
    namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara = datapara
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
                minmfreqexpr = arryminmfreqexpr[l]
            if k == 1:
                maxmfreqexpr = arrymaxmfreqexpr[l]
            if k == 2:
                numbfreqexpr = arrynumbfreqexpr[l]
            if k == 3:
                exprfluxstdvinst = arryexprfluxstdvinst[l]
            if k == 4:
                exprfluxstdvfrac = arryexprfluxstdvfrac[l]

            freqexprstdv = ones(numbfreqexpr) * 0.01
            
            cnfg = retr_cnfg( \
                             numbswep=numbswep, \
                             samptype=samptype, \
                             verbtype=1, \
                             inclcmbrmono=True, \
                             numbfreqexpr=numbfreqexpr, \
                             minmfreqexpr=minmfreqexpr, \
                             maxmfreqexpr=maxmfreqexpr, \
                             freqexprstdv=freqexprstdv, \
                             exprfluxstdvinst=exprfluxstdvinst, \
                             exprfluxstdvfrac=exprfluxstdvfrac \
                            )

            statparagrid[k, l, :, :] = init(cnfg)

    for k in range(numbpertpara):
        
        path = os.environ["CMBR_DIST_DATA_PATH"] + '/png/stdv%d.png' % k
        fig, axgr = plt.subplots(numbpara / 2, 2, figsize=(14, numbpara * 3))

        if k == 0:
            xaxi = arryminmfreqexpr
        if k == 1:
            xaxi = arrymaxmfreqexpr
        if k == 2:
            xaxi = arrynumbfreqexpr
        if k == 3:
            xaxi = arryexprfluxstdvinst
        if k == 4:
            xaxi = arryexprfluxstdvfrac
            
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

