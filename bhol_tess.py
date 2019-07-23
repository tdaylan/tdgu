import numpy as np
import os
import matplotlib.pyplot as plt
import scipy
from scipy import fftpack
from scipy import signal
from scipy import optimize

class gdatstrt():
    pass


def retr_chi2(para, *args):
    
    gdat = args[0]['gdat']

    peri = para[0]
    amplellp = para[1]
    ampldopp = para[2]
    ampllens = para[3]
    phaslens = para[4]
    stdvlens = para[5]
    
    lcur = retr_lcur(peri, amplellp, ampldopp, ampllens, phaslens, stdvlens, gdat.boolfull)
    chi2 = sum((lcur - lcurmodl)**2 / gdat.varinois) 

    return chi2


def retr_lcur(peri, amplellp, ampldopp, ampllens, phaslens, stdvlens, boolfull):
    
    #print 'peri'
    #print peri
    #print 'amplellp'
    #print amplellp
    #print 'ampldopp'
    #print ampldopp
    #print 'ampllens'
    #print ampllens
    #print 'stdvnois'
    #print stdvnois
    #print

    lcurellp = amplellp * np.sin(2. * np.pi * gdat.time / peri)
    lcurdopp = -ampldopp * np.sin(2. * np.pi * gdat.time / peri)
    lcurlens = np.zeros_like(lcurdopp) 
    for k in range(10):
        lcurlens += ampllens / np.sqrt(2. * np.pi) / stdvlens * np.exp(-0.5 * ((k + phaslens) * peri - gdat.time)**2 / stdvlens**2)
    
    lcur = lcurellp + lcurdopp + lcurlens
    lcur += stdvnois * np.random.randn(gdat.numbtime)
    
    if gdat.boolfull:
        return lcurellp, lcurdopp, lcurlens, lcur
    else:
        return lcur
    

def retr_lcurdata(lcur):

    lcurdata = (1. + 1e-1 * np.random.randn(gdat.numbtime)) * lcur
    lcurdata[gdat.indxtimelink] = 0.
    #lcurdata[gdat.indxtimelink] = np.nan
    
    return lcurdata


def eval_modl_wrap(inpt, *args):

    print 'Call number %d to eval_modl_wrap()' % gdat.indxwrap

        
def plot_lcur(path, lcur, lcurellp=None, lcurdopp=None, lcurlens=None, titl=None):
    
    figr, axis = plt.subplots(figsize=(6, 3))
    axis.plot(gdat.time, lcur, ls='', marker='o', markersize=1, label='Tot', color='black')
    if lcurellp is not None:
        axis.plot(gdat.time, lcurellp, ls='', marker='o', markersize=1, label='EV')
        axis.plot(gdat.time, lcurdopp, ls='', marker='o', markersize=1, label='DP')
        axis.plot(gdat.time, lcurlens, ls='', marker='o', markersize=1, label='SL')
        axis.legend()
    axis.set_xlabel('T [BJD]')
    if titl != None:
        axis.set_title(titl)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def plot_psdn(path, psdn, psdnellp=None, psdndopp=None, psdnlens=None, titl=None):
    
    figr, axis = plt.subplots(figsize=(6, 3))
    axis.plot(perisamp, psdn**2, ls='', marker='o', markersize=1, label='Tot', color='black')
    if psdnellp is not None:
        axis.plot(perisamp, psdnellp**2, ls='', marker='o', markersize=1, label='EV')
        axis.plot(perisamp, psdndopp**2, ls='', marker='o', markersize=1, label='DP')
        axis.plot(perisamp, psdnlens**2, ls='', marker='o', markersize=1, label='SL')
        axis.legend()
    axis.axvline(peri, ls='--', alpha=0.3, color='black')
    axis.set_xlabel('P [day]')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    

def samp_pararand(gdat):
    
    para = np.empty(gdat.numbpara)
    para[0] = np.random.uniform() * 9. + 1.
    para[1] = (1. + 1e-1 * np.random.randn()) * 10.
    para[2] = (1. + 1e-1 * np.random.randn()) * 1.
    para[3] = (1. + 1e-1 * np.random.randn()) * 1.
    para[4] = (1. + 1e-1 * np.random.randn()) * 0.5
    para[5] = (1. + 1e-1 * np.random.randn()) * 0.1
    para[6] = (1. + 1e-1 * np.random.randn()) * 1e-2

    return para


# initialize
np.random.seed(0)
gdat = gdatstrt()
gdat.pathdata = '/Users/tansu/Desktop/bhol_tess/'
os.system('mkdir -p %s' % gdat.pathdata)
datatype = 'mock'

gdat.numbtime = int(27.4 * 24. * 2.)
gdat.time = np.linspace(0., 27.4, gdat.numbtime)

gdat.indxtimelink = np.where(abs(gdat.time - 13.7) < 2.)[0]

if datatype == 'mock':
    # mock data setup 

    peri = 3. # [days]
    stdvnois = 1e-2
    amplellp = 10.
    ampldopp = 1.
    ampllens = 1.
    phaslens = 0.5
    stdvlens = 0.1
    
    gdat.numbpara = 7
    
    argsdict = {'gdat':gdat}
    args = (argsdict)
    paratrue = np.empty(gdat.numbpara)
    paratrue[0] = peri
    paratrue[1] = amplellp
    paratrue[2] = ampldopp
    paratrue[3] = ampllens
    paratrue[4] = phaslens
    paratrue[5] = stdvlens
    paratrue[6] = stdvnois

    gdat.boolfull = True
    lcurellp, lcurdopp, lcurlens, lcur = retr_lcur(peri, amplellp, ampldopp, ampllens, phaslens, stdvlens, True)

    delttime = 1. / 24. / 2.
    fs = 1. / delttime
    freq, psdn = scipy.signal.periodogram(lcur, fs=fs, window=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1)
    freq, psdnellp = scipy.signal.periodogram(lcurellp, fs=fs, window=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1)
    freq, psdndopp = scipy.signal.periodogram(lcurdopp, fs=fs, window=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1)
    freq, psdnlens = scipy.signal.periodogram(lcurlens, fs=fs, window=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1)
    perisamp = 1. / freq
    
    pathplot = gdat.pathdata + 'lcur.pdf'
    if not os.path.exists(pathplot):
        plot_lcur(pathplot, lcur, lcurellp, lcurdopp, lcurlens)
    
    pathplot = gdat.pathdata + 'psdn.pdf'
    if not os.path.exists(pathplot):
        plot_psdn(pathplot, psdn, psdnellp=psdnellp, psdndopp=psdndopp, psdnlens=psdnlens)
    
    parainit = paratrue
	
    #bounds = []
    #rhobeg = np.empty(gdat.numbpara)
    #for k in range(gdat.numbpara):
    #    if k == 0:
    #        limt = [0., 1.]
    #    else:
    #        limt = [0., 1.]
    #    bounds.append(limt)
    #bounds = None

    gdat.boolfull = False
    numbsamp = 10
    indxsamp = np.arange(numbsamp)
    boolsigntrue = np.ones(numbsamp, dtype=bool)
    boolsigntrue[0] = False
    boolsignpred = np.empty_like(boolsigntrue)

if datatype == 'mock':
    
    if boolsigntrue[k]:
        paratrue = samp_pararand(gdat)
        gdat.lcur = retr_lcur(paratrue[0], paratrue[1], paratrue[2], paratrue[3], paratrue[4], paratrue[5], False)
        gdat.lcurdata = retr_lcurdata(gdat.lcur)
    else:
        gdat.lcurdata = np.random.randn(gdat.numbtime)
else:
    path = '/scratch/tmp/orbit-9/cam2/ccd1/'
    
for k in indxsamp:
    
    freq, gdat.psdn = scipy.signal.periodogram(gdat.lcurdata, fs=fs, window=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1)

    if np.amax(gdat.psdn) > 10.:
        boolsignpred[k] = True
    
        #np.correlate(gdat.psdn, psdnmodl)
         
        #argsdict = {'gdat':gdat}
        #args = (argsdict)
        #print 'Calling the minimizer...'
        #objt = optimize.minimize(retr_lcur, parainit, args=args, bounds=bounds, \
        #                            #method='SLSQP', \
        #                            method='COBYLA', \
        #                            #options={'gtol': 1e-05, 'norm': np.inf, 'maxiter':10000, 'eps': 1.4901161193847656e-08, 'disp': True, 'return_all': False}, \
        #                            options={'rhobeg': rhobeg, 'norm': np.inf, 'maxiter':10000, 'eps': 1.4901161193847656e-08, 'disp': True, 'return_all': False}, \
        #                           )
        #chi2 = objt.fun
        #parafitt = objt.x
        
        #chi2 = np.random.chisquare(1)
        #parafitt = np.random.randn()
        #boolbholpred = chi2 < 2.
    
        boolbholpred = np.amax(gdat.psdn) > 1.

    else:
        boolsignpred[k] = False
    
    titl = 'Classified as '
    if boolsignpred[k]:
        titl += 'BHC candidate'
    else:
        titl += 'background'
    
    path = gdat.pathdata + 'lcur%04d.pdf' % k
    plot_lcur(path, gdat.lcurdata, titl=titl)
    
    path = gdat.pathdata + 'psdn%04d.pdf' % k
    plot_psdn(path, gdat.psdn, titl=titl)
        
