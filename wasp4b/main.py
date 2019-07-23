"""
Analysis of RV residuals

Expects $WASP4B_DATA_PATH to be set to a path where the input data exists in $WASP4B_DATA_PATH/data/
Generated output plots will be generated under $WASP4B_DATA_PATH/plot/

Tansu Daylan
MIT Kavli Institute, Cambridge, MA, 02109, US
tansu.daylan@gmail.com
www.tansudaylan.com
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy

import scipy 
import os
from astropy.time import Time
import emcee

def summgene(arry):
    
    print np.amin(arry)
    print np.amax(arry)
    print np.mean(arry)
    print arry.shape


def retr_radvmodl(ampl, phasoffs, freqderi, time):
    
    radvmodl = ampl * np.cos(2. * np.pi * (freq * time + 0.5 * freqderi * time**2) + 2. * np.pi * phasoffs)
    #print 'phasoffs'
    #print phasoffs
    #print 'freq * time'
    #print freq * time
    #print 
    #print 'radvmodl'
    #print radvmodl

    return radvmodl


def retr_listresi(para, phasobsv, numbobsv, radvobsv, radvstdvobsv):
    
    ampl = para[0]
    phasoffs = para[1]
    if peritype == 'deca':
        freqderi = para[2]
        totloffs = para[3]
    else:
        offs = para[2:]

    listresi = []
    for kk, k in enumerate(indxdataincl):
        if peritype == 'deca':
            radvmodl = retr_radvmodl(ampl, phasoffs, freqderi, listtimeobsv[k])
            listresi.append(radvmodl - listradv[k] - totloffs)
        if peritype == 'cons':
            radvmodl = retr_radvmodl(ampl, phasoffs, 0., listtimeobsv[k])
            listresi.append(radvmodl - listradv[k] - offs[kk])
    
    return listresi


def retr_lpos(para, phasobsv, numbobsv, radvobsv, radvstdvobsv):

    listresi = retr_listresi(para, phasobsv, numbobsv, radvobsv, radvstdvobsv)
    
    if ((para < limtpara[0, :]) | (para > limtpara[1, :])).any():
        lpos = -np.inf
    else:
        llik = 0.
        for kk, k in enumerate(indxdataincl):
            llik += np.sum(-0.5 * listresi[kk]**2 / listeadv[k]**2)
        lpos = llik

    return lpos


def retr_chi2doff(paramlik, phasobsv, numbobsv, radvobsv, radvstdvobsv):
    
    llik = retr_llik(paramlik, phasobsv, numbobsv, radvobsv, radvstdvobsv)
    
    chi2doff = -2. * llikmlik / (numbobsv - numbpara)
     
    return chi2doff


# initialization
pathdata = os.environ['WASP4B_DATA_PATH'] + '/'
pathalle = pathdata + 'data/'
pathplot = pathdata + 'plot/'
os.system('mkdir -p %s' % pathplot)

# read the data
listlabldata = ['HARPS (Husnoo+2012)', 'CORALIE (Wilson+ 2008)', 'CORALIE (Triaud+2010)', 'HARPS (Triaud+2010)', 'HIRES (Knutson+2014)']
listfiledata = ['HARPS_HUSN', 'CORALIE_WILS', 'CORALIE_TRIA', 'HARPS_TRIA', 'HIRES']
liststrgdata = ['harh', 'corw', 'cort', 'hart', 'hire']
numbdata = len(liststrgdata)
indxdata = np.arange(numbdata)
listtimeobsv = [[] for k in indxdata]
listphasobsv = [[] for k in indxdata]
listradv = [[] for k in indxdata]
listeadv = [[] for k in indxdata]

peri = 1.33823204
freq = 1. / peri

cntr = 0
indxobsvdata = {}
for k, strgdata in enumerate(liststrgdata):
    arry = np.loadtxt(pathalle + listfiledata[k] + '.csv', delimiter=',')
    listtimeobsv[k] = arry[:, 0]
    # take out the mean of each data set
    listradv[k] = arry[:, 1]# - np.mean(arry[:, 1])
    if strgdata.startswith('har') or strgdata.startswith('cor'):
        listradv[k] -= 57.7326
    listeadv[k] = arry[:, 2]
    listradv[k] *= 1e3
    listeadv[k] *= 1e3
    
    indxobsvdata[strgdata] = np.arange(cntr, cntr + arry[:, 2].size)
    cntr += arry[:, 2].size

    listphasobsv[k] = (listtimeobsv[k] / peri) % 1.

timeobsv = np.concatenate(listtimeobsv)
phasobsv = np.concatenate(listphasobsv)
radvobsv = np.concatenate(listradv)
radvstdvobsv = np.concatenate(listeadv)

indxtimeobsvsort = np.argsort(timeobsv)
timeobsv = timeobsv[indxtimeobsvsort]
phasobsv = phasobsv[indxtimeobsvsort]
radvobsv = radvobsv[indxtimeobsvsort]
radvstdvobsv = radvstdvobsv[indxtimeobsvsort]


minmtimeobsv = np.amin(timeobsv)
maxmtimeobsv = np.amax(timeobsv)

# plot the data
## raw RV curve
figr, axis = plt.subplots(figsize=(12, 6))
for k, strgdata in enumerate(liststrgdata):
    ydat = listradv[k]
    temp, listcaps, temp = axis.errorbar(listtimeobsv[k], ydat, yerr=listeadv[k], marker='o', ls='', markersize=0.5, label=listlabldata[k], capsize=1, lw=0.5)
    for caps in listcaps:
        caps.set_markeredgewidth(0.5)
axistwin = axis.twiny()
axistwin.set_xlim(Time(axis.get_xlim(), format='jd', scale='utc').decimalyear)
axis.set_xlabel('Time [BJD]')
axis.set_ylabel('RV [m/s]')
axis.legend()
plt.tight_layout()
plt.savefig(pathplot + 'radvraww.pdf')
plt.close()


## phase-folded RV curve
figr, axis = plt.subplots(figsize=(12, 6))
for k, strgdata in enumerate(liststrgdata):
    ydat = listradv[k]
    temp, listcaps, temp = axis.errorbar(listphasobsv[k], ydat, yerr=listeadv[k], marker='o', ls='', markersize=0.5, label=listlabldata[k], capsize=1, lw=0.5)
    for caps in listcaps:
        caps.set_markeredgewidth(0.5)
axis.set_xlabel('Phase')
axis.set_ylabel('RV [m/s]')
plt.legend()
plt.tight_layout()
plt.savefig(pathplot + 'radvphas.pdf')
plt.close()


listdictcnfg = [ \
                ['cons', 'fulldata'], \
                ['cons', 'husnknut'], \
                ['cons', 'husn'], \
                ['cons', 'husncora'], \
                #['deca', 'husncora'], \
               ]

for dictcnfg in listdictcnfg:
    
    peritype = dictcnfg[0]
    strgcnfg = dictcnfg[1]

    #listlabldata = ['HARPS (Husnoo+2012)', 'CORALIE (Wilson+ 2008)', 'CORALIE (Triaud+2010)', 'HARPS (Triaud+2010)', 'HIRES (Knutson+2014)']
    
    if strgcnfg == 'fulldata':
        indxdataincl = range(numbdata)
    if strgcnfg == 'husnknut':
        indxdataincl = [0, 4]
    if strgcnfg == 'husn':
        indxdataincl = [0]
    if strgcnfg == 'husncora':
        indxdataincl = [0, 1, 2]
    
    rtag = '%s_%s' % (peritype, strgcnfg)

    numbpara = 2
    if peritype == 'deca':
        numbpara += 2
    else:
        numbpara += len(indxdataincl)
    numbsampwalk = 100
    numbsampburn = 20000
    numbwalk = 50
    numbsamp = numbsampwalk * numbwalk
    
    indxsamp = np.arange(numbsamp)
    indxsampwalk = np.arange(numbsampwalk)
    
    indxpara = np.arange(numbpara)
    limtpara = np.empty((2, numbpara))
    # offs
    limtpara[0, 0] = 0.
    limtpara[1, 0] = 1000.
    limtpara[0, 1] = 0.
    limtpara[1, 1] = 0.2
    if peritype == 'deca':
        limtpara[0, 2] = -1e-13
        limtpara[1, 2] = 1e-13
        limtpara[0, 3] = -100.
        limtpara[1, 3] = 100.
    else:
        for k in range(len(indxdataincl)):
            limtpara[0, k+2] = -100.
            limtpara[1, k+2] = 100.
    
    indxwalk = np.arange(numbwalk)
    parainit = []
    for k in indxwalk:
        parainit.append(np.empty(numbpara))
        meannorm = (limtpara[0, :] + limtpara[1, :]) / 2.
        stdvnorm = (limtpara[0, :] - limtpara[1, :]) / 10.
        parainit[k]  = (scipy.stats.truncnorm.rvs((limtpara[0, :] - meannorm) / stdvnorm, (limtpara[1, :] - meannorm) / stdvnorm)) * stdvnorm + meannorm
    numbsampwalk = numbsamp / numbwalk
    numbsampwalkburn = numbsampburn / numbwalk
    
    numbobsv = timeobsv.size
    indxobsv = np.arange(numbobsv)
    
    dictllik = [phasobsv, numbobsv, radvobsv, radvstdvobsv]
    
    objtsamp = emcee.EnsembleSampler(numbwalk, numbpara, retr_lpos, args=dictllik, threads=3)
    parainitburn, prob, state = objtsamp.run_mcmc(parainit, numbsampwalkburn)
    objtsamp.reset()
    objtsamp.run_mcmc(parainitburn, numbsampwalk)
    
    ## parameter
    ### trace
    for k in indxpara:
        
        figr, axis = plt.subplots()
        for i in indxwalk:
            axis.plot(indxsampwalk, objtsamp.chain[i, :, k])
        path = pathplot + '%s_tracwalk_%04d.pdf' % (rtag, k)
        print 'Writing to %s...' % path
        plt.savefig(path)
        plt.close()
            
    ### histogram
    for k in indxpara:
        figr, axis = plt.subplots()
        axis.hist(objtsamp.flatchain[:, k]) 
        path = pathplot + '%s_hist_%04d.pdf' % (rtag, k)
        print 'Writing to %s...' % path
        plt.savefig(path)
        plt.close()
    
    ## log-likelihood
    figr, axis = plt.subplots()
    for i in indxwalk:
        axis.plot(indxsampwalk, objtsamp.lnprobability[:, i])
    path = pathplot + '%s_llik.pdf' % rtag
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()
    
    ### sample model radv
    numbradvmodl = 100
    indxradvmodl = np.arange(numbradvmodl)
    indxsamprand = np.random.choice(indxsamp, numbradvmodl, replace=False)
    yerr = np.empty((2, numbobsv))
    yerr[0, :] = radvstdvobsv
    yerr[1, :] = radvstdvobsv
    numbobsvfine = 100
    #phasfine = np.linspace(0., 1., numbobsvfine)
    #timefine = minmtimeobsv + phasfine * peri
    timefine = np.linspace(minmtimeobsv - 100., maxmtimeobsv + 100., numbobsvfine)
    phasfine = (timefine / peri) % 1.
    indx = np.argsort(phasfine)
    phasfine = phasfine[indx]
    timefine = timefine[indx]
    indxmlik = np.argmax(objtsamp.lnprobability)
    indxmlik = np.unravel_index(indxmlik, objtsamp.lnprobability.shape)
    paramlik = objtsamp.chain[indxmlik[1], indxmlik[0], :]
    llikmlik = retr_lpos(paramlik, phasobsv, numbobsv, radvobsv, radvstdvobsv)
    chi2doff = -2. * llikmlik / (numbobsv - numbpara)
    
    rmsq = {}
    rmsq['totl'] = np.sqrt(np.mean(np.concatenate(retr_listresi(paramlik, phasobsv, numbobsv, radvobsv, radvstdvobsv))**2))

    if peritype != 'deca':
        radvmodlfine = np.empty((numbsamp, numbobsvfine))
        for k in indxradvmodl:
            ampl = objtsamp.flatchain[indxsamprand[k], 0]
            phasoffs = objtsamp.flatchain[indxsamprand[k], 1]
            radvmodlfine[k, :] = retr_radvmodl(ampl, phasoffs, 0., timefine)
        
        ## phase-folded RV curve
        figr, axis = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw = {'height_ratios':[2, 1]})
        ampl = paramlik[0]
        phasoffs = paramlik[1]
        for k in indxradvmodl:
            axis[0].plot(phasfine, radvmodlfine[k, :], alpha=0.05, color='b')
        for kk, k in enumerate(indxdataincl):
            ydat = listradv[k]
            
            if peritype == 'deca':
                if strgdata == 'knut':
                    ydattemp = ydat + paramlik[3] + paramlik[4]
                else:
                    ydattemp = ydat + paramlik[4]
            else:
                ydattemp = ydat + paramlik[2+kk]
            
            temp, listcaps, temp = axis[0].errorbar(listphasobsv[k], ydattemp, yerr=listeadv[k], marker='o', ls='', markersize=0.5, \
                                                                                                    label=listlabldata[k], capsize=1, lw=0.5)
            for caps in listcaps:
                caps.set_markeredgewidth(0.5)
            radvmodlmlik = retr_radvmodl(ampl, phasoffs, 0., listtimeobsv[k])
            temp, listcaps, temp = axis[1].errorbar(listphasobsv[k], ydattemp - radvmodlmlik, yerr=listeadv[k], marker='o', ls='', \
                                                                                                    markersize=0.5, label=listlabldata[k], capsize=1, lw=0.5)
            for caps in listcaps:
                caps.set_markeredgewidth(0.5)
        axis[0].set_ylabel('RV [m/s]')
        axis[0].set_xticklabels([])
        axis[1].axhline(0, ls='--', alpha=0.2, color='black')
        axis[1].set_xlabel('Phase')
        axis[1].set_ylabel('Residual RV [m/s]')
        axis[0].text(0.03, 0.2, 'RMS = %.4g' % rmsq['totl'], transform=axis[0].transAxes)
        axis[0].legend(loc=4)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.)
        path = pathplot + '%s_radvphasmodl.pdf' % rtag
        print 'Writing to %s...' % path
        plt.savefig(path)
        plt.close()
    
    print


