
# coding: utf-8

# In[70]:

get_ipython().run_cell_magic(u'javascript', u'', u'IPython.OutputArea.auto_scroll_threshold = 9999;')


# In[71]:

get_ipython().magic(u'matplotlib inline')
get_ipython().magic(u"config InlineBackend.figure_format = 'retina'")


# In[72]:

import os, time

from numpy import *
from numpy.random import *
from numpy.random import choice

import scipy as sp
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline

import pyfits as pf

import tdpy_util

# plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)


# In[73]:

# determine the thermal state of the IGM due to photoheating
## blazar heating rate
    
    #plheatbl = 10**(0.0315 * (1. + red)**3 - 0.512 * (1. + red)**2 + 2.27 * (1. + red) - 2.38) / 3.154e22 # [MeV/s]
    #plheatbl[where(reds gt 5.7)] = 0.
    
    
## photoheating rate
    
    #  phheatdm = zeros(nreds)
    #  phheatqg = zeros(nreds)
    #  phheatqghm = zeros(nreds)
    #  phheatqs = zeros(nreds)
    #  for c=0, nred-1 do begin
    #    phheatdm[c] = trapz(enph, pucr_phheatint(phdmflux[:,c]), xr=[ryd, enph[nenph-1]]) # [MeV/s]
    #    phheatqg[c] = trapz(enph, pucr_phheatint(phqgflux[:,c]), xr=[ryd, enph[nenph-1]]) # [MeV/s]
    #    phheatqghm[c] = trapz(enph, pucr_phheatint(phqgfluxhm[:,c]), xr=[ryd, enph[nenph-1]]) # [MeV/s]
    #    phheatqs[c] = trapz(enph, pucr_phheatint(phqsflux[:,c]), xr=[ryd, enph[nenph-1]]) # [MeV/s]
    #  endfor
    

    #if makeplot:
    #    puc_heat_plot
    
    
    # density - temperature relation
    
    #  npatch = 100
    #  initreds = 20.
    #  mindif = min(abs(reds - initred), jredtemp)
    #  temp = cmbt
    #
    #  puc_lambda
    #  heat = phheatdm
    #  temppatchdm = zeros((npatch, nreds))
    # dpatchdensdzarr = zeros((npatch, nreds))
    #  patchdensarr = zeros((npatch, nreds))
    #  for i=0, npatch-1 do begin
    #    patchdensarr[i,:] = 1. / ((1. + lambda[i,0] * grwf) * (1. + lambda[i,1] * grwf) * (1. + lambda[i,2] * grwf)) - 1.
    #    dpatchdensdzarr[i,:] = deriv(red, patchdensarr[i,:])
    #    patchdens = patchdensarr[i,:]
    #    dpatchdensdz = dpatchdensdzarr[i,:]
    #    temppatchdm[i,:] = puc_temp_solve
    #  endfor    

    #if makeplot:
    #    pucr_patch_plot 
    #
    #
    #  # thermal evolution of the IGM at the mean overdensity patch
    #
    #  patchdens = patchdensarr[0,:]
    #  dpatchdensdz = dpatchdensdzarr[0,:]
    #
    #  ionfrac = 1e-6
    #
    #  heat = phheatdm * ionfrac
    #  tempdm = puc_temp_solve()
    #
    #  heat = phheatqs * ionfrac
    #  tempqs = puc_temp_solve()
    #
    #  heat = plheatbl * ionfrac
    #  tempbl = puc_temp_solve()
    #
    #if makeplot:
    #    puc_temp_plot 


# In[74]:

def eden_egbl():
    
    path = os.environ["PHOT_IONZ_DATA_PATH"] + '/egbl.csv'
    egbldata = loadtxt(path)
    wlenegbl = egbldata[:, 0] * 1e-4 # [cm]
    freqegbl = flipud(velolght / wlenegbl) # [Hz]
    enpiegbl = plnkcons * freqegbl # [MeV]
    fluxegbl = egbldata[:, 1] # [W/m^2/sr]
    edenegbltemp = flipud(fluxegbl) * 2. * pi / velolght / 1e4 * 1.6e13 / enpiegbl**2 # [1/cm^3/MeV]
    
    edenegbl = zeros((nenpi, nreds))
    for c in range(nreds):
        enpitemp = enpi / (1. + reds[c]) # [MeV]
        jenpi = where((enpitemp < max(enpiegbl)) & (enpitemp > min(enpiegbl)))
        edenegbl[jenpi,c] = interp1d(enpiegbl, edenegbltemp)(enpitemp[jenpi]) * (1. + reds[c])**2 # [1/cm^3/MeV]

    return edenegbl
        


def eden_cmbr():
    
    edencmbr = 8. * pi * enpi[:, None]**2 / (velolght * plnkcons)**3 / (exp(enpi[:,None] / tempcmbrnunc / boltcons / (1. + reds[None, :])) - 1.) # [1/cm^3/MeV]
    
    return edencmbr


# In[75]:

def edot_plot():

    sec2gyr = 3.171e-17 # [s/Gyrs]
    
    ndengas = [0.01, 1., 100.] # [1/cm^3]
    edengas = [1e-8, 1e-6, 1e-4] # [MeV/cm^3]
    intsmag = [1e-10, 1e-8, 1e-6] # [muG]
    
    time = enel[:,None] / retr_edot(intsmag, ndengas, edenrad) * sec2gyr
    timediff = zscldif**2 / diffnor / (enel * 1e-3)**diffind
    

    labels = []
    for a in range(nndengas):
        labels.append([r'Brem, $\rho_{gas}$ = %.2g  cm$^{-3}$' % ndengas[a]])
    for a in range(nedenrad):
        labels.append(['ICS, $n_{isrf}$ = %.2g  1/cm$^{-3}$' % edenrad[a]])
    for a in range(nintsmag):
        labels.append([r'Synch, $B_$ = %.2g $\mu$G' % intsmag[a]])
        
    fig, ax = plt.subplots(figsize=(7,7))
    ax.loglog(enel, time)
    ax.set_xlabel('$E_e$ [MeV]')
    ax.set_ylabel(r'$\tau$ [Gyrs]')
    ax.legend()
    
    plt.savefig(plotpath + 'edotelec.png')
    plt.close()


# In[76]:

def retr_psec(thiswnum):

    q = thiswnum / 0.15
    cq = 14.4 + 325 / (1. + 60.5 * q**1.11)
    lq = log(exp(1.) + 1.84 * q)
    tranfunc = lq / (lq + cq * q**2)
    psecprim = 2e-9 * thiswnum**psecindx
    
    psec = psecprim * tranfunc**2
    
    return psec, tranfunc


# In[77]:

def retr_llik(sampvarb):

    csecvelo = sampvarb[0]
    csecfrac = sampvarb[1]
    masspart = sampvarb[2]
    dmatslop = sampvarb[3]
    
    
    fluxphotdmatintp = fluxphotdmat[:, :, 0] * csecvelo / csecvelopivt + fluxphotdmat[:, :, 1] * csecfrac / csecfracpivt
    
    
    fluxphotmodl = fluxphothm12[fenph, 0] + fluxphotdmatintp[fenph, 0]
    lpos = sum(-log(sqrt(2. * pi * fluxphotexprvari)) - 0.5 * (fluxphotexpr - fluxphotmodl)**2 / fluxphotexprvari)
    
    
    if False:
        
        print 'retr_llik'
        print 'csecvelo'
        print csecvelo
        print 'csecfrac'
        print csecfrac
        print 'csecvelo / csecvelopivt'
        print csecvelo / csecvelopivt
        print 'csecfrac / csecfracpivt'
        print csecfrac / csecfracpivt
        print 'fluxphotdmat[fenph, 0, 0]'
        print fluxphotdmat[fenph, 0, 0]
        print 'fluxphotdmat[fenph, 0, 1]'
        print fluxphotdmat[fenph, 0, 1]
        print 'fluxphotdmatintp'
        print fluxphotdmatintp
        print 'fluxphothm12[fenph, 0]'
        print fluxphothm12[fenph, 0]
        print 'fluxphotmodl'
        print fluxphotmodl
        print 'fluxphothm12[fenph, 0]'
        print fluxphothm12[fenph, 0]
        print 'fluxphotexpr'
        print fluxphotexpr
        print 'fluxphotexprvari'
        print fluxphotexprvari
        print 'lpos'
        print lpos
        print
        print
        print
        print
        print
        print
        print
        print
        
        
    sampcalc = [fluxphotdmatintp]
    
    return lpos, sampcalc


# In[78]:

def retr_hmfn():
      
    # growth factor
    grwf = zeros(nreds) # growth factor
    for c in range(nreds):
        diffgrwfdiffreds = funchubb[c] * (1. + reds) / funchubb**3
        grwf[c] = trapz(diffgrwfdiffreds[c:], reds[c:])
    grwf /= grwf[0]
    
    # radius, wavelength and wavenumber corresponding to the halo mass
    rsphhalo = (3. * massp * solm2mgev / 4. / pi / omegdmat / edencritnunc / odenviri)**(1./3.) / kprc2cmet / 1e3 # [Mpc]
    wlenhalo = 4. * rsphhalo
    wnumhalo = 2. * pi / wlenhalo

    # power spectrum of density fluctuations
    psec, tranfunc = retr_psec(wnum)

    # RMS density fluctuations
    fluc = zeros(nmass + 1)
    diffflucdiffwnum = zeros((nwnum, nmass + 1))
    funcwndw = zeros((nwnum, nmass + 1))
    for d in range(nmass + 1):
        wang = wnum * wlenhalo[d]
        funcwndw[:, d] = 3. * (sin(wang) - wang * cos(wang)) / wang**3
        diffflucdiffwnum[:, d] = wnum**2 * psec * funcwndw[:, d]**2 / 2. / pi**2
        fluc[d] = sqrt(trapz(diffflucdiffwnum[:, d], wnum, axis=0))
    # temp
    fluc *= 0.55 / interp1d(massp, fluc)(1e15)
    #fluc *= odenrmsq8mpc / interp1d(rsphhalo, fluc)(8. / hubbcons)
    fluc = fluc[:, None] * grwf[None, :]

    # halo mass function
    difflogtflucdiffmass = -diff(log(fluc), axis=0) / diff(massp)[:, None]
    peakhght = odencoll / fluc[:-1, :]
    funcfluc = shtrnorm * sqrt(2. * shtrwgth / pi) * (1. + (1. / peakhght**2 / shtrwgth)**shtrindx) *         peakhght * exp(-shtrwgth * peakhght**2 / 2.)
        
    # temp
    fudgfact = 1.8
    diffnhaldiffmass = fudgfact * funcfluc * edendmat[None, :] * kprc2cmet**3 / solm2mgev / mass[:, None] *         difflogtflucdiffmass # [1/kpc^3/Msun]
        
    if makeplot:
        
        fig, ax = plt.subplots()
        ax.loglog(reds, grwf)
        ax.set_xlabel(r'$z$')
        ax.set_ylabel('$D(z)$')
        ax.set_title('Growth factor')
        ax.legend()
        plt.savefig(plotpath + 'grwf.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.loglog(massp, rsphhalo)
        ax.set_xlabel(r'$M [M_\odot]$')
        ax.set_ylabel('$r_H$ [Mpc]')
        ax = ax.twinx()
        ax.loglog(massp, wnumhalo, ls='--')
        ax.set_ylabel('$k_H$ [Mpc$^{-1}$]')
        ax.legend()
        plt.savefig(plotpath + 'rsphhalo.png')
        plt.close()
        
        
        fig, ax = plt.subplots()
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('P(k) [Mpc$^3$]')
        ax.loglog(wnum, psec)
        ax.set_title('Primordial Matter Power Spectrum')
        plt.savefig(plotpath + 'psec.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_title('Transfer function')
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('$T(k)$')
        ax.loglog(wnum, tranfunc)
        plt.savefig(plotpath + 'tranfunc.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_title('Window function')
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('$W(kR)$')
        for d in range(njmass):
            ax.loglog(wnum, funcwndw[:, jmass[d]]**2, label=strgmass[d])
        ax.legend(loc=3)
        plt.savefig(plotpath + 'funcwndw.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_title('Contribution of spatial scales to the RMS density fluctuations')
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('$d\sigma^2/dk(M)$')
        for d in range(njmass):
            ax.loglog(wnum, diffflucdiffwnum[:, jmass[d]], label=strgmass[d])
        ax.legend(loc=3)
        plt.savefig(plotpath + 'diffflucdiffwnum.png')
        plt.close()
        
        
        fig, ax = plt.subplots()
        for c in range(njredslate):
            ax.loglog(mass, difflogtflucdiffmass[:, jredslate[c]], label=strgredslate[c])
        ax.set_xlabel(r'$M [M_\odot]$')
        ax.set_ylabel(r'd$\log\sigma$/d$M$')
        ax.legend()
        plt.savefig(plotpath + 'difflogtflucdiffmass.png')
        plt.close()
        
        
        fig, ax = plt.subplots()
        for c in range(njredslate):
            ax.loglog(massp, fluc[:, jredslate[c]], label=strgredslate[c])
        ax.set_xlabel(r'$M [M_\odot]$')
        ax.set_ylabel(r'$\sigma$')
        ax.axhline(odencoll, ls='--', label='Critical linear overdensity at collapse')
        ax.legend()
        plt.savefig(plotpath + 'fluc.png')
        plt.close()

        
        
        fig, ax = plt.subplots()
        ax.loglog(fluc[:-1, 0], funcfluc[:, 0])
        ax.set_xlabel(r'$\sigma$')
        ax.set_ylabel('$f(\sigma)$')
        plt.savefig(plotpath + 'funcfluc.png')
        plt.close()

    
        datamsm1 = loadtxt(os.environ["PHOT_IONZ_DATA_PATH"] + '/msm1.csv')
        datamsm2 = loadtxt(os.environ["PHOT_IONZ_DATA_PATH"] + '/msm2.csv')
        

        fig, ax = plt.subplots()
        fig.suptitle('Halo mass function', fontsize=18)
        ax.errorbar(datamsm1[:, 0], datamsm1[:, 1] / datamsm1[:, 0], ls='',                     yerr=sqrt(datamsm1[:, 1] / datamsm1[:, 0]), # / (500 / 0.7)**3), \
                    marker='o', markersize=5, color='black', label=r'MS-I, $z = 0$')
        ax.errorbar(datamsm2[:, 0], datamsm2[:, 1] / datamsm2[:, 0], ls='',                     yerr=sqrt(datamsm2[:, 1] / datamsm2[:, 0]), # / (100 / 0.7)**3), \
                    marker='D', markersize=5, color='grey', label=r'MS-II, $z = 0$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        for c in range(njredslate):
            plot = ax.loglog(mass, mass * diffnhaldiffmass[:, jredslate[c]] * 1e9, label=strgredslate[c])
        ax.set_ylim([1e-9, 1e5])
        ax.set_ylabel('$dN_h/d\log M$ [1/Mpc$^3$]')
        ax.set_xlabel('$M [M_\odot]$')     
        ax.legend()
        plt.savefig(plotpath + 'diffnhaldiffmass.png')
        plt.close()



    return diffnhaldiffmass, peakhght


# In[79]:

def retr_fluxphothm12():
    
    # Haardt & Madau 2012 quasar + galaxy UV/X-ray background flux
    name = os.environ["PHOT_IONZ_DATA_PATH"] + '/photfluxhm12.dat'
    tabl = loadtxt(name)

    wlenhm12 = tabl[1:576, 0] # [A]
    freqhm12 = velolght / flipud(wlenhm12) * cmet2angs # [Hz]
    enphhm12 = plnkcons * freqhm12 # [MeV]
    redshm12 = tabl[0, 0:60]
    fluxphothm12 = flipud(tabl[1:576, 1:61]) # [erg/cm^2/s/Hz/sr]
    fluxphothm12 *= ergs2mgev / plnkcons / enphhm12[:, None] # [1/cm^2/s/sr/MeV]
    
    enphhm12, jenphhm12 = unique(enphhm12, return_index=True)
    fluxphothm12 = fluxphothm12[jenphhm12, :]

    fluxphothm12 = interp2d(enphhm12, redshm12, fluxphothm12.T)(enph, reds).T

    return fluxphothm12
    
    # quasar photon emissivity
    #  phqsemis1 = 10**24.6 * erg2mev / (kprc2cmet*1d3)**3 * (1. + reds)**7.68 * exp(-0.28 * reds) / (exp(1.77 * reds) + 26.3) # [1/cm^3/s/Hz/MeV]
    #  phqsemis = ((enph / enerrydb)**(-0.44) * frph) # phqsemis1 # [/cm^3/s/MeV]
    #  ienph = where(enph gt plnkcons * velolght / 1.3e-5)
    #  phqsemis[ienph,:] = ((enph[ienph] / enerrydb)**(-1.57) * frph) # phqsemis1 # [1/cm^3/s/MeV]
    #
    #
    # Haardt & Madau 2012 quasar + galaxy UV/X-ray background emissivity

    #hmemis = read_ascii('$PHOT_IONZ_DATA_PATH/dat/qgemis.dat') 
    #wlemishm = phqgemisdat.field01[0,1:378] # [Angstrom]
    #fremishm = reverse(reform(velolght / (wlemishm * 1e-8))) # [Hz]
    #enemishm = plnkcons * fremishm # [MeV]
    #phqgemis = phqgemisdat.field01[1:59,1:378] # [erg/Mpc^3/s/Hz]
    #phqgemis = transpose(phqgemis)
    #phqgemis = reverse(phqgemis,1)
    #phqgemis *=  2.12e-68 * fremishm * (1. + zeros(60)) # [1/cm^3/s/MeV]
    #phemis[:,:,1] = interpolate(phqgemis, interpol(indgen(378), enemishm, enph), interpol(indgen(59), redhm, red)) 
    #phemisryd[:,1] = interpolate(phqgemis, interpol(indgen(378), enemishm, 1. * enerrydb), interpol(indgen(59), redhm, red)) 


# In[80]:

def plot_matr(ax, xdat, ydat, labl, loc=1):
    
    listlinestyl = [':', '--', '-']
    listcolr = ['b', 'r', 'g']
    
    for i in range(3):
        for  j in range(3):
            if len(xdat.shape) == 3:
                ax.plot(xdat[i, j, :], ydat[i, j, :], color=listcolr[j], ls=listlinestyl[i])
            else:
                ax.plot(xdat, ydat[i, j, :], color=c[j], ls=ls[i])

    line = []
    for k in arange(3):
        line.append(plt.Line2D((0,1),(0,0), color='k', ls=listlinestyl[k]))
    for l in range(3):
        line.append(plt.Line2D((0,1),(0,0), color=listcolr[l]))
    ax.legend(line, labl, loc=loc, ncol=2) 


# In[81]:

def retr_fluxphotdmat(csecvelo, csecfrac, masspart, dmatslop):
    
    # run tag
    rtag = propmodl + '_' + concmodl + '_' + subsmodl + '_' + igmamodl +     '_massandm%.3g_' % masspart + anch + '_csecvelo%.3g'  % -log10(csecvelo) +     '_csecfrac%.3g' % -log10(csecfrac)
    
    # DM annihilation spectrum
    multintptemp = interp1d(massp4dm, multp4dm, axis=1)(masspart) / masspart / enelscalp4dm
    
    menel = where((enel > amin(enelscalp4dm * masspart)) & (enel < amax(enelscalp4dm * masspart)))[0]
    multintp = zeros(nenel)
    multintp[menel] = interp1d(enelscalp4dm * masspart, multintptemp)(enel[menel])
    
    
    # energy density and velocity variance in DM halos
    global edendmathalo, velovarihalo, rvelvarihalo
    velovarihalo = zeros((nmass, nreds, nrsph))
    edendmathalo = zeros((nmass, nreds, nrsph))
    for d in range(nmass):
        for c in range(nreds):
            edendmathalonorm = odenviri * edendmat[c] *                 conc[d, c]**3 / 3. / (log(1. + conc[d, c]) - conc[d, c] / (1. + conc[d, c])) # [MeV/cm^3]
            rsphscal = rsphviri[d, c] / conc[d, c] # [kpc]
            rsphnorm = rsph[d, c, :] / rsphscal # [1]
            edendmathalo[d, c, :] = edendmathalonorm / rsphnorm / (1. + rsphnorm)**2 # [MeV/cm^3]
            
            edendemcnorm = (2. * demccons - 3.)**3 / 4. / (5. - 2. * demccons) * edendmathalonorm # [MeV/cm^3]
            rsphdemc = (5. - 2. * demccons) / (2. * demccons - 3.) * rsphscal # [kpc]
            velovarihalo[d, c, :] = 4. * pi * gravconsredu / 200. * 81. * edendemcnorm**(1. / 3.) *                 rsphdemc**2 * (edendmathalo[d, c, :] * (rsph[d, c, :] / rsphdemc)**demccons)**(2. / 3.) /                 velolght**2 * kprc2cmet**2 # [1] 
    rvelvarihalo = demcsigm * velovarihalo
          
    # DM annihilation cross section
    csecvelohalo = zeros((nmass, nreds, nrsph, nmpol))
    csecvelohalo[:, :, :, 0] = csecvelo
    csecvelohalo[:, :, :, 1] = csecvelo * csecfrac * rvelvarihalo


    # electron emissivity in DM halo
    global emiselechaloenel, diffemiselechalodiffrsphenel, lumielechaloenel,         diffemiselecclmpdiffmassenel, emiselecenel, emiselecclmpenel,         emiselecsmthenel, rvelvarihavg,         rvelvariiavg, rvelvariiavgsmth, rvelvariiavgclmp,         ndenelecsavg, ndenelecsavgenel, emisphot,         emisphotrydb, fluxphot, fluxphotrydb, fluxphotenph, diffionrdiffenph, ionr
        
        
    emiselechalohost = zeros((nenel, nmass, nreds, nrsph, nmpol))
    emiselechalosubh = zeros((nenel, nmass, nreds, nrsph, nmpol))
    emiselechalotemp = zeros((nenel, nmass, nmass, nreds, nrsph, nmpol))
    emiselechalo = zeros((nenel, nmass, nreds, nrsph, nmpol))
    diffemiselechalodiffrsph = zeros((nenel, nmass, nreds, nrsph, nmpol))
    lumielechalo = zeros((nenel, nmass, nreds, nmpol))
    for a in range(nenel):
        for d in range(nmass):
            for c in range(nreds):
                for m in range(nmpol):
                    emiselechalohost[a, d, c, :, m] = csecvelohalo[d, c, :, m] / masspart**2 *                         multintp[a] * fesc[a, d, c, :] * edendmathalo[d, c, :]**2 # [1/s/cm^3/MeV]

                    for f in range(nmass):
                        if f < d:
                            emiselechalotemp[a, d, f, c, :, m] = 0.
                    jmass = where(mass[f] < mass[d])[0]
                    emiselechalosubh[a, d, c, :, m] = trapz(emiselechalotemp[a, d, :, c, :, m], mass, axis=0)
                  
                    emiselechalo[a, d, c, :, m] = emiselechalohost[a, d, c, :, m] + emiselechalosubh[a, d, c, :, m]
                    
                    diffemiselechalodiffrsph[a, d, c, :, m] = 4. * pi * rsph[d, c, :]**2 *                         emiselechalo[a, d, c, :, m] * kprc2cmet**2 # [1/s/kpc/MeV]
                    lumielechalo[a, d, c, m] =                         trapz(diffemiselechalodiffrsph[a, d, c, :, m] * kprc2cmet, rsph[d, c, :]) # [1/s/MeV]

    emiselechaloenel = trapz(emiselechalo, enel, axis=0) # [1/s/cm^3]
    lumielechaloenel = trapz(lumielechalo, enel, axis=0) # [1/s]
    diffemiselechalodiffrsphenel = trapz(diffemiselechalodiffrsph, enel, axis=0) # [1/s/cm^3]
    

    
    
    # mean DM relative velocity
    diffrvelvarihavgdiffrsph = zeros((nmass, nreds, nrsph))
    rvelvarihavg = zeros((nmass, nreds))
    for d in range(nmass):
        for c in range(nreds):
            diffrvelvarihavgdiffrsph[d, c, :] = 3. * rsph[d, c, :]**2 * rvelvarihalo[d, c, :] / amax(rsph[d, c, :])**3
            rvelvarihavg[d, c] = trapz(diffrvelvarihavgdiffrsph[d, c, :], rsph[d, c, :]) # [1]

    # spatially averaged electron emissivity    
    ## smooth component
    # temp -- add p-wave to the smth component
    emiselecsmth = zeros((nenel, nreds, nmpol))
    for c in range(nreds):
        for m in range(nmpol):
            emiselecsmth[:, c, m] = csecvelo * edendmat[c]**2 / masspart**2 * multintp # [1/cm^3/s/MeV]
    emiselecsmthenel = trapz(emiselecsmth, enel, axis=0) # [1/cm^3/s]
    

    ## clumpy component
    diffemiselecclmpdiffmass = zeros((nenel, nmass, nreds, nmpol))
    for a in range(nenel):
        for d in range(nmass):
            for c in range(nreds):
                for m in range(nmpol):
                    diffemiselecclmpdiffmass[a, :, c, m] = diffnhaldiffmass[:, c] *                         lumielechalo[a, :, c, m] / kprc2cmet**3 * (1. + reds[c])**3 # [1/cm^3/s/MeV/Msun]
    emiselecclmp = trapz(diffemiselecclmpdiffmass, mass, axis=1) # [1/cm^3/s/MeV]
    emiselecclmpenel = trapz(emiselecclmp, enel, axis=0) # [1/cm^3/s]
    diffemiselecclmpdiffmassenel = trapz(diffemiselecclmpdiffmass, enel, axis=0) # [1/cm^3/s/Msun]
    ## total
    emiselec = emiselecsmth + emiselecclmp
    emiselecenel = emiselecsmthenel + emiselecclmpenel
    

    # dark matter velocity variance
    diffrvelvariiavgdiffmass = diffnhaldiffmass * rvelvarihavg * 4. * pi * amax(rsph, axis=2)**3 / 3. # [1/Msun]
    rvelvariiavgclmp = trapz(diffrvelvariiavgdiffmass, mass, axis=0) # [1]
    tempiavgsmth = tempcmbrnunc * (1. + reds) / (1. + 1e9  / (1. + reds) / (1. + ((1. + reds) / 1e9)**2.5))
    rvelvariiavgsmth = 3. * boltcons * tempiavgsmth / masspart
    rvelvariiavg = rvelvariiavgsmth + rvelvariiavgclmp


    
    # mean electron number density in the IGM
    ndenelecsavg = zeros((nenel, nreds, nmpol))
    for c in range(nreds):
        for m in range(nmpol):
            for i in range(nenel-1):
                ndenelecsavg[i, c, m] = trapz(emiselec[i:nenel, c, m], enel[i:nenel], axis=0) / edotegbl[i, c] # [1/cm^3/MeV]
    ndenelecsavgenel = trapz(ndenelecsavg, enel, axis=0) # [1/cm^3]

    # photon emissivity
    emisphot = trapz(specinvc[:, :, :, None] * ndenelecsavg[None, :, :, :], enel, axis=1) # [1/cm^3/s/MeV]
    emisphotrydb = interp1d(enph, emisphot, axis=0)(enerrydb) # [1/cm^3/s/MeV]
    emisphotenph = trapz(emisphot, enph, axis=0) # [1/cm^3/s]

    # photon flux 
    difffluxphotdiffreds = zeros((nenph, nreds, nreds, nmpol))
    fluxphot = zeros((nenph, nreds, nmpol))
    for a in range(nenph):
        for c in range(nreds):

            enphshft = enph[a] * (1. + reds[c]) / (1. + reds)
            jredsshft = where((minmenph < enphshft) & (enphshft < maxmenph))[0]
            
            if jredsshft.size == 0:
                continue
                
            for m in range(nmpol):
                
                emisphotintp = interp1d(enph, emisphot[:, c, m])(enphshft[jredsshft]) # [1/cm^3/s/MeV]
                difffluxphotdiffreds[a, c, jredsshft, m] = velolght * timehubb[jredsshft] / 4. / pi                     * emisphotintp * exp(-odep[a, c, jredsshft]) *                     (1. + reds[c])**3 / (1. + reds[jredsshft])**4 # [1/cm^2/s/sr/MeV]

                fluxphot[a, c, m] = trapz(difffluxphotdiffreds[a, c, :, m], reds) # [1/cm^2/s/sr/MeV]  
                
    fluxphotrydb = interp1d(enph, fluxphot, axis=0)(enerrydb) # [1/cm^2/s/sr/MeV]
    fluxphotenph = trapz(fluxphot, enph, axis=0) # [1/cm^2/s/sr]


    # photo-ionization rate
    diffionrdiffenph = 4. * pi * csecionz[:, None, None] * fluxphot # [1/s/MeV]
    ionr = trapz(diffionrdiffenph, enph, axis=0) # [1/s]
    
    return ionr, fluxphot
        


# In[82]:

def plot_halo():

    listlinestyl = [':', '--', '-']
    listcolr = ['b', 'g', 'r']
    
    # concentration parameter
    fig, ax = plt.subplots()
    for m in range(nconcmodl):
        for c in range(njredslate):
            ax.loglog(mass, conccatl[:, jredslate[c], m], ls=listlinestyl[m], color=listcolr[c])
    ax.set_xlabel('$M [M_\odot]$')
    ax.set_ylabel('$c(M,z)$')
    ax.set_title('Mean halo concentration', fontsize=18)
    ax.set_ylim([1., None])

    listline = []
    listline.append(mpl.lines.Line2D([], [], ls=listlinestyl[0], color='black'))
    listline.append(mpl.lines.Line2D([], [], ls=listlinestyl[1], color='black'))
    listline.append(mpl.lines.Line2D([], [], ls=listlinestyl[2], color='black'))
    listline.append(mpl.lines.Line2D([], [], color=listcolr[0]))
    listline.append(mpl.lines.Line2D([], [], color=listcolr[1]))
    listline.append(mpl.lines.Line2D([], [], color=listcolr[2]))
    
    listlabl = []
    listlabl.append('Sanchez-Conde & Prada 2014')
    listlabl.append('Duffy et al. 2008')
    listlabl.append('Comerford & Natarajan 2007')
    listlabl.append(strgredslate[0])
    listlabl.append(strgredslate[1])
    listlabl.append(strgredslate[2])
    
    ax.legend(listline, listlabl)

    plt.savefig(plotpath + 'halo.png')
    plt.close()
    

    # substructure
    #fig, ax = plt.subplots()
    #ax.loglog(rsph, subfs)
    #ax.loglog(rsph, subbfrsp)           
    #ax.set_xlabel('$r$ [kpc]')
    #ax.set_ylabel('$f_s(r)$')


    # DM energy density and relative velocity variance in halos
    temp = meshgrid(jmass, jredslate, indexing='ij')
    labl = strgmass + strgredslate

    fig, ax = plt.subplots()
    fig.suptitle('DM energy density in halos', fontsize=18)                
    plot_matr(ax, rsph[temp[0], temp[1], :], edendmathalo[temp[0], temp[1], :],               labl=labl, loc=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$r$ [kpc]')
    ax.set_ylabel(r'$\rho(r,M,z)$ [MeV/cm$^3$]')
    plt.savefig(plotpath + 'edendmathalo.png')
    plt.close()



    fig, ax = plt.subplots()
    fig.suptitle('DM relative velocity variance in halos', fontsize=18)
    plot_matr(ax, rsph[temp[0], temp[1], :], rvelvarihalo[temp[0], temp[1], :], labl=labl, loc=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$r$ [kpc]')
    ax.set_ylabel('$v^2 (M,z,r)$')
    plt.savefig(plotpath + 'rvelvarihalo.png')
    plt.close()


    # virial and scale radii
    fig, ax = plt.subplots()
    for c in range(njredslate):
        ax.loglog(mass, rsphviri[:, jredslate[c]], label=strgredslate[c])
    ax.set_xlabel('$M [M_\odot]$')
    ax.set_ylabel('$r_{vir}(M,z)$ [kpc]')
    ax.set_title('Halo Virial Radius')
    ax.legend(loc=2)
    plt.savefig(plotpath + 'rsphviri.png')
    plt.close()



    # e^-/e^+ emissivity in a halo
    fig, axgrd = plt.subplots(njmass, nmpol, sharey='row', sharex='all', figsize=(14,21))
    fig.suptitle('$e^-/e^+$ emissivity in a halo', fontsize=18)
    for d, axrow in enumerate(axgrd):
        for p, ax in enumerate(axrow):
            for c in range(njredslate):
                ax.loglog(rsph[jmass[d], jredslate[c], :], emiselechaloenel[jmass[d], jredslate[c], :, p],                           label=strgredslate[c])
            ax.set_xlabel('$r$ [kpc]')
            if d == 0:
                ax.set_title(strgmpol[p])
            ax.text(0.23, 0.15, strgmass[d], ha='center', va='center', transform = ax.transAxes)
            if d == njmass / 2:
                ax.set_ylabel(r'$dN_e/dtdV (M,z,r)$ [1/s/cm$^3$]')
    ax.legend()
    plt.savefig(plotpath + 'emiselechaloenel.png')
    plt.close()



    # e^-/e^+ differential luminosity of a halo
    fig, axgrd = plt.subplots(njmass, nmpol, sharey='row', sharex='all', figsize=(14,21))
    fig.suptitle('$e^-/e^+$ differential luminosity in a halo', fontsize=18)
    for d, axrow in enumerate(axgrd):
        for p, ax in enumerate(axrow):
            for c in range(njredslate):
                ax.loglog(rsph[jmass[d], jredslate[c],:], diffemiselechalodiffrsphenel[jmass[d], jredslate[c],:,p],                           label=strgredslate[c])
            ax.set_xlabel('$r$ [kpc]')
            if d == 0:
                ax.set_title(strgmpol[p])
            ax.text(0.23, 0.12, strgmass[d], ha='center', va='center', transform = ax.transAxes)
            if d == njmass / 2:
                ax.set_ylabel(r'$dN_e/dtd\log r$ [1/s]')
    ax.legend()
    plt.savefig(plotpath + 'diffemiselechalodiffrsphenel.png')
    plt.close()



    # e^-/e^+ luminosity of a halo
    fig, axrow = plt.subplots(1, nmpol, sharey='row', sharex='all', figsize=(14,7))
    fig.suptitle('$e^-/e^+$ luminosity of a DM halo', fontsize=18)
    for p, ax in enumerate(axrow):
        for c in range(njredslate):
            ax.loglog(mass, lumielechaloenel[:, jredslate[c],p], label=strgredslate[c])
        ax.set_xlabel(r'$M$ $[M_\odot]$')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel('$dN_e/dt (M,z)$ [1/s]')
    ax.legend(loc=2)
    plt.savefig(plotpath + 'lumielechaloenel.png')
    plt.close()


    


    # spatially averaged e^-/e^+ differential emissivity
    fig, axrow = plt.subplots(1, nmpol, sharey='row', figsize=(14,7))
    fig.suptitle('$e^-/e^+$ emissivity due to DM annihilation per decade in halo mass', fontsize=18)
    for p, ax in enumerate(axrow):
        for c in range(njredslate):
            ax.loglog(mass, mass * diffemiselecclmpdiffmassenel[:, jredslate[c],p], label=strgredslate[c])
        ax.set_xlabel(r'$M$ $[M_\odot]$')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel(r'$dN_e/dtdVd\log M$ [1/s/cm$^3$]')
            ax.legend()
        ax.set_ylim([1e-36, 1e-28])
    plt.savefig(plotpath + 'diffemiselecclmpdiffmassenel.png')
    plt.close()


    ylabel = '$dN_e/dVdt (z)$ [1/cm$^3$/s]'
    listvarb = [emiselecenel, emiselecclmpenel, emiselecsmthenel]
    listlabl = ['Total', 'Clumpy', 'Smooth']
    fig, axrow = plt.subplots(1, 2, figsize=(14,7))
    fig.suptitle('Spatially averaged $e^-/e^+$ emissivity', fontsize=18)
    for p, ax in enumerate(axrow):
        ax.set_xlabel('$z$')
        for g in range(3):
            ax.loglog(reds, listvarb[g].squeeze()[:, p], label=listlabl[g])
        if p == 0:
            ax.set_ylabel(ylabel)
        ax.set_title(strgmpol[p])
        ax.set_ylim([1e-36, 1e-24])
        ax.loglog(redssfrd, emiselecsfrd, markersize=10,                   label='Star Formation', color='black', ls='none', marker='o')
        if p == 0:
            ax.legend(loc=1)
    plt.savefig(plotpath + 'emiselecenel.png')
    plt.close()



    ylabel = '$dN_e/dVdt (z)$ [1/cm$^3$/s]'
    listvarb = [emiselecenel, emiselecclmpenel, emiselecsmthenel]
    fig, axrow = plt.subplots(1, 2, figsize=(14,7))
    fig.suptitle('Spatially averaged $e^-/e^+$ comoving emissivity', fontsize=18)
    for p, ax in enumerate(axrow):
        ax.set_xlabel('$z$')
        for g in range(3):
            ax.loglog(reds, listvarb[g].squeeze()[:, p] / (1. + reds)**3, label=listlabl[g])
        if p == 0:
            ax.set_ylabel(ylabel)
        ax.set_title(strgmpol[p])
        ax.set_ylim([1e-36, 1e-24])
        ax.loglog(redssfrd, emiselecsfrd / (1. + redssfrd)**3, markersize=10,                   label='Star Formation', color='black', ls='none', marker='o')
        if p == 0:
            ax.legend(loc=1)
    plt.savefig(plotpath + 'emiselecenelcomo.png')
    plt.close()



    # Mean DM relative velocity variance in a halo
    fig, ax = plt.subplots(figsize=(7, 7))
    for c in range(njredslate):
        ax.loglog(mass, sqrt(rvelvarihavg[:, jredslate[c]]), label=strgredslate[c])
    ax.set_title('Mean DM relative velocity variance in a halo', fontsize=18)
    ax.set_ylabel('$\sigma_v$')
    ax.set_xlabel(r'$M$ $[M_\odot]$')
    ax.legend(loc=2)
    plt.savefig(plotpath + 'rvelvarihavg.png')
    plt.close()


    # Mean DM relative velocity variance
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_title('Spatially averaged RMS relative velocity of DM', fontsize=18)
    ax.set_xlabel('$z$')
    ax.set_ylim([1e-12, 1.])
    ax.loglog(reds, sqrt(rvelvariiavgclmp))
    ax.loglog(reds, sqrt(rvelvariiavgsmth))
    ax.loglog(reds, sqrt(rvelvariiavg))
    ax.set_ylabel('$\sigma_v$')
    plt.savefig(plotpath + 'rvelvariiavg.png')
    plt.close()



    
def plot_elec_flux():
    
    # e^-/e^+ number density in the IGM
    fig, axrow = plt.subplots(1, nmpol, sharey='all', figsize=(14,7))
    fig.suptitle('Differential $e^-/e^+$ number density in the IGM', fontsize=18)
    for p, ax in enumerate(axrow):
        for c in range(njredslate):
            ax.loglog(enel[0:nenel-1], enel[0:nenel-1] * ndenelecsavg[0:nenel-1,jredslate[c],p],                       label=strgredslate[c])
        ax.set_xlabel(r'$E_e$ [GeV]')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel(r'$dN_e/dV/d\log E$ [1/cm$^3$]')
    ax.legend(loc=2)
    plt.savefig(plotpath + 'ndenelecsavg.png')
    plt.close()


    # e^-/e^+ number density in the IGM
    fig, axrow = plt.subplots(1, nmpol, sharey='all', figsize=(14,7))
    fig.suptitle('$e^-/e^+$ number density in the IGM', fontsize=18)
    for p, ax in enumerate(axrow):
        ax.loglog(reds, ndenelecsavgenel[:, p])
        ax.set_xlabel('$z$')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel('$dN_e/dV (z)$ [1/cm$^3$]')
    plt.savefig(plotpath + 'ndenelecsavgenel.png')
    plt.close() 
    
    
        
def plot_invrcomp():

    listlinestyl = ['-', '--', '-.']
    listcolr = ['b', 'g', 'r']
    
    fig, ax = plt.subplots()
    for b in range(njenel):
        for c in range(njredslate):
            ax.loglog(enph * 1e6, enph * sum(specinvccatl[:, jenel[b], jredslate[c], :], 1), ls=listlinestyl[b], c=listcolr[c])
    ax.set_xlabel('$E_\gamma$ [eV]')
    ax.set_ylabel(r'$dN_\gamma/dtd\log E_\gamma$ [1/s]')
    ax.set_ylim([1e-15, 1e-7])
    ax.set_title('Spectrum of IC scattered CMB photons')

    listlabl = []
    listlabl.append('$E_e$ = ' + tdpy_util.mexp(enel[jenel[0]] * 1e-3) + ' GeV')
    listlabl.append('$E_e$ = ' + tdpy_util.mexp(enel[jenel[1]] * 1e-3) + ' GeV')
    listlabl.append('$E_e$ = ' + tdpy_util.mexp(enel[jenel[2]] * 1e-3) + ' GeV')
    listlabl.extend(strgredslate)
    
    listline = []
    for k in arange(3):
        listline.append(plt.Line2D((0,1),(0,0), color='black', ls=listlinestyl[k]))
    for k in range(3):
        listline.append(plt.Line2D((0,1),(0,0), color=listcolr[k]))
    ax.legend(listline, listlabl, loc=9, ncol=2) 
    
    plt.savefig(plotpath + 'specinvccatl.png')
    plt.close() 
    

def plot_igma():
      
    extt = [amin(reds), amax(reds), amin(cden), amax(cden)]
    
    fig, ax = plt.subplots(figsize=(14, 14))
    fig.suptitle('Absorber abundance per decade in column density and redshift at 1 Rydberg')
    zdat = reds[None, :] * cden[:, None] * diffnabsdiffcdendiffreds
    imag = ax.imshow(zdat, extent=extt, cmap='Reds', norm=mpl.colors.LogNorm(), aspect='auto')
    for k in range(cdenbrek.size):
        ax.axhline(cdenbrek[k], color='grey', ls='--')
    for k in range(redsbrek.size):
        ax.axvline(redsbrek[k], color='grey', ls='--')
    ax.set_ylabel('$N_{HI}$ [cm$^{-2}$]')
    ax.set_xlabel('$z$')
    ax.set_title('$d^2N_{abs}/d\log N_{HI}/d\log z$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(imag, ax=ax, fraction=0.04)
    plt.savefig(plotpath + 'diffnabsdiffcdendiffreds.png')
    plt.close() 
    
    fig, ax = plt.subplots(figsize=(14, 14))
    fig.suptitle('Optical depth per decade in column density and redshift at 1 Rydberg')
    zdat = reds[None, :] * cden[:, None] * diffoptddiffcdendiffredsrydb
    zdat[where(zdat > 1.)] = 1.
    imag = ax.imshow(zdat, extent=extt, cmap='Reds', norm=mpl.colors.LogNorm(), aspect='auto')
    for k in range(cdenbrek.size):
        ax.axhline(cdenbrek[k], color='grey', ls='--')
    for k in range(redsbrek.size):
        ax.axvline(redsbrek[k], color='grey', ls='--')
    ax.set_ylabel('$N_{HI}$ [cm$^{-2}$]')
    ax.set_xlabel('$z$')
    ax.set_title(r'$d^2\tau_{abs}/d\log N_{HI}/d\log z$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.colorbar(imag, ax=ax, fraction=0.04)
    plt.savefig(plotpath + 'diffoptddiffcdendiffredsrydb.png')
    plt.close() 
    
    fig, ax = plt.subplots()
    ax.loglog(reds, reds * diffoptddiffredsrydb)
    for k in range(redsbrek.size):
        ax.axvline(redsbrek[k], color='grey', ls='--')
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$d\tau / d\log z$')
    ax.set_ylim([None, 1.])
    ax.set_title('Optical depth per decade in redshift at 1 Rydberg')
    plt.savefig(plotpath + 'diffoptddiffredsrydb.png')
    plt.close() 
    
        
    fig, ax = plt.subplots()
    for c in range(njredslate):
        ax.loglog(reds[jredslate[c]:], odeprydb[jredslate[c], jredslate[c]:], label=strgredslate[c])
    for k in range(redsbrek.size):
        ax.axvline(redsbrek[k], color='grey', ls='--')
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$\tau$')
    ax.set_ylim([None, 1.])
    ax.set_title('Optical depth at 1 Rydberg')
    plt.savefig(plotpath + 'odep.png')
    plt.close() 
    
    
def plot_fluxphot():
    
    
    # differential photon emissivity 
    fig, axrow = plt.subplots(1, 2, sharey='all', figsize=(14,7))
    fig.suptitle('Differential photon emissivity', fontsize=18)
    for i, ax in enumerate(axrow):
        for c in range(njredslate):
            plot = ax.loglog(enph * 1e6, enph * emisphot[:, jredslate[c], i], label=strgredslate[c])
        ax.set_xlabel('$E_\gamma$ [eV]')
        ax.set_ylabel(r'$dN_\gamma/dtdVd\log E_\gamma$ [1/s]')
    ax.legend()
    plt.savefig(plotpath + 'emisphot.png')
    plt.close() 
    
    # photon emissivity 
    fig, axrow = plt.subplots(1, 2, sharey='all', figsize=(14,7))
    fig.suptitle('Photon emissivity', fontsize=18)
    for i, ax in enumerate(axrow):
        plot = ax.loglog(reds, enerrydb * emisphotrydb[:, i])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$dN_\gamma/dtdV$ [1/cm$^3$/s]')
    ax.legend()
    plt.savefig(plotpath + 'emisphotrydb.png')
    plt.close() 
    

    # differential photon flux 
    fig, axrow = plt.subplots(1, nmpol, sharey='all', figsize=(14, 7))
    fig.suptitle('Cosmic UV/X-ray background flux per decade in photon energy')
    for m, ax in enumerate(axrow):
        for c in range(njredslate):
            plot = ax.loglog(enph * 1e6, enph * fluxphot[:, jredslate[c], m], label=strgredslate[c])
        ax.set_xlabel('$E_\gamma$ [eV]')
        ax.set_ylabel(r'$dN_\gamma/dtdAd\log E_\gamma$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(plotpath + 'fluxphot.png')
    plt.close() 
    
    # differential photon flux at 1 Ryd
    fig, axrow = plt.subplots(1, nmpol, sharey='all', figsize=(14, 7))
    fig.suptitle('Photon flux at 1 Rydberg')
    for m, ax in enumerate(axrow):
        plot = ax.loglog(reds, enerrydb * fluxphotrydb[:, m])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$dN_\gamma/dtdAd\log E_\gamma$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(plotpath + 'fluxphotrydb.png')
    plt.close() 
    
    # integrated photon flux 
    fig, axrow = plt.subplots(1, 2, sharey='all', figsize=(14,7))
    fig.suptitle('Integral photon flux')
    for m, ax in enumerate(axrow):
        plot = ax.loglog(reds, fluxphotenph[:, m])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$dN_\gamma/dtdA (z)$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(plotpath + 'fluxphotenph.png')
    plt.close() 
        
def plot_ionr():
    

    fig, axrow = plt.subplots(1, nmpol, sharey='all', figsize=(14,7))
    fig.suptitle('Differential meta-galactic ionization rate per decade in photon energy', fontsize=18)
    for i, ax in enumerate(axrow):
        for c in range(njreds):
            plot = ax.loglog(enph * 1e6, diffionrdiffenph[:,jredslate[c],i], label=strgredslate[c])
        ax.set_xlabel('$E_\gamma [eV]$')
        ax.set_ylabel(r'$d\Gamma/d\log E_\gamma$ [1/s/H]')
        ax.axvline(enerrydb * 1e6, ls='--', color='grey')
    ax.legend()
    plt.savefig(plotpath + 'diffionrdiffenph.png')
    plt.close() 
        
    listname = ['ionrbolt.csv', 'ionrbeck.csv', 'ionrfauc.csv']
    listlabl = ['Bolton & Haehnelt (2007)', 'Becker et al. (2007)', 'Faucher-Giguere (2008)']
    listmrkr = ['o', 'D', 'x']
    
    fig, axrow = plt.subplots(1, nmpol, sharey='all', figsize=(14,7))
    for i, ax in enumerate(axrow):
        plot = ax.loglog(reds, ionr[:,i])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$\Gamma$ [1/s/H]')
        
        for k, name in enumerate(listname):
            path = os.environ["PHOT_IONZ_DATA_PATH"] + '/' + name
            data = loadtxt(path)
            ndata = data.shape[0] / 3
            yerr = zeros((2, ndata))
            xdat = data[0:ndata, 0]
            ydat = data[0:ndata, 1] * 1e-12
            yerr[0, :] = data[ndata:2*ndata, 1]
            yerr[1, :] = data[2*ndata:3*ndata, 1]
            yerr = abs(yerr - ydat)
            ax.errorbar(xdat, ydat, yerr=yerr, ls='', marker=listmrkr[k])

    ax.legend()
    plt.savefig(plotpath + 'ionr.png')
    plt.close() 
        


# In[83]:

def retr_fluxphotexpr():
 
    name = os.environ["PHOT_IONZ_DATA_PATH"] + '/xray_background.dat'
    tabl = loadtxt(name)
    
    enphexpr = tabl[:,0]
    enphexpr *= 1e-6 # [MeV]
    fluxphotexpr = tabl[:, 1]
    fluxphotexpr *= 1e-3 / enphexpr**2 # [1/cm^2/s/sr/MeV]

    fenph = where((min(enphexpr) < enph) & (enph < max(enphexpr)))[0]
    
    fluxphotexpr = interp1d(enphexpr, fluxphotexpr)(enph[fenph])
    fluxphotexprvari = (fluxphotexpr * 1.)**2
    
    return fluxphotexpr, fluxphotexprvari, fenph


# In[84]:

def cdfn_logu(data, minmdata, maxmdata):
    dataunit = log(data / minmdata) / log(maxmdata / minmdata)
    return dataunit


def icdf_logu(dataunit, minmdata, maxmdata):
    data = minmdata * exp(dataunit * log(maxmdata / minmdata))
    return data


# In[85]:

def retr_mocksampvarb():
    
    npara = 4
    
    mocksampvarb = empty(npara)
    mocksampvarb[0] = 1e-26
    mocksampvarb[1] = 1e4
    mocksampvarb[2] = 1e5
    mocksampvarb[3] = 1.
    
    return mocksampvarb


def retr_datapara():
    
    npara = 4
    
    dictpara = dict()
    minmpara = zeros(npara)
    maxmpara = zeros(npara)
    namepara = empty(npara, dtype=object)
    scalpara = empty(npara, dtype=object)
    lablpara = empty(npara, dtype=object)
    unitpara = empty(npara, dtype=object)
    varipara = zeros(npara)

    dictpara['csecvelo'] = 0
    namepara[0] = 'csecvelo'
    minmpara[0] = 3e-32
    maxmpara[0] = 3e-20
    scalpara[0] = 'logt'
    lablpara[0] = '$a$'
    unitpara[0] = '[cm$^3$/s]'
    varipara[0] = 3e-1
    
    dictpara['csecfrac'] = 1
    namepara[1] = 'csecfrac'
    minmpara[1] = 1e-2
    maxmpara[1] = 1e10
    scalpara[1] = 'logt'
    lablpara[1] = '$b/a$'
    unitpara[1] = ''
    varipara[1] = 3e-1
     
    dictpara['masspart'] = 2
    namepara[2] = 'masspart'
    minmpara[2] = 1e4
    maxmpara[2] = 1e6
    scalpara[2] = 'logt'
    lablpara[2] = '$M$'
    unitpara[2] = '[MeV]'
    varipara[2] = 3e-1

    dictpara['dmatslop'] = 3
    namepara[3] = 'dmatslop'
    minmpara[3] = 0.8
    maxmpara[3] = 1.5
    scalpara[3] = 'self'
    lablpara[3] = r'$\gamma$'
    unitpara[3] = ''
    varipara[3] = 3e-1

    datapara = namepara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara
    
    return datapara


# In[86]:

def retr_fluxphotdmatintp(csecvelo, csecfrac, masspart, dmatslop):
    

    fluxphotdmatintp = fluxphotdmat[:, :, 0] * csecvelo / csecvelopivt + fluxphotdmat[:, :, 1] * csecfrac / csecfracpivt
    
    if False:
        
        print 'retr_fluxphotdmatintp'
        print 'csecvelo'
        print csecvelo
        print 'csecvelopivt'
        print csecvelopivt
        print 'csecfrac'
        print csecfrac
        print 'csecfracpivt'
        print csecfracpivt
        print 'fluxphotdmat[:, :, 0]'
        print fluxphotdmat[:, :, 0]
        print 'fluxphotdmat[:, :, 1]'
        print fluxphotdmat[:, :, 1]
        print 'fluxphotdmatintp[:, 0]'
        print fluxphotdmatintp[:, 0]
        
        
        
    return fluxphotdmatintp


# In[87]:

def plot_sfrd():
    
    fig, ax = plt.subplots()
    ax.set_title('Star Formation Rate Density')
    
    ax.plot(redssfrd, sfrd)
    ax.set_yscale('log')
    
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'SFRD [erg/Mpc$^3$/yr]')
    plt.savefig(plotpath + 'sfrd.png')
    plt.close()

    


# In[88]:

def plot_hm12(fluxphotdmat=None, listfluxphotdmat=None):
    fig, ax = plt.subplots()
    ax.set_title('UV/X-ray photon background')
    
    ax.fill_between(enph[fenph] * 1e6, enph[fenph] * (fluxphotexpr + sqrt(fluxphotexprvari)),                     enph[fenph] * (fluxphotexpr - sqrt(fluxphotexprvari)), color='lightblue')
    ax.loglog(enph[fenph] * 1e6, enph[fenph] * fluxphotexpr, label='ROSAT')
    
    listname = ['moretotl', 'moreothr', 'more']
    listlabl = ['Total XRB, Moretti et al. (2009)', 'Unresolved XRB, Worsley et al. (2006)', 'Unresolved XRB, Moretti et al. (2012)']
    listcolr = ['black', 'yellow', 'green']
    for k, name in enumerate(listname):
        path = os.environ["PHOT_IONZ_DATA_PATH"] + '/' + name + '.csv'
        datamore = loadtxt(path)
        enphmore = datamore[:, 0] * 1e-3 # [MeV]
        fluxphotmore = ergs2mgev * (180. / pi)**2 * datamore[:, 1] / enphmore**2 # [1/cm^2/s/sr/MeV]
        #ax.fill_between(enphmore * 1e6, 2. * enphmore * fluxphotmore, 0.5 * enphmore * fluxphotmore, color='lightgreen')
        ax.loglog(enphmore * 1e6, enphmore * fluxphotmore, label=listlabl[k], color=listcolr[k])


    listcolr = ['b', 'g', 'r']
    for c in range(njredslate):
        ax.loglog(enph * 1e6, enph * fluxphothm12[:, jredslate[c]],                          label='Haardt & Madau (2012), ' + strgredslate[c])

        if fluxphotdmat != None:
            ax.loglog(enph * 1e6, enph * fluxphotdmat[:, jredslate[c]], label='DM, ' + strgredslate[c])

        if listfluxphotdmat != None:
            tdpy_util.plot_braz(ax, enph * 1e6, enph[None, :] * listfluxphotdmat[:, :, jredslate[c]],                                 lcol=listcolr[c], alpha=0.5, dcol=listcolr[c], mcol='black')
            
    ax.set_xlabel(r'$E_\gamma$ [eV]')
    ax.set_ylabel(r'$EdN/dE$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(plotpath + 'fluxphothm12.png')
    plt.close()

    


# In[89]:

def init(cnfg):
    
    global makeplot, verbtype
    makeplot = cnfg['makeplot']
    verbtype = cnfg['verbtype']
    
    nswep = cnfg['nswep']
    nburn = cnfg['nburn']
    nthin = cnfg['nthin']

    global plotpath
    plotpath = os.environ["PHOT_IONZ_DATA_PATH"] + '/png/'

    global minmpara, maxmpara, scalpara, namepara, lablpara, unitpara, varipara, dictpara, ipara, npara
    mocksampvarb = retr_mocksampvarb()
    datapara = retr_datapara()
    namepara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara = datapara    
    npara = len(lablpara)
    ipara = arange(npara)
    
    global strgmpol
    strgmpol = ['s-wave', 'p-wave']
    
    # conversion factors
    global cmet2angs, cmet2angs, ergs2mgev, myrs2secd, solm2mgev,     kprc2cmet, magf2mage, evee2rydb
    
    # fundamental constants
    global masselec, velolght, strtcons, hubbcons, plnkcons, boltcons, gravconsredu
    global csecthom, csecionz, enerrydb
    
    # astrophyical constants 
    global massmilk, radisolr
    
    # cosmological constants
    global omegdene, omegcurv, omegdmat, edencritnunc, tempcmbrnunc
    global shtrnorm, shtrwgth, shtrindx
    
    global odenviri, odencoll
     
    # cosmological parameters
    global odenrmsq8mpc, psecindx
    
    global anch
    anch = 'b'
    
    # axes
    global nreds, nradi, nrsph, nzaxi, nenph, nenel, nenpi, ncden, nmass, nwnum
    global reds, radi, rsph, zaxi, enph, enel, enpi, cden, mass, massp, wnum
    global jreds, jredslate, strgredslate, jmass, jenel
    global njreds, njredslate, njmass, njenel
    

    global demccons, demcsigm
    demccons = 35. / 18.
    demcsigm = 6.
    gravlght = 1.19e-34 # [cm^3/MeV/s^2]
    
    # constants
    ## cosmological constants

    omegbmat = 0.049 # baryonic matter abundance today
    omegdmat = 0.26 # dark matter abundance today
    omegradi = 4.8e-5 # radiation abundance today
    omegdene = 0.69 # dark energy abundance today
    odenrmsq8mpc = 0.83 # rms density fluctuation in spheres of radius 8/h Mpc
    psecindx = 0.96 # spectral index of the primordial power spectrum
    hubbcons = 0.704 # reduced Hubble constant
    tempcmbrnunc = 2.725 # CMB temperature today [K]
    shtrwgth = 0.707 # Sheth-Tormen
    shtrnorm = 0.3222 # Sheth-Tormen
    shtrindx = 0.3 # Sheth-Tormen
    
    omegmatt = omegbmat + omegdmat
    
    masselec = 0.511 # electron mass [MeV]
    velolght = 2.998e10 # speed of light [cm/s]
    strtcons = 7.297e-3 # fine structure constant
 

    
    odencoll = 1.686 # linear overdensity at collapse
    odenviri = 18. * pi**2 # overdensity at virialization
    
    radisolr = 8.5 # radial distance from the GC to the Sun [kpc]
    edencritnunc = 5.3e-3 # critical density of the Universe [MeV/cm^3]

    csecthom = 6.65e-25 # [cm^2] Thompson cross section 
    csecionz = 6.3e-18 # [cm^2] neutral Hydrogen photon-ionization cross section at 13.6 eV 
    enerrydb = 13.5984e-6 # Rydberg energy [MeV]

    plnkcons = 4.136e-21 # Planck constant [MeV s]
    boltcons = 8.6173e-11 # Boltzmann constant [MeV/K]

    massmilk = 1e12 # mass of the Milky Way [Msun]

    ergs2mgev = 6.241509e5 # conversion factor from erg to MeV
    myrs2secd = 3.154e13 # million year per second
    solm2mgev = 1.115e60 # Solar mass in MeVc^2
    kprc2cmet = 3.086e21 # conversion factor from kpc to cm
    magffac = 4.966835e-8 # [(MeV/cm^3)/(muG^2/mu0)] 

    
    demccons = 35. / 18. # Dehnen McLaughlin profile index
    sigm = 6. # ratio of mean squared relative velocity to 1-particle velocity variance
    gravconsredu = 1.19e-34 # [cm^3/MeV/s^2]
    cmet2angs = 1e8
    
    
    # axes
    # temp
    nradi = 50
    nrsph = 50
    nzaxi = 50
    nenph = 50
    nenel = 50
    nenpi = 50
    ncden = 50
    nreds = 50
    nmass = 50
    nwnum = 500
    
    global minmenph, maxmenph
    minmenph = 1e-6 # [MeV]
    maxmenph = 1e-1 # [MeV]
    enph = logspace(log10(minmenph), log10(maxmenph), nenph)
    jenph = [0, nenph/2, nenph-1]
    njenph = len(jenph)

    global minmenel, maxmenel
    minmenel = 5e1 # [MeV]
    maxmenel = 1e5 # [MeV]
    enel = logspace(log10(minmenel), log10(maxmenel), nenel)
    denel = enel[1:nenel-1] - enel[0:nenel-2]
    jenel = [0, nenel / 4, nenel / 2]
    njenel = len(jenel)

    

    global minmreds, maxmreds
    minmreds = 1e-1
    maxmreds = 1e2
    reds = logspace(log10(minmreds), log10(maxmreds), nreds)
    
    redsprox = [0., 1e3, 1e6]
    njreds = len(redsprox)
    jreds = []
    for k in range(njreds):
        jreds.append(argmin(abs(reds - redsprox[k])))
        
    redslateprox = [0.1, 2., 6.]
    njredslate = len(redslateprox)
    jredslate = []
    strgredslate = []
    for k in range(njredslate):
        jredslate.append(argmin(abs(reds - redslateprox[k])))
        strgredslate.append('$z = %.3g$' % redslateprox[k])
        
    global minmmass, maxmmass
    minmmass = 1e8 # [Solar Mass]
    maxmmass = 1e16 # [Solar Mass]
    massp = logspace(log10(minmmass), log10(maxmmass), nmass + 1)
    mass = massp[0:-1]
    
    massprox = [1e10, 1e12, 1e15]
    njmass = len(massprox)
    jmass = []
    global strgmass
    strgmass = []
    for d in range(njmass):
        jmass.append(argmin(abs(mass - massprox[d])))
        strgmass.append('$M$ = ' + tdpy_util.mexp(massprox[d]) + r' $M_\odot$')

            
    minmcden = 10**11. # [1/cm^2]
    maxmcden = 10**22. # [1/cm^2]
    cden = logspace(log10(minmcden), log10(maxmcden), ncden)

    minmenpi = 1e-12 # [MeV]
    maxmenpi = 1e-5 # [MeV]
    enpi = logspace(log10(minmenpi), log10(maxmenpi), nenpi)

    # wavenumber axis
    minwnum = 1e-4
    maxwnum = 1e4
    wnum = logspace(log10(minwnum), log10(maxwnum), nwnum)
    


    
    
    global edenbmat, edendmat, edenmatt, edenradi, edendene
    edenbmat = omegbmat * edencritnunc * (1. + reds)**3
    edendmat = omegdmat * edencritnunc * (1. + reds)**3
    edenmatt = omegmatt * edencritnunc * (1. + reds)**3
    edenradi = omegradi * edencritnunc * (1. + reds)**4
    edendene = omegdene * edencritnunc
    
    global edendmatnunc
    edendmatnunc = omegdmat * edencritnunc
    
    global rsphviri
    rsphviri = (3. * mass[:, None] * solm2mgev / 4. / pi / odenviri / edendmat / kprc2cmet**3)**(1. / 3.) # [kpc]
    minmrsph = rsphviri * 1e-3
    maxmrsph = rsphviri * 1e1
    rsph = zeros((nmass, nreds, nrsph))
    for c in range(nreds):
        for d in range(nmass):             
            rsph[d, c, :] = logspace(log10(minmrsph[d,c]), log10(maxmrsph[d,c]), nrsph) # [kpc]

    
    frph = enph / plnkcons # [Hz]
    csecionz = csecionz * (enerrydb / enph)**3
    
    
    global timehubbnunc, funchubb, timehubb
    timehubbnunc = 1. / (100. * hubbcons / (kprc2cmet / 1e2))
    funchubb = sqrt(omegdene + omegdmat * (1. + reds**3) + omegradi * (1. + reds**4))
    timehubb = timehubbnunc / funchubb
    
    # plot annihilation spectrum
    if makeplot:
        anchlabl = ['$e^-e^+$', r'$\mu\bar{\mu}$', r'$\tau^-\tau^+$', r'$b\bar{b}$']
        masspartlabl = ['10 GeV', '100 GeV', '1 TeV']
        anchlist = ['e', 'mu', 'tau', 'b']
        colorlist = ['b', 'g', 'r', 'm']
        linelist = ['-', '--', '-.']
        masspart = array([1e4, 1e5, 1e6])
        nmasspart = masspart.size

        global multp4dm, enelscalp4dm, massp4dm
        nanch = len(anchlist)
        fig, ax = plt.subplots()
        for a, anch in enumerate(anchlist):
            multp4dm, enelscalp4dm, massp4dm = tdpy_util.retr_p4dm_spec(anch)
            for k in range(nmasspart):
                jmassp4dm = argmin(abs(massp4dm - masspart[k]))
                ax.loglog(enelscalp4dm * masspart[k] * 1e-3, multp4dm[:, jmassp4dm], color=colorlist[a],                           ls=linelist[k], label=anchlabl[a] + ', ' + masspartlabl[k])
        ax.set_xlim([0.05, 1e3]) 
        ax.set_ylim([1e-1, 5e2])
        ax.set_xlabel('$E$ [GeV]')
        ax.set_ylabel(r'$dN/d\log E$')
        ax.set_title('Prompt product spectrum per annihilation')
        ax.legend(ncol=4, loc=2)
        plt.savefig(plotpath + 'multp4dm.png')
        plt.close()


    
    multp4dm, enelscalp4dm, massp4dm  = tdpy_util.retr_p4dm_spec(anch)

    #if makeplot:
        #edot_plot
        #galprop_plot

    # Haardt Madau 2012 quasar and galaxy background model
    global fluxphothm12
    fluxphothm12 = retr_fluxphothm12()
    
    # experimental background
    global fluxphotexpr, fluxphotexprvari, fenph
    fluxphotexpr, fluxphotexprvari, fenph = retr_fluxphotexpr()
    
    if makeplot:
        plot_hm12()
        
        
    
    # CMB energy density
    global edencmbr
    edencmbr = eden_cmbr() # [1/cm^3/MeV]
    edencmbrenpi = trapz(enpi[:, None] * edencmbr, enpi, axis=0) # [MeV/cm^3]
    
    # EGBL energy density
    global edenegbl
    edenegbl = eden_egbl() # [1/cm^3/MeV]
    edenegblenpi = trapz(enpi[:, None] * edenegbl, enpi, axis=0) # [MeV/cm^3]
    
    # Energy loss rate on EGBL
    global edotegbl
    # temp
    edotegbl = 4. * csecthom * velolght / 3. / masselec**2 * edencmbrenpi[None, :] * enel[:,None]**2 # [MeV/s]
 
    if makeplot:
        fig, ax = plt.subplots()
        for c in range(njredslate):
            plot = ax.loglog(enpi * 1e6, enpi * edencmbr[:, jredslate[c]], label=strgredslate[c])
            plot_ = ax.loglog(enpi * 1e6, enpi * edenegbl[:, jredslate[c]], '--', color=plot[0].get_color())
        ax.set_xlabel('$E_\gamma$ [eV]')
        ax.set_ylabel(r'$dN_\gamma/dVd\log E_\gamma$ [1/cm$^3$]')
        ax.set_title('Extragalactic Background Radiation')
        ax.set_ylim([1e-1, 1e6])
        leg = ax.legend(loc=2)
        plt.savefig(plotpath + 'edencmbr.png')
        plt.close()

    # halo mass function
    global diffnhaldiffmass, peakhght
    diffnhaldiffmass, peakhght = retr_hmfn()
    

    global nmpol    
    nmpol = 2
    
    global propmodl, concmodl, subsmodl, igmamodl, conc

    propmodl = 'effi'
    concmodl = 'duff'
    subsmodl = 'none'
    igmamodl = 'clum'
    

    # halo concentration model
    global conccatl, conc, nconcmodl
    nconcmodl = 3
    conccatl = zeros((nmass, nreds, nconcmodl))

    ## Sanchez-Conde & Prada 2014
    cona = [37.5153, -1.5093, 1.636e-2, 3.66e-4, -2.89237e-5, 5.32e-7]
    for k in range(6):
        conccatl[:, :, 0] += cona[k] * log(mass[:, None] / hubbcons)**k * funchubb[None, :]**(-2./3.)
    if concmodl == 'sanc':
        conc = conccatl[:, :, 0]

    ## Duffy et al. 2007
    cmind = -0.081
    crind = -0.71
    cmpiv = 2e12 / hubbcons
    cmnor = 7.8
    conccatl[:,:,1] = cmnor * outer((mass / cmpiv)**cmind, (1. + reds)**crind)
    if concmodl == 'duff':
        conc = conccatl[:,:,1]

    ## Comerford & Natarajan 2007
    cmind = -0.15
    crind = -1.
    cmpiv = 1.3e13 / hubbcons
    cmnor = 14.5
    conccatl[:, :, 2] = cmnor * (mass[:, None] / cmpiv)**cmind * (1. + reds[None, :])**crind
    if concmodl == 'come':
        conc = conccatl[:,:,2]

    # substructure model
    global subs
    subs = ones((nmass, nreds, nrsph))
    if subsmodl == 'none':
        subs[:] = 1.
    
    # cosmic ray escape model
    uisrf = 6e-7 # [MeV/cm^3]
    magfd = 5. # [muG]
    gasdn = 1. # [1/cm^3]

    # cosmic raw propagation model
    global fesc
    fesc = zeros((nenel, nmass, nreds, nrsph))
    if propmodl == 'effi':
        fesc[:] = 1.
    else:
        if propmodl == 'minm':
            magfd *= 0.2 # [muG]
            gasdn *= 0.2 # [1/cm^3]
            uisrf *= 0.2 # [MeV/cm^3]
        if propmodl == 'maxm':
            magfd *= 5. # [muG]
            gasdn *= 5. # [1/cm^3]
            uisrf *= 5. # [MeV/cm^3]

        e0 = 1. # [GeV]
        if propmodl == 'minm':
            k0 = 5.074e-17 # [kpc^2/s]
            delta = 0.85
            lz = 1. # [kpc]
        if propmodl == 'medi':
            k0 = 3.551e-16 # [kpc^2/s]
            delta = 0.7
            lz = 4. # [kpc]
        if propmodl == 'maxm':
            k0 = 2.426e-15 # [kpc^2/s]
            delta = 0.46
            lz = 8. # [kpc]
        
        fesc = retr_fesc()
        
    # IGM effective optical depth model
    global cdenbrek, redsbrek
    cdenbrek = array([10**15., 10**17.5, 10**19., 10**20.3])
    redsbrek = array([1.56, 5.5])
    global odep, odeprydb
    odep = zeros((nenph, nreds, nreds))
    odeprydb = zeros((nreds, nreds))
    if igmamodl == 'clum':
        
        global diffnabsdiffcdendiffreds, diffoptddiffcdendiffredsrydb, diffoptddiffredsrydb
        diffnabsdiffcdendiffreds = zeros((ncden, nreds))

        # Lyman - alpha forest
        jcden = where(cden < 10**15.)[0]

        jreds = where((reds > 1.56) & (reds < 5.5))[0]
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**7.079 * cden[jcden,None]**(-1.5) * (1. + reds[None,jreds])**3

        jreds = where(reds < 1.56)[0]
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**8.238 * cden[jcden,None]**(-1.5) * (1. + reds[None,jreds])**0.16

        jreds = where(reds > 5.5)[0]
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1, amin(jreds):amax(jreds)+2] = 10**1.470 * cden[jcden,None]**(-1.5) * (1. + reds[None,jreds])**9.9

        jcden = where((cden > 10**15.) & (cden < 10**17.5))

        jreds = where((reds > 1.56) & (reds < 5.5))
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**14.58 * cden[jcden,None]**(-2.) * (1. + reds[None,jreds])**3.

        jreds = where(reds < 1.56)
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**15.74 * cden[jcden,None]**(-2.) * (1. + reds[None,jreds])**0.16

        jreds = where(reds > 5.5)
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**8.97 * cden[jcden,None]**(-2.) * (1. + reds[None,jreds])**9.9

        # Super Lyman limit systems
        jcden = where((cden > 10**19.) & (cden < 10**20.3))

        jreds = where(reds > 1.56)
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**(-0.347) * cden[jcden,None]**(-1.05) * (1. + reds[None,jreds])**1.27

        jreds = where(reds < 1.56)
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**(0.107) * cden[jcden,None]**(-1.05) * (1. + reds[None,jreds])**0.16

        # Damped Lyman alpha systems
        jcden = where(cden > 10**20.3)

        jreds = where(reds > 1.56)
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**18.94 * cden[jcden,None]**(-2.) * (1. + reds[None,jreds])**1.27

        jreds = where(reds < 1.56)
        diffnabsdiffcdendiffreds[amin(jcden):amax(jcden)+1,amin(jreds):amax(jreds)+1] = 10**19.393 * cden[jcden,None]**(-2.) * (1. + reds[None,jreds])**0.16

        diffoptddiffcdendiffreds = diffnabsdiffcdendiffreds[:, None, :] *             (1. - exp(-cden[:, None, None] * csecionz[None, :, None]))
            
        diffoptddiffreds = trapz(diffoptddiffcdendiffreds, cden, axis=0)
        
        diffoptddiffredsrydb = interp1d(enph, diffoptddiffreds, axis=0)(enerrydb)
        diffoptddiffcdendiffredsrydb = interp1d(enph, diffoptddiffcdendiffreds, axis=1)(enerrydb)
        
        odep = zeros((nenph, nreds, nreds))
        for c in range(nreds):
            for g in range(c + 1, nreds):
                odep[:, c, g] = trapz(diffoptddiffreds[:, c:g], reds[c:g], axis=1)

        odeprydb = interp1d(enph, odep, axis=0)(enerrydb)
        
        if makeplot:
            plot_igma()
                
    # compute the differential ICS power as a function of final photon and electron energies
    global specinvc, specinvccatlenph, specinvccatl
    specinvccatl = zeros((nenph, nenel, nreds, 2))
    nqics = 100
    maxmqics = 1.
    for b in range(nenel):
        for a in range(nenph):
            elecrelefact = enel[b] / masselec
            eps = enph[a] / enel[b]
            if eps >= 1.:
                continue
            for c in range(nreds):
                minmqics = 1. / 4. / elecrelefact**2
                qics = logspace(log10(minmqics), log10(maxmqics), nqics)
                enpitemp = masselec**2 / 4. / enel[b] * eps / (1. - eps) / qics
                qfac = 2. * qics * log(qics) + qics + 1. - 2. * qics**2 + eps**2 * (1. - qics) / 2. / (1. - eps)
                specinvctemp = 3. * csecthom * velolght / 4. / elecrelefact**2 *                     (enph[a] - enpitemp) * qfac / qics / enph[a] # [cm^3/s]
                jqics = where((enpitemp < max(enpi)) & (enpitemp > min(enpi)))[0]
                specinvccatl[a, b, c, 0] = trapz(specinvctemp[jqics] * interp1d(enpi, edencmbr[:,c])(enpitemp[jqics]), qics[jqics]) # [1/s/MeV] 
                specinvccatl[a, b, c, 1] = trapz(specinvctemp[jqics] * interp1d(enpi, edenegbl[:,c])(enpitemp[jqics]), qics[jqics]) # [1/s/MeV]
        
    # temp
    specinvc = specinvccatl[:, :, :, 0]
    
    specinvccatlenph = trapz(specinvccatl, enph, axis=0) # [1/s]
         
    
    # star formation rate density
    global sfrd, redssfrd, emiselecsfrd
    path = os.environ["PHOT_IONZ_DATA_PATH"] + '/sfrd.csv'
    sfrddata = loadtxt(path)
    redssfrd = sfrddata[:, 0]
    sfrd = sfrddata[:, 1]
    emiselecsfrd = 1e51 * ergs2mgev / (1e3 * kprc2cmet)**3 / (1e-6 * myrs2secd) * 0.15 * sfrd * 1e-2 * (1. + redssfrd)**3
    

         
    global csecvelopivt, csecfracpivt, masspartpivt, dmatsloppivt
    csecvelopivt = 3e-26
    csecfracpivt = 1e4
    masspartpivt = 1e4
    dmatsloppivt = 1.
    global ionrdmat, fluxphotdmat
    
    path = os.environ["PHOT_IONZ_DATA_PATH"] + '/data.fits'
    # temp
    if False and os.path.isfile(path):
        ionrdmat = pf.getdata(path, 0)
        fluxphotdmat = pf.getdata(path, 1)
    else:
        ionrdmat, fluxphotdmat = retr_fluxphotdmat(csecvelopivt, csecfracpivt, masspartpivt, dmatsloppivt)
        pf.append(path, ionrdmat)
        pf.append(path, fluxphotdmat)

        if makeplot:
            plot_halo()
            plot_elec_flux()
            plot_invrcomp()
            plot_fluxphot()
            plot_ionr()
            plot_sfrd()


    optiprop = False
    nsampcalc = 1

    thissampvarb = zeros(npara)
    thissampvarb[0] = csecvelopivt
    thissampvarb[1] = csecfracpivt
    thissampvarb[2] = masspartpivt
    thissampvarb[3] = dmatsloppivt
    thissamp = tdpy_util.cdfn_samp(thissampvarb, datapara)

    nsamp = (nswep - nburn) / nthin
    

    sampbund = tdpy_util.mcmc(nswep, retr_llik, datapara, optiprop=optiprop,                               thissamp=thissamp, nburn=nburn,                               nthin=nthin, verbtype=verbtype,                               nsampcalc=nsampcalc, plotpath=plotpath)
    
    listsampvarb = sampbund[0]
    listsamp = sampbund[1]
    listsampcalc = sampbund[2]
    listllik = sampbund[3]
    listaccp = sampbund[4]
    listjsampvari = sampbund[5]
    
    
    listfluxphotdmat = empty((nsamp, nenph, nreds))
    for k in range(nsamp):
        listfluxphotdmat[k, :, :] = listsampcalc[0][k]
          
            
    fig, ax = plt.subplots()
    ax.set_title('UV/X-ray photon background')
    
    ax.fill_between(enph[fenph] * 1e6, enph[fenph] * (fluxphotexpr + sqrt(fluxphotexprvari)),                     enph[fenph] * (fluxphotexpr - sqrt(fluxphotexprvari)), color='lightblue')
    
    ax.loglog(enph[fenph] * 1e6, enph[fenph] * fluxphotexpr, label='ROSAT')
    
    ax.loglog(enph * 1e6, enph * fluxphothm12[:, 0], label='Haardt & Madau (2012)')
    tdpy_util.plot_braz(ax, enph * 1e6, enph[None, :] * listfluxphotdmat[:, :, 0], alpha=0.5, mcol='black')

    ax.set_xlabel(r'$E_\gamma$ [eV]')
    ax.set_ylabel(r'$EdN/dE$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(plotpath + 'fluxphotpost.png')
    plt.close()
    
    strgpara = lablpara + ' ' + unitpara
    path = os.environ["PHOT_IONZ_DATA_PATH"] + '/png/mcmc'
    tdpy_util.plot_mcmc(listsampvarb, strgpara, scalpara=scalpara, path=path, quan=True)

    for k in range(npara):
        path = plotpath + namepara[k] + '.png'
        tdpy_util.plot_trac(listsampvarb[:, k], strgpara[k], scalpara=scalpara[k], path=path, quan=True)


# In[90]:

def retr_cnfg(               datatype='inpt',               datalabl='PIXIE',               nswep=100,               nburn=None,               nthin=None,               nfreqexpr=100,               minmfreqexpr=1e9,               maxmfreqexpr=1e13,               plotperd=10000,               verbtype=1,               makeplot=False,               ):
        
    cnfg = dict()
    
    # data type and label
    cnfg['datatype'] = datatype
    cnfg['datalabl'] = datalabl

    # sampler setup
    cnfg['nswep'] = nswep
    cnfg['nburn'] = nburn
    cnfg['nthin'] = nthin

    # frequency axis
    cnfg['nfreqexpr'] = nfreqexpr
    cnfg['minmfreqexpr'] = minmfreqexpr
    cnfg['maxmfreqexpr'] = maxmfreqexpr

    #
    cnfg['plotperd'] = plotperd
    cnfg['verbtype'] = verbtype
    cnfg['makeplot'] = makeplot
    
    return cnfg



# In[91]:

def cnfg_nomi():
    
    cnfg = retr_cnfg(                      nswep=100000,                      verbtype=0,                      nburn=0,                      nthin=1,                      makeplot=True                     )
    
    init(cnfg)
    


# In[92]:

if __name__ == '__main__':
    
    cnfg_nomi()


# In[ ]:



