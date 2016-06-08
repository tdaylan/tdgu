import os, time

from numpy import *
from numpy.random import *
from numpy.random import choice

import scipy as sp
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline

import pyfits as pf

import tdpy.util
import tdpy.mcmc

# plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)

# determine the thermal state of the IGM due to photoheating
## blazar heating rate
    
    #plheatbl = 10**(0.0315 * (1. + red)**3 - 0.512 * (1. + red)**2 + 2.27 * (1. + red) - 2.38) / 3.154e22 # [MeV/s]
    #plheatbl[where(reds gt 5.7)] = 0.
    
## photoheating rate
    
    #  phheatdm = zeros(numbreds)
    #  phheatqg = zeros(numbreds)
    #  phheatqghm = zeros(numbreds)
    #  phheatqs = zeros(numbreds)
    #  for c=0, nred-1 do begin
    #    phheatdm[c] = trapz(enph, pucr_phheatint(phdmflux[:,c]), xr=[ryd, enph[numbenph-1]]) # [MeV/s]
    #    phheatqg[c] = trapz(enph, pucr_phheatint(phqgflux[:,c]), xr=[ryd, enph[numbenph-1]]) # [MeV/s]
    #    phheatqghm[c] = trapz(enph, pucr_phheatint(phqgfluxhm[:,c]), xr=[ryd, enph[numbenph-1]]) # [MeV/s]
    #    phheatqs[c] = trapz(enph, pucr_phheatint(phqsflux[:,c]), xr=[ryd, enph[numbenph-1]]) # [MeV/s]
    #  endfor
    

    #if makeplot:
    #    plot_heat
    
    
    # density - temperature relation
    
    #  npatch = 100
    #  initreds = 20.
    #  mindif = min(abs(reds - initred), jredtemp)
    #  temp = cmbt
    #
    #  puc_lambda
    #  heat = phheatdm
    #  temppatchdm = zeros((npatch, numbreds))
    # dpatchdensdzarr = zeros((npatch, numbreds))
    #  patchdensarr = zeros((npatch, numbreds))
    #  for i=0, npatch-1 do begin
    #    patchdensarr[i,:] = 1. / ((1. + lambda[i,0] * grwf) * (1. + lambda[i,1] * grwf) * (1. + lambda[i,2] * grwf)) - 1.
    #    dpatchdensdzarr[i,:] = deriv(red, patchdensarr[i,:])
    #    patchdens = patchdensarr[i,:]
    #    dpatchdensdz = dpatchdensdzarr[i,:]
    #    temppatchdm[i,:] = puc_temp_solve
    #  endfor    

    #if makeplot:
    #    plot_ptch 
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
    #    plot_temp 


def retr_edenegbl():
    
    path = os.environ["PHOT_IONZ_DATA_PATH"] + '/egbl.csv'
    egbldata = loadtxt(path)
    wlenegbl = egbldata[:, 0] * 1e-4 # [cm]
    freqegbl = flipud(velolght / wlenegbl) # [Hz]
    enpiegbl = plnkcons * freqegbl # [MeV]
    fluxegbl = egbldata[:, 1] # [W/m^2/sr]
    edenegbltemp = flipud(fluxegbl) * 2. * pi / velolght / 1e4 * 1.6e13 / enpiegbl**2 # [1/cm^3/MeV]
    edenegbl = zeros((numbenpi, numbreds))
    for c in indxreds:
        enpitemp = enpi / (1. + reds[c]) # [MeV]
        indxenpitemp = where((enpitemp < max(enpiegbl)) & (enpitemp > min(enpiegbl)))
        edenegbl[indxenpitemp,c] = interp1d(enpiegbl, edenegbltemp)(enpitemp[indxenpitemp]) * (1. + reds[c])**2 # [1/cm^3/MeV]

    return edenegbl
        

def retr_edencmbr():
    
    edencmbr = 8. * pi * enpi[:, None]**2 / (velolght * plnkcons)**3 / (exp(enpi[:,None] / tempcmbrnunc / boltcons / (1. + reds[None, :])) - 1.) # [1/cm^3/MeV]
    
    return edencmbr


def plot_edot():

    sec2gyr = 3.171e-17 # [s/Gyrs]
    
    ndengas = [0.01, 1., 100.] # [1/cm^3]
    edengas = [1e-8, 1e-6, 1e-4] # [MeV/cm^3]
    intsmag = [1e-10, 1e-8, 1e-6] # [muG]
    
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
    
    plt.savefig(pathplot + 'edotelec.png')
    plt.close()


def retr_psec(thiswnum):

    q = thiswnum / 0.15
    cq = 14.4 + 325 / (1. + 60.5 * q**1.11)
    lq = log(exp(1.) + 1.84 * q)
    tranfunc = lq / (lq + cq * q**2)
    psecprim = 2e-9 * thiswnum**psecindx
    psec = psecprim * tranfunc**2
    
    return psec, tranfunc


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
        
    sampcalc = [fluxphotdmatintp[:, 0]]
    
    return lpos, sampcalc


def retr_hmfn():
      
    # growth factor
    grwf = zeros(numbreds) # growth factor
    for c in indxreds:
        diffgrwfdiffreds = funchubb[c] * (1. + reds) / funchubb**3
        grwf[c] = trapz(diffgrwfdiffreds[c:], reds[c:])
    grwf /= grwf[0]
    
    # radius, wavelength and wavenumber corresponding to the halo mass
    rsphhalo = (3. * massprim * solm2mgev / 4. / pi / omegdmat / edencritnunc / odenviri)**(1./3.) / kprc2cmet / 1e3 # [Mpc]
    wlenhalo = 4. * rsphhalo
    wnumhalo = 2. * pi / wlenhalo

    # power spectrum of density fluctuations
    psec, tranfunc = retr_psec(wnum)

    # RMS density fluctuations
    fluc = zeros(numbmass + 1)
    diffflucdiffwnum = zeros((numbwnum, numbmass + 1))
    funcwndw = zeros((numbwnum, numbmass + 1))
    for d in range(numbmass + 1):
        wang = wnum * wlenhalo[d]
        funcwndw[:, d] = 3. * (sin(wang) - wang * cos(wang)) / wang**3
        diffflucdiffwnum[:, d] = wnum**2 * psec * funcwndw[:, d]**2 / 2. / pi**2
        fluc[d] = sqrt(trapz(diffflucdiffwnum[:, d], wnum, axis=0))
    # temp
    fluc *= 0.55 / interp1d(massprim, fluc)(1e15)
    #fluc *= odenrmsq8mpc / interp1d(rsphhalo, fluc)(8. / hubbcons)
    fluc = fluc[:, None] * grwf[None, :]

    # halo mass function
    difflogtflucdiffmass = -diff(log(fluc), axis=0) / diff(massprim)[:, None]
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
        plt.savefig(pathplot + 'grwf.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.loglog(massprim, rsphhalo)
        ax.set_xlabel(r'$M [M_\odot]$')
        ax.set_ylabel('$r_H$ [Mpc]')
        ax = ax.twinx()
        ax.loglog(massprim, wnumhalo, ls='--')
        ax.set_ylabel('$k_H$ [Mpc$^{-1}$]')
        ax.legend()
        plt.savefig(pathplot + 'rsphhalo.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('P(k) [Mpc$^3$]')
        ax.loglog(wnum, psec)
        ax.set_title('Primordial Matter Power Spectrum')
        plt.savefig(pathplot + 'psec.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_title('Transfer function')
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('$T(k)$')
        ax.loglog(wnum, tranfunc)
        plt.savefig(pathplot + 'tranfunc.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_title('Window function')
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('$W(kR)$')
        for d in range(numbmassplot):
            ax.loglog(wnum, funcwndw[:, indxmassplot[d]]**2, label=strgmass[d])
        ax.legend(loc=3)
        plt.savefig(pathplot + 'funcwndw.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.set_title('Contribution of spatial scales to the RMS density fluctuations')
        ax.set_xlabel('$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('$d\sigma^2/dk(M)$')
        for d in range(numbmassplot):
            ax.loglog(wnum, diffflucdiffwnum[:, indxmassplot[d]], label=strgmass[d])
        ax.legend(loc=3)
        plt.savefig(pathplot + 'diffflucdiffwnum.png')
        plt.close()
        
        fig, ax = plt.subplots()
        for c in range(numbredsplotlate):
            ax.loglog(mass, difflogtflucdiffmass[:, indxredsplotlate[c]], label=strgredslate[c])
        ax.set_xlabel(r'$M [M_\odot]$')
        ax.set_ylabel(r'd$\log\sigma$/d$M$')
        ax.legend()
        plt.savefig(pathplot + 'difflogtflucdiffmass.png')
        plt.close()
        
        fig, ax = plt.subplots()
        for c in range(numbredsplotlate):
            ax.loglog(massprim, fluc[:, indxredsplotlate[c]], label=strgredslate[c])
        ax.set_xlabel(r'$M [M_\odot]$')
        ax.set_ylabel(r'$\sigma$')
        ax.axhline(odencoll, ls='--', label='Critical linear overdensity at collapse')
        ax.legend()
        plt.savefig(pathplot + 'fluc.png')
        plt.close()
        
        fig, ax = plt.subplots()
        ax.loglog(fluc[:-1, 0], funcfluc[:, 0])
        ax.set_xlabel(r'$\sigma$')
        ax.set_ylabel('$f(\sigma)$')
        plt.savefig(pathplot + 'funcfluc.png')
        plt.close()
        
        datamsm1 = loadtxt(os.environ["PHOT_IONZ_DATA_PATH"] + '/msm1.csv')
        datamsm2 = loadtxt(os.environ["PHOT_IONZ_DATA_PATH"] + '/msm2.csv')

        fig, ax = plt.subplots()
        fig.suptitle('Halo mass function', fontsize=18)
        ax.errorbar(datamsm1[:, 0], datamsm1[:, 1] / datamsm1[:, 0], ls='', yerr=sqrt(datamsm1[:, 1] / datamsm1[:, 0]), # / (500 / 0.7)**3), \
                    marker='o', markersize=5, color='black', label=r'MS-I, $z = 0$')
        ax.errorbar(datamsm2[:, 0], datamsm2[:, 1] / datamsm2[:, 0], ls='', yerr=sqrt(datamsm2[:, 1] / datamsm2[:, 0]), # / (100 / 0.7)**3), \
                    marker='D', markersize=5, color='grey', label=r'MS-II, $z = 0$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        for c in range(numbredsplotlate):
            plot = ax.loglog(mass, mass * diffnhaldiffmass[:, indxredsplotlate[c]] * 1e9, label=strgredslate[c])
        ax.set_ylim([1e-9, 1e5])
        ax.set_ylabel('$dN_h/d\log M$ [1/Mpc$^3$]')
        ax.set_xlabel('$M [M_\odot]$')     
        ax.legend()
        plt.savefig(pathplot + 'diffnhaldiffmass.png')
        plt.close()

    return diffnhaldiffmass, peakhght


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


def retr_fluxphotdmat(csecvelo, csecfrac, masspart, dmatslop):
    
    # run tag
    global rtag
    rtag = propmodl + '_' + concmodl + '_' + subsmodl + '_' + igmamodl + '_massandm%.3g_' % masspart + anch + \
        '_csecvelo%.3g'  % -log10(csecvelo) + '_csecfrac%.3g' % -log10(csecfrac)
    
    # DM annihilation spectrum
    multintptemp = interp1d(massp4dm, multp4dm, axis=1)(masspart) / masspart / enelscalp4dm
    indxeneltemp = where((enel > amin(enelscalp4dm * masspart)) & (enel < amax(enelscalp4dm * masspart)))[0]
    multintp = zeros(numbenel)
    multintp[indxeneltemp] = interp1d(enelscalp4dm * masspart, multintptemp)(enel[indxeneltemp])
    
    # energy density and velocity variance in DM halos
    global edendmathalo, velovarihalo, rvelvarihalo
    velovarihalo = zeros((numbmass, numbreds, numbrsph))
    edendmathalo = zeros((numbmass, numbreds, numbrsph))
    for d in range(numbmass):
        for c in indxreds:
            edendmathalonorm = odenviri * edendmat[c] * conc[d, c]**3 / 3. / (log(1. + conc[d, c]) - conc[d, c] / (1. + conc[d, c])) # [MeV/cm^3]
            rsphscal = rsphviri[d, c] / conc[d, c] # [kpc]
            rsphnorm = rsph[d, c, :] / rsphscal # [1]
            edendmathalo[d, c, :] = edendmathalonorm / rsphnorm / (1. + rsphnorm)**2 # [MeV/cm^3]
            edendemcnorm = (2. * demccons - 3.)**3 / 4. / (5. - 2. * demccons) * edendmathalonorm # [MeV/cm^3]
            rsphdemc = (5. - 2. * demccons) / (2. * demccons - 3.) * rsphscal # [kpc]
            velovarihalo[d, c, :] = 4. * pi * gravconsredu / 200. * 81. * edendemcnorm**(1. / 3.) * rsphdemc**2 * (edendmathalo[d, c, :] \
                * (rsph[d, c, :] / rsphdemc)**demccons)**(2. / 3.) /                 velolght**2 * kprc2cmet**2 # [1] 
    rvelvarihalo = demcsigm * velovarihalo
          
    # DM annihilation cross section
    csecvelohalo = zeros((numbmass, numbreds, numbrsph, numbmpol))
    csecvelohalo[:, :, :, 0] = csecvelo
    csecvelohalo[:, :, :, 1] = csecvelo * csecfrac * rvelvarihalo

    # electron emissivity in DM halo
    global emiselechaloenel, diffemiselechalodiffrsphenel, lumielechaloenel, diffemiselecclmpdiffmassenel, emiselecenel, emiselecclmpenel, \
        emiselecsmthenel, rvelvarihavg, rvelvariiavg, rvelvariiavgsmth, rvelvariiavgclmp, ndiffenelecsavg, ndiffenelecsavgenel, emisphot, \
        emisphotrydb, fluxphot, fluxphotrydb, fluxphotenph, diffionrdiffenph, ionr
        
    emiselechalohost = zeros((numbenel, numbmass, numbreds, numbrsph, numbmpol))
    emiselechalosubh = zeros((numbenel, numbmass, numbreds, numbrsph, numbmpol))
    emiselechalotemp = zeros((numbenel, numbmass, numbmass, numbreds, numbrsph, numbmpol))
    emiselechalo = zeros((numbenel, numbmass, numbreds, numbrsph, numbmpol))
    diffemiselechalodiffrsph = zeros((numbenel, numbmass, numbreds, numbrsph, numbmpol))
    lumielechalo = zeros((numbenel, numbmass, numbreds, numbmpol))
    for a in range(numbenel):
        for d in range(numbmass):
            for c in indxreds:
                for m in range(numbmpol):
                    emiselechalohost[a, d, c, :, m] = csecvelohalo[d, c, :, m] / masspart**2 *    multintp[a] * fesc[a, d, c, :] * edendmathalo[d, c, :]**2 # [1/s/cm^3/MeV]
                    for f in range(numbmass):
                        if f < d:
                            emiselechalotemp[a, d, f, c, :, m] = 0.
                    indxmassplot = where(mass[f] < mass[d])[0]
                    emiselechalosubh[a, d, c, :, m] = trapz(emiselechalotemp[a, d, :, c, :, m], mass, axis=0)
                    emiselechalo[a, d, c, :, m] = emiselechalohost[a, d, c, :, m] + emiselechalosubh[a, d, c, :, m]
                    diffemiselechalodiffrsph[a, d, c, :, m] = 4. * pi * rsph[d, c, :]**2 *    emiselechalo[a, d, c, :, m] * kprc2cmet**2 # [1/s/kpc/MeV]
                    lumielechalo[a, d, c, m] =    trapz(diffemiselechalodiffrsph[a, d, c, :, m] * kprc2cmet, rsph[d, c, :]) # [1/s/MeV]
    emiselechaloenel = trapz(emiselechalo, enel, axis=0) # [1/s/cm^3]
    lumielechaloenel = trapz(lumielechalo, enel, axis=0) # [1/s]
    diffemiselechalodiffrsphenel = trapz(diffemiselechalodiffrsph, enel, axis=0) # [1/s/cm^3]
    
    # mean DM relative velocity
    diffrvelvarihavgdiffrsph = zeros((numbmass, numbreds, numbrsph))
    rvelvarihavg = zeros((numbmass, numbreds))
    for d in range(numbmass):
        for c in indxreds:
            diffrvelvarihavgdiffrsph[d, c, :] = 3. * rsph[d, c, :]**2 * rvelvarihalo[d, c, :] / amax(rsph[d, c, :])**3
            rvelvarihavg[d, c] = trapz(diffrvelvarihavgdiffrsph[d, c, :], rsph[d, c, :]) # [1]

    # spatially averaged electron emissivity    
    ## smooth component
    # temp -- add p-wave to the smth component
    emiselecsmth = zeros((numbenel, numbreds, numbmpol))
    for c in indxreds:
        for m in range(numbmpol):
            emiselecsmth[:, c, m] = csecvelo * edendmat[c]**2 / masspart**2 * multintp # [1/cm^3/s/MeV]
    emiselecsmthenel = trapz(emiselecsmth, enel, axis=0) # [1/cm^3/s]

    ## clumpy component
    diffemiselecclmpdiffmass = zeros((numbenel, numbmass, numbreds, numbmpol))
    for a in range(numbenel):
        for d in range(numbmass):
            for c in indxreds:
                for m in range(numbmpol):
                    diffemiselecclmpdiffmass[a, :, c, m] = diffnhaldiffmass[:, c] * lumielechalo[a, :, c, m] / kprc2cmet**3 * (1. + reds[c])**3 # [1/cm^3/s/MeV/Msun]
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
    ndiffenelecsavg = zeros((numbenel, numbreds, numbmpol))
    for c in indxreds:
        for m in range(numbmpol):
            for i in range(numbenel-1):
                ndiffenelecsavg[i, c, m] = trapz(emiselec[i:numbenel, c, m], enel[i:numbenel], axis=0) / edotegbl[i, c] # [1/cm^3/MeV]
    ndiffenelecsavgenel = trapz(ndiffenelecsavg, enel, axis=0) # [1/cm^3]

    # photon emissivity
    emisphot = trapz(specinvc[:, :, :, None] * ndiffenelecsavg[None, :, :, :], enel, axis=1) # [1/cm^3/s/MeV]
    emisphotrydb = interp1d(enph, emisphot, axis=0)(enerrydb) # [1/cm^3/s/MeV]
    emisphotenph = trapz(emisphot, enph, axis=0) # [1/cm^3/s]

    # photon flux 
    difffluxphotdiffreds = zeros((numbenph, numbreds, numbreds, numbmpol))
    fluxphot = zeros((numbenph, numbreds, numbmpol))
    for a in range(numbenph):
        for c in indxreds:
            enphshft = enph[a] * (1. + reds[c]) / (1. + reds)
            indxredsshft = where((minmenph < enphshft) & (enphshft < maxmenph))[0]
            if indxredsshft.size == 0:
                continue
            for m in range(numbmpol):
                emisphotintp = interp1d(enph, emisphot[:, c, m])(enphshft[indxredsshft]) # [1/cm^3/s/MeV]
                difffluxphotdiffreds[a, c, indxredsshft, m] = velolght * timehubb[indxredsshft] / 4. / pi * emisphotintp * exp(-odep[a, c, indxredsshft]) * \
                    (1. + reds[c])**3 / (1. + reds[indxredsshft])**4 # [1/cm^2/s/sr/MeV]
                fluxphot[a, c, m] = trapz(difffluxphotdiffreds[a, c, :, m], reds) # [1/cm^2/s/sr/MeV]  
    fluxphotrydb = interp1d(enph, fluxphot, axis=0)(enerrydb) # [1/cm^2/s/sr/MeV]
    fluxphotenph = trapz(fluxphot, enph, axis=0) # [1/cm^2/s/sr]

    # photo-ionization rate
    diffionrdiffenph = 4. * pi * csecionz[:, None, None] * fluxphot # [1/s/MeV]
    ionr = trapz(diffionrdiffenph, enph, axis=0) # [1/s]
    
    return ionr, fluxphot
        

def plot_halo():

    listlinestyl = [':', '--', '-']
    listcolr = ['b', 'g', 'r']
    
    # concentration parameter
    fig, ax = plt.subplots()
    for m in range(nconcmodl):
        for c in range(numbredsplotlate):
            ax.loglog(mass, conccatl[:, indxredsplotlate[c], m], ls=listlinestyl[m], color=listcolr[c])
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

    plt.savefig(pathplot + 'halo.png')
    plt.close()

    # substructure
    #fig, ax = plt.subplots()
    #ax.loglog(rsph, subfs)
    #ax.loglog(rsph, subbfrsp)           
    #ax.set_xlabel('$r$ [kpc]')
    #ax.set_ylabel('$f_s(r)$')

    # DM energy density and relative velocity variance in halos
    temp = meshgrid(indxmassplot, indxredsplotlate, indexing='ij')
    labl = strgmass + strgredslate

    fig, ax = plt.subplots()
    fig.suptitle('DM energy density in halos', fontsize=18)                
    plot_matr(ax, rsph[temp[0], temp[1], :], edendmathalo[temp[0], temp[1], :],               labl=labl, loc=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$r$ [kpc]')
    ax.set_ylabel(r'$\rho(r,M,z)$ [MeV/cm$^3$]')
    plt.savefig(pathplot + 'edendmathalo.png')
    plt.close()

    fig, ax = plt.subplots()
    fig.suptitle('DM relative velocity variance in halos', fontsize=18)
    plot_matr(ax, rsph[temp[0], temp[1], :], rvelvarihalo[temp[0], temp[1], :], labl=labl, loc=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$r$ [kpc]')
    ax.set_ylabel('$v^2 (M,z,r)$')
    plt.savefig(pathplot + 'rvelvarihalo.png')
    plt.close()

    # virial and scale radii
    fig, ax = plt.subplots()
    for c in range(numbredsplotlate):
        ax.loglog(mass, rsphviri[:, indxredsplotlate[c]], label=strgredslate[c])
    ax.set_xlabel('$M [M_\odot]$')
    ax.set_ylabel('$r_{vir}(M,z)$ [kpc]')
    ax.set_title('Halo Virial Radius')
    ax.legend(loc=2)
    plt.savefig(pathplot + 'rsphviri.png')
    plt.close()

    # e^-/e^+ emissivity in a halo
    fig, axgrd = plt.subplots(numbmassplot, numbmpol, sharey='row', sharex='all', figsize=(14,21))
    fig.suptitle('$e^-/e^+$ emissivity in a halo', fontsize=18)
    for d, axrow in enumerate(axgrd):
        for p, ax in enumerate(axrow):
            for c in range(numbredsplotlate):
                ax.loglog(rsph[indxmassplot[d], indxredsplotlate[c], :], emiselechaloenel[indxmassplot[d], indxredsplotlate[c], :, p],      label=strgredslate[c])
            ax.set_xlabel('$r$ [kpc]')
            if d == 0:
                ax.set_title(strgmpol[p])
            ax.text(0.23, 0.15, strgmass[d], ha='center', va='center', transform = ax.transAxes)
            if d == numbmassplot / 2:
                ax.set_ylabel(r'$dN_e/dtdV (M,z,r)$ [1/s/cm$^3$]')
    ax.legend()
    plt.savefig(pathplot + 'emiselechaloenel.png')
    plt.close()

    # e^-/e^+ differential luminosity of a halo
    fig, axgrd = plt.subplots(numbmassplot, numbmpol, sharey='row', sharex='all', figsize=(14,21))
    fig.suptitle('$e^-/e^+$ differential luminosity in a halo', fontsize=18)
    for d, axrow in enumerate(axgrd):
        for p, ax in enumerate(axrow):
            for c in range(numbredsplotlate):
                ax.loglog(rsph[indxmassplot[d], indxredsplotlate[c],:], diffemiselechalodiffrsphenel[indxmassplot[d], indxredsplotlate[c],:,p],      label=strgredslate[c])
            ax.set_xlabel('$r$ [kpc]')
            if d == 0:
                ax.set_title(strgmpol[p])
            ax.text(0.23, 0.12, strgmass[d], ha='center', va='center', transform = ax.transAxes)
            if d == numbmassplot / 2:
                ax.set_ylabel(r'$dN_e/dtd\log r$ [1/s]')
    ax.legend()
    plt.savefig(pathplot + 'diffemiselechalodiffrsphenel.png')
    plt.close()

    # e^-/e^+ luminosity of a halo
    fig, axrow = plt.subplots(1, numbmpol, sharey='row', sharex='all', figsize=(14,7))
    fig.suptitle('$e^-/e^+$ luminosity of a DM halo', fontsize=18)
    for p, ax in enumerate(axrow):
        for c in range(numbredsplotlate):
            ax.loglog(mass, lumielechaloenel[:, indxredsplotlate[c],p], label=strgredslate[c])
        ax.set_xlabel(r'$M$ $[M_\odot]$')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel('$dN_e/dt (M,z)$ [1/s]')
    ax.legend(loc=2)
    plt.savefig(pathplot + 'lumielechaloenel.png')
    plt.close()

    # spatially averaged e^-/e^+ differential emissivity
    fig, axrow = plt.subplots(1, numbmpol, sharey='row', figsize=(14,7))
    fig.suptitle('$e^-/e^+$ emissivity due to DM annihilation per decade in halo mass', fontsize=18)
    for p, ax in enumerate(axrow):
        for c in range(numbredsplotlate):
            ax.loglog(mass, mass * diffemiselecclmpdiffmassenel[:, indxredsplotlate[c],p], label=strgredslate[c])
        ax.set_xlabel(r'$M$ $[M_\odot]$')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel(r'$dN_e/dtdVd\log M$ [1/s/cm$^3$]')
            ax.legend()
        ax.set_ylim([1e-36, 1e-28])
    plt.savefig(pathplot + 'diffemiselecclmpdiffmassenel.png')
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
    plt.savefig(pathplot + 'emiselecenel.png')
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
    plt.savefig(pathplot + 'emiselecenelcomo.png')
    plt.close()

    # Mean DM relative velocity variance in a halo
    fig, ax = plt.subplots(figsize=(7, 7))
    for c in range(numbredsplotlate):
        ax.loglog(mass, sqrt(rvelvarihavg[:, indxredsplotlate[c]]), label=strgredslate[c])
    ax.set_title('Mean DM relative velocity variance in a halo', fontsize=18)
    ax.set_ylabel('$\sigma_v$')
    ax.set_xlabel(r'$M$ $[M_\odot]$')
    ax.legend(loc=2)
    plt.savefig(pathplot + 'rvelvarihavg.png')
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
    plt.savefig(pathplot + 'rvelvariiavg.png')
    plt.close()

    
def plot_elec_flux():
    
    # e^-/e^+ number density in the IGM
    fig, axrow = plt.subplots(1, numbmpol, sharey='all', figsize=(14,7))
    fig.suptitle('Differential $e^-/e^+$ number density in the IGM', fontsize=18)
    for p, ax in enumerate(axrow):
        for c in range(numbredsplotlate):
            ax.loglog(enel[0:numbenel-1], enel[0:numbenel-1] * ndiffenelecsavg[0:numbenel-1,indxredsplotlate[c],p],  label=strgredslate[c])
        ax.set_xlabel(r'$E_e$ [GeV]')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel(r'$dN_e/dV/d\log E$ [1/cm$^3$]')
    ax.legend(loc=2)
    plt.savefig(pathplot + 'ndiffenelecsavg.png')
    plt.close()

    # e^-/e^+ number density in the IGM
    fig, axrow = plt.subplots(1, numbmpol, sharey='all', figsize=(14,7))
    fig.suptitle('$e^-/e^+$ number density in the IGM', fontsize=18)
    for p, ax in enumerate(axrow):
        ax.loglog(reds, ndiffenelecsavgenel[:, p])
        ax.set_xlabel('$z$')
        ax.set_title(strgmpol[p])
        if p == 0:
            ax.set_ylabel('$dN_e/dV (z)$ [1/cm$^3$]')
    plt.savefig(pathplot + 'ndiffenelecsavgenel.png')
    plt.close() 
    
        
def plot_invrcomp():

    listlinestyl = ['-', '--', '-.']
    listcolr = ['b', 'g', 'r']
    
    fig, ax = plt.subplots()
    for b in range(numbenelplot):
        for c in range(numbredsplotlate):
            ax.loglog(enph * 1e6, enph * sum(specinvccatl[:, indxenelplot[b], indxredsplotlate[c], :], 1), ls=listlinestyl[b], c=listcolr[c])
    ax.set_xlabel('$E_\gamma$ [eV]')
    ax.set_ylabel(r'$dN_\gamma/dtd\log E_\gamma$ [1/s]')
    ax.set_ylim([1e-15, 1e-7])
    ax.set_title('Spectrum of IC scattered CMB photons')

    listlabl = []
    listlabl.append('$E_e$ = ' + tdpy.util.mexp(enel[indxenelplot[0]] * 1e-3) + ' GeV')
    listlabl.append('$E_e$ = ' + tdpy.util.mexp(enel[indxenelplot[1]] * 1e-3) + ' GeV')
    listlabl.append('$E_e$ = ' + tdpy.util.mexp(enel[indxenelplot[2]] * 1e-3) + ' GeV')
    listlabl.extend(strgredslate)
    
    listline = []
    for k in arange(3):
        listline.append(plt.Line2D((0,1),(0,0), color='black', ls=listlinestyl[k]))
    for k in range(3):
        listline.append(plt.Line2D((0,1),(0,0), color=listcolr[k]))
    ax.legend(listline, listlabl, loc=9, ncol=2) 
    
    plt.savefig(pathplot + 'specinvccatl.png')
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
    plt.savefig(pathplot + 'diffnabsdiffcdendiffreds.png')
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
    plt.savefig(pathplot + 'diffoptddiffcdendiffredsrydb.png')
    plt.close() 
    
    fig, ax = plt.subplots()
    ax.loglog(reds, reds * diffoptddiffredsrydb)
    for k in range(redsbrek.size):
        ax.axvline(redsbrek[k], color='grey', ls='--')
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$d\tau / d\log z$')
    ax.set_ylim([None, 1.])
    ax.set_title('Optical depth per decade in redshift at 1 Rydberg')
    plt.savefig(pathplot + 'diffoptddiffredsrydb.png')
    plt.close() 
        
    fig, ax = plt.subplots()
    for c in range(numbredsplotlate):
        ax.loglog(reds[indxredsplotlate[c]:], odeprydb[indxredsplotlate[c], indxredsplotlate[c]:], label=strgredslate[c])
    for k in range(redsbrek.size):
        ax.axvline(redsbrek[k], color='grey', ls='--')
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$\tau$')
    ax.set_ylim([None, 1.])
    ax.set_title('Optical depth at 1 Rydberg')
    plt.savefig(pathplot + 'odep.png')
    plt.close() 
    
    
def plot_fluxphot():
    
    # differential photon emissivity 
    fig, axrow = plt.subplots(1, 2, sharey='all', figsize=(14,7))
    fig.suptitle('Differential photon emissivity', fontsize=18)
    for i, ax in enumerate(axrow):
        for c in range(numbredsplotlate):
            plot = ax.loglog(enph * 1e6, enph * emisphot[:, indxredsplotlate[c], i], label=strgredslate[c])
        ax.set_xlabel('$E_\gamma$ [eV]')
        ax.set_ylabel(r'$dN_\gamma/dtdVd\log E_\gamma$ [1/s]')
    ax.legend()
    plt.savefig(pathplot + 'emisphot.png')
    plt.close() 
    
    # photon emissivity 
    fig, axrow = plt.subplots(1, 2, sharey='all', figsize=(14,7))
    fig.suptitle('Photon emissivity', fontsize=18)
    for i, ax in enumerate(axrow):
        plot = ax.loglog(reds, enerrydb * emisphotrydb[:, i])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$dN_\gamma/dtdV$ [1/cm$^3$/s]')
    ax.legend()
    plt.savefig(pathplot + 'emisphotrydb.png')
    plt.close() 
    
    # differential photon flux 
    fig, axrow = plt.subplots(1, numbmpol, sharey='all', figsize=(14, 7))
    fig.suptitle('Cosmic UV/X-ray background flux per decade in photon energy')
    for m, ax in enumerate(axrow):
        for c in range(numbredsplotlate):
            plot = ax.loglog(enph * 1e6, enph * fluxphot[:, indxredsplotlate[c], m], label=strgredslate[c])
        ax.set_xlabel('$E_\gamma$ [eV]')
        ax.set_ylabel(r'$dN_\gamma/dtdAd\log E_\gamma$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(pathplot + 'fluxphot.png')
    plt.close() 
    
    # differential photon flux at 1 Ryd
    fig, axrow = plt.subplots(1, numbmpol, sharey='all', figsize=(14, 7))
    fig.suptitle('Photon flux at 1 Rydberg')
    for m, ax in enumerate(axrow):
        plot = ax.loglog(reds, enerrydb * fluxphotrydb[:, m])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$dN_\gamma/dtdAd\log E_\gamma$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(pathplot + 'fluxphotrydb.png')
    plt.close() 
    
    # integrated photon flux 
    fig, axrow = plt.subplots(1, 2, sharey='all', figsize=(14,7))
    fig.suptitle('Integral photon flux')
    for m, ax in enumerate(axrow):
        plot = ax.loglog(reds, fluxphotenph[:, m])
        ax.set_xlabel('$z$')
        ax.set_ylabel(r'$dN_\gamma/dtdA (z)$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(pathplot + 'fluxphotenph.png')
    plt.close() 

        
def plot_ionr():

    fig, axrow = plt.subplots(1, numbmpol, sharey='all', figsize=(14,7))
    fig.suptitle('Differential meta-galactic ionization rate per decade in photon energy', fontsize=18)
    for i, ax in enumerate(axrow):
        for c in indxredsplot:
            plot = ax.loglog(enph * 1e6, diffionrdiffenph[:,indxredsplotlate[c],i], label=strgredslate[c])
        ax.set_xlabel('$E_\gamma [eV]$')
        ax.set_ylabel(r'$d\Gamma/d\log E_\gamma$ [1/s/H]')
        ax.axvline(enerrydb * 1e6, ls='--', color='grey')
    ax.legend()
    plt.savefig(pathplot + 'diffionrdiffenph.png')
    plt.close() 
        
    listname = ['ionrbolt.csv', 'ionrbeck.csv', 'ionrfauc.csv']
    listlabl = ['Bolton & Haehnelt (2007)', 'Becker et al. (2007)', 'Faucher-Giguere (2008)']
    listmrkr = ['o', 'D', 'x']
    
    fig, axrow = plt.subplots(1, numbmpol, sharey='all', figsize=(14,7))
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
    plt.savefig(pathplot + 'ionr.png')
    plt.close() 
        

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


def retr_mocksampvarb():
    
    numbpara = 4
    mocksampvarb = empty(numbpara)
    mocksampvarb[0] = 1e-26
    mocksampvarb[1] = 1e4
    mocksampvarb[2] = 1e5
    mocksampvarb[3] = 1.
    
    return mocksampvarb


def retr_datapara():
    
    numbpara = 4
    
    dictpara = dict()
    minmpara = zeros(numbpara)
    maxmpara = zeros(numbpara)
    namepara = empty(numbpara, dtype=object)
    scalpara = empty(numbpara, dtype=object)
    lablpara = empty(numbpara, dtype=object)
    unitpara = empty(numbpara, dtype=object)
    varipara = zeros(numbpara)

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

    strgpara = namepara + ' ' + lablpara
    datapara = namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara
    
    return datapara


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


def plot_sfrd():
    
    fig, ax = plt.subplots()
    ax.set_title('Star Formation Rate Density')
    
    ax.plot(redssfrd, sfrd)
    ax.set_yscale('log')
    
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'SFRD [erg/Mpc$^3$/yr]')
    plt.savefig(pathplot + 'sfrd.png')
    plt.close()


def plot_hm12(fluxphotdmat=None, listfluxphotdmat=None):
    fig, ax = plt.subplots()
    ax.set_title('UV/X-ray photon background')
    
    ax.fill_between(enph[fenph] * 1e6, enph[fenph] * (fluxphotexpr + sqrt(fluxphotexprvari)), enph[fenph] * (fluxphotexpr - sqrt(fluxphotexprvari)), color='lightblue')
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
    for c in range(numbredsplotlate):
        ax.loglog(enph * 1e6, enph * fluxphothm12[:, indxredsplotlate[c]],     label='Haardt & Madau (2012), ' + strgredslate[c])
        if fluxphotdmat != None:
            ax.loglog(enph * 1e6, enph * fluxphotdmat[:, indxredsplotlate[c]], label='DM, ' + strgredslate[c])
        if listfluxphotdmat != None:
            tdpy.mcmc.plot_braz(ax, enph * 1e6, enph[None, :] * listfluxphotdmat[:, :, indxredsplotlate[c]],            lcol=listcolr[c], alpha=0.5, dcol=listcolr[c], mcol='black')
    ax.set_xlabel(r'$E_\gamma$ [eV]')
    ax.set_ylabel(r'$EdN/dE$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(pathplot + 'fluxphothm12.png')
    plt.close()


def init(cnfg):
    
    global makeplot, verbtype
    makeplot = cnfg['makeplot']
    verbtype = cnfg['verbtype']
    
    numbswep = cnfg['numbswep']
    numbburn = cnfg['numbburn']
    factthin = cnfg['factthin']

    global pathplot
    pathbase = os.environ["PHOT_IONZ_DATA_PATH"]
    pathplot = pathbase + '/png/'

    global minmpara, maxmpara, scalpara, namepara, lablpara, unitpara, varipara, dictpara, ipara, numbpara
    mocksampvarb = retr_mocksampvarb()
    datapara = retr_datapara()
    namepara, strgpara, minmpara, maxmpara, scalpara, lablpara, unitpara, varipara, dictpara = datapara
    numbpara = len(lablpara)
    indxpara = arange(numbpara)
    
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
    global numbreds, numbradi, numbrsph, numbzaxi, numbenph, numbenel, numbenpi, numbcden, numbmass, numbwnum
    global reds, radi, rsph, zaxi, enph, enel, enpi, cden, mass, massprim, wnum
    global indxredsplot, indxredsplotlate, indxmassplot, indxenelplot
    global numbredsplot, numbredsplotlate, numbmassplot, numbenelplot
    global strgredslate 
    
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
    tempcmbrnunc = 2.725 # CMB temperature today [K]
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
    numbradi = 50
    numbrsph = 50
    numbzaxi = 50
    numbenph = 50
    numbenel = 50
    numbenpi = 50
    numbcden = 50
    numbreds = 50
    numbmass = 50
    numbwnum = 500
  
    global indxradi, indxrsph, indxzaxi, indxenph, indxenel, indxenpi, indxcden, indxreds, indxmass, indxwnum
    indxradi = arange(numbradi)
    indxrsph = arange(numbrsph)
    indxzaxi = arange(numbzaxi)
    indxenph = arange(numbenph)
    indxenel = arange(numbenel)
    indxenpi = arange(numbenpi)
    indxcden = arange(numbcden)
    indxreds = arange(numbreds)
    indxmass = arange(numbmass)
    indxwnum = arange(numbwnum)

    global minmenph, maxmenph
    minmenph = 1e-6 # [MeV]
    maxmenph = 1e-1 # [MeV]
    enph = logspace(log10(minmenph), log10(maxmenph), numbenph)
    jenph = [0, numbenph/2, numbenph-1]
    njenph = len(jenph)

    global minmenel, maxmenel
    minmenel = 5e1 # [MeV]
    maxmenel = 1e5 # [MeV]
    enel = logspace(log10(minmenel), log10(maxmenel), numbenel)
    diffenel = enel[1:numbenel-1] - enel[0:numbenel-2]
    indxenelplot = [0, numbenel / 4, numbenel / 2]
    numbenelplot = len(indxenelplot)

    global minmreds, maxmreds
    minmreds = 1e-1
    maxmreds = 1e2
    reds = logspace(log10(minmreds), log10(maxmreds), numbreds)
    
    redsprox = [0., 1e3, 1e6]
    numbredsplot = len(redsprox)
    indxredsplot = []
    for k in indxredsplot:
        indxredsplot.append(argmin(abs(reds - redsprox[k])))
        
    redslateprox = [0.1, 2., 6.]
    numbredsplotlate = len(redslateprox)
    indxredsplotlate = []
    strgredslate = []
    for k in range(numbredsplotlate):
        indxredsplotlate.append(argmin(abs(reds - redslateprox[k])))
        strgredslate.append('$z = %.3g$' % redslateprox[k])
        
    global minmmass, maxmmass
    minmmass = 1e8 # [Solar Mass]
    maxmmass = 1e16 # [Solar Mass]
    massprim = logspace(log10(minmmass), log10(maxmmass), numbmass + 1)
    mass = massprim[:-1]
    
    massprox = [1e10, 1e12, 1e15]
    numbmassplot = len(massprox)
    indxmassplot = []
    global strgmass
    strgmass = []
    for d in range(numbmassplot):
        indxmassplot.append(argmin(abs(mass - massprox[d])))
        strgmass.append('$M$ = ' + tdpy.util.mexp(massprox[d]) + r' $M_\odot$')
            
    minmcden = 10**11. # [1/cm^2]
    maxmcden = 10**22. # [1/cm^2]
    cden = logspace(log10(minmcden), log10(maxmcden), numbcden)

    minmenpi = 1e-12 # [MeV]
    maxmenpi = 1e-5 # [MeV]
    enpi = logspace(log10(minmenpi), log10(maxmenpi), numbenpi)

    # wavenumber axis
    minmwnum = 1e-4
    maxmwnum = 1e4
    wnum = logspace(log10(minmwnum), log10(maxmwnum), numbwnum)
    
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
    rsph = zeros((numbmass, numbreds, numbrsph))
    for c in indxreds:
        for d in range(numbmass):             
            rsph[d, c, :] = logspace(log10(minmrsph[d,c]), log10(maxmrsph[d,c]), numbrsph) # [kpc]
    
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
        numbmasspart = masspart.size

        global multp4dm, enelscalp4dm, massp4dm
        numbanch = len(anchlist)
        fig, ax = plt.subplots()
        for a, anch in enumerate(anchlist):
            multp4dm, enelscalp4dm, massp4dm = tdpy.util.retr_p4dm_spec(anch)
            for k in range(numbmasspart):
                indxmassplotp4dm = argmin(abs(massp4dm - masspart[k]))
                ax.loglog(enelscalp4dm * masspart[k] * 1e-3, multp4dm[:, indxmassplotp4dm], color=colorlist[a], ls=linelist[k], label=anchlabl[a] + ', ' + masspartlabl[k])
        ax.set_xlim([0.05, 1e3]) 
        ax.set_ylim([1e-1, 5e2])
        ax.set_xlabel('$E$ [GeV]')
        ax.set_ylabel(r'$dN/d\log E$')
        ax.set_title('Prompt product spectrum per annihilation')
        ax.legend(ncol=4, loc=2)
        plt.savefig(pathplot + 'multp4dm.png')
        plt.close()
    
    multp4dm, enelscalp4dm, massp4dm = tdpy.util.retr_p4dm_spec(anch)

    #if makeplot:
        #plot_edot
        #plot_galprop

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
    edencmbr = retr_edencmbr() # [1/cm^3/MeV]
    edencmbrenpi = trapz(enpi[:, None] * edencmbr, enpi, axis=0) # [MeV/cm^3]
    
    # EGBL energy density
    global edenegbl
    edenegbl = retr_edenegbl() # [1/cm^3/MeV]
    edenegblenpi = trapz(enpi[:, None] * edenegbl, enpi, axis=0) # [MeV/cm^3]
    
    # Energy loss rate on EGBL
    global edotegbl
    # temp
    edotegbl = 4. * csecthom * velolght / 3. / masselec**2 * edencmbrenpi[None, :] * enel[:,None]**2 # [MeV/s]
 
    if makeplot:
        fig, ax = plt.subplots()
        for c in range(numbredsplotlate):
            plot = ax.loglog(enpi * 1e6, enpi * edencmbr[:, indxredsplotlate[c]], label=strgredslate[c])
            plot_ = ax.loglog(enpi * 1e6, enpi * edenegbl[:, indxredsplotlate[c]], '--', color=plot[0].get_color())
        ax.set_xlabel('$E_\gamma$ [eV]')
        ax.set_ylabel(r'$dN_\gamma/dVd\log E_\gamma$ [1/cm$^3$]')
        ax.set_title('Extragalactic Background Radiation')
        ax.set_ylim([1e-1, 1e6])
        leg = ax.legend(loc=2)
        plt.savefig(pathplot + 'edencmbr.png')
        plt.close()

    # halo mass function
    global diffnhaldiffmass, peakhght
    diffnhaldiffmass, peakhght = retr_hmfn()
    
    global numbmpol    
    numbmpol = 2
    
    global propmodl, concmodl, subsmodl, igmamodl, conc

    propmodl = 'effi'
    concmodl = 'duff'
    subsmodl = 'none'
    igmamodl = 'clum'

    # halo concentration model
    global conccatl, conc, nconcmodl
    nconcmodl = 3
    conccatl = zeros((numbmass, numbreds, nconcmodl))

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
    subs = ones((numbmass, numbreds, numbrsph))
    if subsmodl == 'none':
        subs[:] = 1.
    
    # cosmic ray escape model
    uisrf = 6e-7 # [MeV/cm^3]
    magfd = 5. # [muG]
    gasdn = 1. # [1/cm^3]

    # cosmic raw propagation model
    global fesc
    fesc = zeros((numbenel, numbmass, numbreds, numbrsph))
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
    odep = zeros((numbenph, numbreds, numbreds))
    odeprydb = zeros((numbreds, numbreds))
    if igmamodl == 'clum':
        
        global diffnabsdiffcdendiffreds, diffoptddiffcdendiffredsrydb, diffoptddiffredsrydb
        diffnabsdiffcdendiffreds = zeros((numbcden, numbreds))

        # lower Lyman - alpha forest
        indxcdentemp = where(cden < 10**15.)[0]

        indxredstemp = where((reds > 1.56) & (reds < 5.5))[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**7.079 * cden[indxcdentemp, None]**(-1.5) * (1. + reds[None, indxredstemp])**3

        indxredstemp = where(reds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**8.238 * cden[indxcdentemp, None]**(-1.5) * (1. + reds[None, indxredstemp])**0.16

        indxredstemp = where(reds > 5.5)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**1.470 * cden[indxcdentemp, None]**(-1.5) * (1. + reds[None, indxredstemp])**9.9

        # upper Lyman - alpha forest
        indxcdentemp = where((cden > 10**15.) & (cden < 10**17.5))[0]
        
        indxredstemp = where((reds > 1.56) & (reds < 5.5))[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**14.58 * cden[indxcdentemp, None]**(-2.) * (1. + reds[None, indxredstemp])**3.

        indxredstemp = where(reds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**15.74 * cden[indxcdentemp, None]**(-2.) * (1. + reds[None, indxredstemp])**0.16

        indxredstemp = where(reds > 5.5)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**8.97 * cden[indxcdentemp, None]**(-2.) * (1. + reds[None, indxredstemp])**9.9

        # Super Lyman limit systems
        indxcdentemp = where((cden > 10**19.) & (cden < 10**20.3))[0]

        indxredstemp = where(reds > 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**(-0.347) * cden[indxcdentemp, None]**(-1.05) * (1. + reds[None, indxredstemp])**1.27

        indxredstemp = where(reds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**(0.107) * cden[indxcdentemp, None]**(-1.05) * (1. + reds[None, indxredstemp])**0.16

        # Damped Lyman alpha systems
        indxcdentemp = where(cden > 10**20.3)[0]

        indxredstemp = where(reds > 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**18.94 * cden[indxcdentemp, None]**(-2.) * (1. + reds[None, indxredstemp])**1.27

        indxredstemp = where(reds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, indxredstemp, indexing='ij')
        diffnabsdiffcdendiffreds[indxtemp] = 10**19.393 * cden[indxcdentemp, None]**(-2.) * (1. + reds[None, indxredstemp])**0.16

        diffoptddiffcdendiffreds = diffnabsdiffcdendiffreds[:, None, :] * (1. - exp(-cden[:, None, None] * csecionz[None, :, None]))
        diffoptddiffreds = trapz(diffoptddiffcdendiffreds, cden, axis=0)
        diffoptddiffredsrydb = interp1d(enph, diffoptddiffreds, axis=0)(enerrydb)
        diffoptddiffcdendiffredsrydb = interp1d(enph, diffoptddiffcdendiffreds, axis=1)(enerrydb)
        odep = zeros((numbenph, numbreds, numbreds))
        for c in indxreds:
            for g in range(c + 1, numbreds):
                odep[:, c, g] = trapz(diffoptddiffreds[:, c:g], reds[c:g], axis=1)
        odeprydb = interp1d(enph, odep, axis=0)(enerrydb)
        
        if makeplot:
            plot_igma()
                
    # compute the differential ICS power as a function of final photon and electron energies
    global specinvc, specinvccatlenph, specinvccatl
    specinvccatl = zeros((numbenph, numbenel, numbreds, 2))
    nqics = 100
    maxmqics = 1.
    for b in range(numbenel):
        for a in range(numbenph):
            elecrelefact = enel[b] / masselec
            eps = enph[a] / enel[b]
            if eps >= 1.:
                continue
            for c in indxreds:
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

    numbproc = 1
    thissampvarb = zeros(numbpara)
    thissampvarb[0] = csecvelopivt
    thissampvarb[1] = csecfracpivt
    thissampvarb[2] = masspartpivt
    thissampvarb[3] = dmatsloppivt
    thissamp = empty((numbproc, numbpara))
    thissamp[0, :] = tdpy.mcmc.cdfn_samp(thissampvarb, datapara)

    numbsamp = (numbswep - numbburn) / factthin
    numbplotside = numbpara
    sampbund = tdpy.mcmc.init(numbproc, numbswep, retr_llik, datapara, thissamp=thissamp, numbburn=numbburn, \
                factthin=factthin, optiprop=optiprop, verbtype=verbtype, pathbase=pathbase, rtag=rtag, numbplotside=numbplotside)
    
    listsampvarb = sampbund[0]
    listsamp = sampbund[1]
    listsampcalc = sampbund[2]
    listllik = sampbund[3]
    listaccp = sampbund[4]
    listjsampvari = sampbund[5]
    
    listfluxphotdmat = listsampcalc[0]
            
    fig, ax = plt.subplots()
    ax.set_title('UV/X-ray photon background')
    
    ax.fill_between(enph[fenph] * 1e6, enph[fenph] * (fluxphotexpr + sqrt(fluxphotexprvari)), enph[fenph] * (fluxphotexpr - sqrt(fluxphotexprvari)), color='lightblue')
    
    ax.loglog(enph[fenph] * 1e6, enph[fenph] * fluxphotexpr, label='ROSAT')
    
    ax.loglog(enph * 1e6, enph * fluxphothm12[:, 0], label='Haardt & Madau (2012)')
    tdpy.mcmc.plot_braz(ax, enph * 1e6, enph[None, :] * listfluxphotdmat[:, :, 0], alpha=0.5, mcol='black')

    ax.set_xlabel(r'$E_\gamma$ [eV]')
    ax.set_ylabel(r'$EdN/dE$ [1/cm$^2$/s/sr]')
    ax.legend()
    plt.savefig(pathplot + 'fluxphotpost.png')
    plt.close()
    
    strgpara = lablpara + ' ' + unitpara
    path = os.environ["PHOT_IONZ_DATA_PATH"] + '/png/mcmc'
    tdpy.mcmc.plot_grid(listsampvarb, strgpara, scalpara=scalpara, path=path, quan=True)

    for k in indxpara:
        path = pathplot + namepara[k] + '.png'
        tdpy.mcmc.plot_trac(listsampvarb[:, k], strgpara[k], scalpara=scalpara[k], path=path, quan=True)


def retr_cnfg( \
              datatype='inpt', \
              datalabl='PIXIE', \
              numbswep=100, \
              numbburn=None, \
              factthin=None, \
              plotperd=10000, \
              verbtype=1, \
              makeplot=False, \
             ):
        
    cnfg = dict()
    
    # data type and label
    cnfg['datatype'] = datatype
    cnfg['datalabl'] = datalabl

    # sampler setup
    cnfg['numbswep'] = numbswep
    cnfg['numbburn'] = numbburn
    cnfg['factthin'] = factthin

    #
    cnfg['plotperd'] = plotperd
    cnfg['verbtype'] = verbtype
    cnfg['makeplot'] = makeplot
    
    return cnfg


def cnfg_nomi():
    
    cnfg = retr_cnfg( \
                     numbswep=100000, \
                     verbtype=0, \
                     numbburn=0, \
                     factthin=1, \
                     makeplot=True, \
                    )
    
    init(cnfg)
    

if __name__ == '__main__':
    
    cnfg_nomi()

