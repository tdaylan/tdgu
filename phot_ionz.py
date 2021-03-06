from __init__ import *

# determine the thermal state of the IGM due to photoheating
## blazar heating rate
    
    #plheatbl = 10**(0.0315 * (1. + red)**3 - 0.512 * (1. + red)**2 + 2.27 * (1. + red) - 2.38) / 3.154e22 # [MeV/s]
    #plheatbl[where(gdat.meanreds gt 5.7)] = 0.
    
## photoheating rate
    
    #  phheatdm = zeros(gdat.numbreds)
    #  phheatqg = zeros(gdat.numbreds)
    #  phheatqghm = zeros(gdat.numbreds)
    #  phheatqs = zeros(gdat.numbreds)
    #  for c=0, nred-1 do begin
    #    phheatdm[c] = trapz(gdat.meanenph, pucr_phheatint(phdmflux[:,c]), xr=[ryd, gdat.meanenph[gdat.numbenph-1]]) # [MeV/s]
    #    phheatqg[c] = trapz(gdat.meanenph, pucr_phheatint(phqgflux[:,c]), xr=[ryd, gdat.meanenph[gdat.numbenph-1]]) # [MeV/s]
    #    phheatqghm[c] = trapz(gdat.meanenph, pucr_phheatint(phqgfluxhm[:,c]), xr=[ryd, gdat.meanenph[gdat.numbenph-1]]) # [MeV/s]
    #    phheatqs[c] = trapz(gdat.meanenph, pucr_phheatint(phqsflux[:,c]), xr=[ryd, gdat.meanenph[gdat.numbenph-1]]) # [MeV/s]
    #  endfor
    

    #if gdat.makeplot:
    #    plot_heat
    
    
    # density - temperature relation
    
    #  npatch = 100
    #  initreds = 20.
    #  mindif = min(abs(gdat.meanreds - initred), jredtemp)
    #  temp = cmbt
    #
    #  puc_lambda
    #  heat = phheatdm
    #  temppatchdm = zeros((npatch, gdat.numbreds))
    # dpatchdensdzarr = zeros((npatch, gdat.numbreds))
    #  patchdensarr = zeros((npatch, gdat.numbreds))
    #  for i=0, npatch-1 do begin
    #    patchdensarr[i,:] = 1. / ((1. + lambda[i,0] * grwf) * (1. + lambda[i,1] * grwf) * (1. + lambda[i,2] * grwf)) - 1.
    #    dpatchdensdzarr[i,:] = deriv(red, patchdensarr[i,:])
    #    patchdens = patchdensarr[i,:]
    #    dpatchdensdz = dpatchdensdzarr[i,:]
    #    temppatchdm[i,:] = puc_temp_solve
    #  endfor    

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
    #if gdat.makeplot:
    #    plot_temp 


def plot_edot(gdat):

    sec2gyr = 3.171e-17 # [s/Gyrs]
    
    ndengas = [0.01, 1., 100.] # [1/cm^3]
    edengas = [1e-8, 1e-6, 1e-4] # [MeV/cm^3]
    intsmag = [1e-10, 1e-8, 1e-6] # [muG]
    
    time = gdat.meanenel[:,None] / retr_edot(intsmag, ndengas, edenrad) * sec2gyr
    timediff = zscldif**2 / diffnor / (enel * 1e-3)**diffind
    
    labels = []
    for a in range(nndengas):
        labels.append([r'Brem, $\rho_{gas}$ = %.2g  cm$^{-3}$' % ndengas[a]])
    for a in range(nedenrad):
        labels.append(['ICS, $n_{isrf}$ = %.2g  1/cm$^{-3}$' % edenrad[a]])
    for a in range(nintsmag):
        labels.append([r'Synch, $B_$ = %.2g $\mu$G' % intsmag[a]])
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.loglog(enel, time)
    axis.set_xlabel('$E_e$ [MeV]')
    axis.set_ylabel(r'$\tau$ [Gyrs]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'edotelec.pdf')
    plt.close()


def retr_edenegbl(gdat):
    
    path = gdat.pathdata + 'egbl.csv'
    egbldata = loadtxt(path)
    wlenegbl = egbldata[:, 0] * 1e-4 # [cm]
    freqegbl = flipud(gdat.velolght / wlenegbl) # [Hz]
    enpiegbl = gdat.plnkcons * freqegbl # [MeV]
    fluxegbl = egbldata[:, 1] # [W/m^2/sr]
    edenegbltemp = flipud(fluxegbl) * 2. * pi / gdat.velolght / 1e4 * 1.6e13 / enpiegbl**2 # [1/cm^3/MeV]
    edenegbl = zeros((gdat.numbenpi, gdat.numbreds))
    for c in gdat.indxreds:
        enpitemp = gdat.meanenpi / (1. + gdat.meanreds[c]) # [MeV]
        indxenpitemp = where((enpitemp < max(enpiegbl)) & (enpitemp > min(enpiegbl)))[0]
        edenegbl[indxenpitemp, c] = interp1d(enpiegbl, edenegbltemp)(enpitemp[indxenpitemp]) * (1. + gdat.meanreds[c])**2 # [1/cm^3/MeV]

    return edenegbl
        

def retr_edencmbr(gdat):
    
    edencmbr = 8. * pi * gdat.meanenpi[:, None]**2 / (gdat.velolght * gdat.plnkcons)**3 / (exp(gdat.meanenpi[:, None] / \
                        gdat.tempcmbrnunc / gdat.boltcons / (1. + gdat.meanreds[None, :])) - 1.) # [1/cm^3/MeV]
    
    return edencmbr


def retr_psec(gdat, thiswnum):

    q = thiswnum / 0.15
    cq = 14.4 + 325 / (1. + 60.5 * q**1.11)
    lq = log(exp(1.) + 1.84 * q)
    tranfunc = lq / (lq + cq * q**2)
    psecprim = 2e-9 * thiswnum**gdat.psecindx
    psec = psecprim * tranfunc**2
    
    return psec, tranfunc


def retr_llik(sampvarb, gdat):

    gdat.csecvelo = sampvarb[0]
    gdat.csecfrac = sampvarb[1]
    gdat.masspart = sampvarb[2]
    #gdat.dmatslop = sampvarb[3]
  
    if gdat.saveflux:
        fluxphotdmattotl = gdat.fluxphotdmat[:, :, 0] * gdat.csecvelo / gdat.csecvelopivt + gdat.fluxphotdmat[:, :, 1] * gdat.csecfrac / gdat.csecfracpivt
    else:
        retr_fluxphotdmat(gdat)
        fluxphotdmattotl = sum(gdat.fluxphotdmat, 2)
    fluxphotmodl = gdat.fluxphothm12[gdat.indxenphexpr, 0] + fluxphotdmattotl[gdat.indxenphexpr, 0]
    lpos = sum(-log(sqrt(2. * pi * gdat.fluxphotexprvari)) - 0.5 * (gdat.fluxphotexpr - fluxphotmodl)**2 / gdat.fluxphotexprvari)
    
    if False:
        print 'hey'
        print 'csecvelo'
        print csecvelo
        print 'csecfrac'
        print csecfrac
        print 'masspart'
        print masspart
        print 'gdat.fluxphothm12'
        print gdat.fluxphothm12
        print 'gdat.fluxphotdmat[:, :, 0]'
        print gdat.fluxphotdmat[:, :, 0]
        print 'fluxphotmodl'
        print fluxphotmodl
        print 'gdat.fluxphotexpr'
        print gdat.fluxphotexpr
        print 'lpos'
        print lpos
        print

    sampcalc = [fluxphotdmattotl[:, 0]]

    return lpos, sampcalc


def retr_hmfn(gdat):
      
    # growth factor
    gdat.grwf = zeros(gdat.numbreds) # growth factor
    for c in gdat.indxreds:
        diffgrwfdiffreds = gdat.funchubb[c] * (1. + gdat.meanreds) / gdat.funchubb**3
        gdat.grwf[c] = trapz(diffgrwfdiffreds[c:], gdat.meanreds[c:])
    # temp -- avoid later division by zero
    gdat.grwf[-1] = gdat.grwf[-2]
    gdat.grwf /= gdat.grwf[0]
    
    # radius, wavelength and wavenumber corresponding to the halo mass
    gdat.rsphhalo = (3. * gdat.meanmassprim * gdat.solm2mgev / 4. / pi / gdat.omegdmat / gdat.edencritnunc / gdat.odenviri)**(1./3.) / gdat.kprc2cmet / 1e3 # [Mpc]
    wlenhalo = 4. * gdat.rsphhalo
    wnumhalo = 2. * pi / wlenhalo

    # power spectrum of density fluctuations
    psec, tranfunc = retr_psec(gdat, gdat.meanwnum)

    # RMS density fluctuations
    fluc = zeros(gdat.numbmass + 1)
    diffflucdiffwnum = zeros((gdat.numbwnum, gdat.numbmass + 1))
    funcwndw = zeros((gdat.numbwnum, gdat.numbmass + 1))
    for d in range(gdat.numbmass + 1):
        wang = gdat.meanwnum * wlenhalo[d]
        funcwndw[:, d] = 3. * (sin(wang) - wang * cos(wang)) / wang**3
        diffflucdiffwnum[:, d] = gdat.meanwnum**2 * psec * funcwndw[:, d]**2 / 2. / pi**2
        fluc[d] = sqrt(trapz(diffflucdiffwnum[:, d], gdat.meanwnum, axis=0))
    fluc *= 0.7 / interp(1e14, gdat.meanmassprim, fluc)
    
    # halo gdat.meanmass function
    difflogtflucdiffmass = -diff(log(fluc)) / diff(gdat.meanmassprim)
    
    peakhght = gdat.odencoll / gdat.grwf[None, :] / fluc[:-1, None]
    
    funcflucprsh = sqrt(2. / pi) * peakhght * exp(-gdat.shtrwgth * peakhght**2 / 2.)
    funcfluc = gdat.shtrnorm * sqrt(2. * gdat.shtrwgth / pi) * (1. + (1. / peakhght**2 / gdat.shtrwgth)**gdat.shtrindx) * \
                                                                                                                    peakhght * exp(-gdat.shtrwgth * peakhght**2 / 2.)
    
    gdat.halodistprsh = funcflucprsh * gdat.edendmat[None, :] * gdat.kprc2cmet**3 / gdat.solm2mgev / gdat.meanmass[:, None] * difflogtflucdiffmass[:, None] # [1/kpc^3/Msun]
    gdat.halodist = funcfluc * gdat.edendmat[None, :] * gdat.kprc2cmet**3 / gdat.solm2mgev / gdat.meanmass[:, None] * difflogtflucdiffmass[:, None] # [1/kpc^3/Msun]
    
    if gdat.makeplot:
        
        # growth factor from Loeb2006
        path = gdat.pathdata + '/odencollloeb.csv'
        odencollloeb = loadtxt(path)
        redsloeb = odencollloeb[:, 0]
        odencollloeb = odencollloeb[:, 1]
        grwfloeb = gdat.odencoll / odencollloeb
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.loglog(gdat.meanreds, gdat.grwf)
        axis.set_xlabel(r'$z$')
        axis.set_ylabel('$D(z)$')
        axis.loglog(redsloeb, grwfloeb, label='Loeb2006')
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'grwf.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, 0.5 * gdat.plotsize))
        axis.loglog(gdat.meanmassprim, gdat.rsphhalo)
        axis.set_xlabel(r'$M [M_\odot]$')
        axis.set_ylabel('$r_H$ [Mpc]')
        axistwin = axis.twinx()
        axistwin.loglog(gdat.meanmassprim, wnumhalo, ls='--')
        axistwin.set_ylabel('$k_H$ [Mpc$^{-1}$]')
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'rsphhalo.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.set_xlabel('$k$ [Mpc$^{-1}$]')
        axis.set_ylabel('$P(k)$ [Mpc$^3$]')
        axis.loglog(gdat.meanwnum, psec)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'psec.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.set_xlabel('$k$ [Mpc$^{-1}$]')
        axis.set_ylabel('$T(k)$')
        axis.loglog(gdat.meanwnum, tranfunc)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'tranfunc.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.set_xlabel('$k$ [Mpc$^{-1}$]')
        axis.set_ylabel('$W(kR)$')
        for d in range(gdat.numbmassplot):
            axis.loglog(gdat.meanwnum, funcwndw[:, gdat.indxmassplot[d]]**2, label=gdat.strgmass[d])
        axis.legend(loc=3)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'funcwndw.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.set_xlabel('$k$ [Mpc$^{-1}$]')
        axis.set_ylabel('$d\sigma^2/dk$')
        for d in range(gdat.numbmassplot):
            axis.loglog(gdat.meanwnum, diffflucdiffwnum[:, gdat.indxmassplot[d]], label=gdat.strgmass[d])
        axis.legend(loc=3)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'diffflucdiffwnum.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.loglog(gdat.meanmass, difflogtflucdiffmass)
        axis.set_xlabel(r'$M [M_\odot]$')
        axis.set_ylabel(r'd$\log\sigma$/d$M$')
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'difflogtflucdiffmass.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.loglog(gdat.meanmassprim, fluc)
        axis.set_xlabel(r'$M [M_\odot]$')
        axis.set_ylabel(r'$\sigma$')
        for c in range(gdat.numbredsplotlate):
            axis.axhline(gdat.odencoll / gdat.grwf[gdat.indxredsplotlate[c]], ls='--', label='$\delta_c(z=%.3g)$' % gdat.meanreds[gdat.indxredsplotlate[c]])
        axis.legend(loc=3)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'fluc.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        for c in range(gdat.numbredsplotlate):
            axis.loglog(fluc[:-1], funcfluc[:, gdat.indxredsplotlate[c]], label=gdat.strgredsplotlate[c])
        axis.set_xlabel(r'$\sigma$')
        axis.set_ylabel('$f(\sigma)$')
        axis.legend(loc=4)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'funcfluc.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
        
        # plot halo mass functions for comparison
        ## Guo2010
        path = gdat.pathdata + '/halodist_msm1.csv'
        datamsm1 = loadtxt(path)
        axis.errorbar(datamsm1[:, 0], datamsm1[:, 1] / datamsm1[:, 0], ls='', yerr=sqrt(datamsm1[:, 1] / datamsm1[:, 0]), # / (500 / 0.7)**3), \
                    marker='o', markersize=5, color='black', label=r'MS-I, $z = 0$')
        path = gdat.pathdata + '/halodist_msm2.csv'
        datamsm2 = loadtxt(path)
        axis.errorbar(datamsm2[:, 0], datamsm2[:, 1] / datamsm2[:, 0], ls='', yerr=sqrt(datamsm2[:, 1] / datamsm2[:, 0]), # / (100 / 0.7)**3), \
                    marker='D', markersize=5, color='grey', label=r'MS-II, $z = 0$')
        axis.set_xscale('log')
        axis.set_yscale('log')
        for c in range(gdat.numbredsplotlate):
            plot = axis.loglog(gdat.meanmass, gdat.meanmass * gdat.halodist[:, gdat.indxredsplotlate[c]] * 1e9, label=gdat.strgredsplotlate[c])
        
        ## Loeb2006
        redstemp = array([0., 5., 10., 20., 30.])
        for k in range(redstemp.size):
            patht = gdat.pathdata + '/halodistloebz%03d.csv' % redstemp[k]
            datatemp = loadtxt(path)
            masstemp = datatemp[:, 0]
            halodisttemp = datatemp[:, 1] / datatemp[:, 0]
            plot = axis.loglog(masstemp, masstemp * halodisttemp, label='Loeb2006, z = %.3g' % redstemp[k], ls='--')
        
        axis.set_ylim([1e-9, 1e5])
        axis.set_ylabel('$dN_h/d\log M$ [1/Mpc$^3$]')
        axis.set_xlabel('$M [M_\odot]$')     
        axis.legend(ncol=2)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'halodist.pdf')
        plt.close()


def retr_fluxphothm12(gdat):
    
    # Haardt & Madau 2012 quasar + galaxy UV/X-ray background flux
    path = gdat.pathdata + 'fluxphothm12.dat'
    tabl = loadtxt(path)

    wlenhm12 = tabl[1:576, 0] # [A]
    freqhm12 = gdat.velolght / flipud(wlenhm12) * gdat.cmet2angs # [Hz]
    gdat.enphhm12 = gdat.plnkcons * freqhm12 # [MeV]
    gdat.redshm12 = tabl[0, :60]
    gdat.fluxphothm12 = flipud(tabl[1:576, 1:61]) # [erg/cm^2/s/Hz/sr]
    gdat.fluxphothm12 *= gdat.ergs2mgev / gdat.plnkcons / gdat.enphhm12[:, None] # [1/cm^2/s/sr/MeV]
    
    gdat.enphhm12, gdat.indxenphhm12 = unique(gdat.enphhm12, return_index=True)
    gdat.fluxphothm12 = gdat.fluxphothm12[gdat.indxenphhm12, :]

    gdat.fluxphothm12 = interp2d(gdat.enphhm12, gdat.redshm12, gdat.fluxphothm12.T)(gdat.meanenph, gdat.meanreds).T

    
    # quasar photon emissivity
    #  phqsemis1 = 10**24.6 * erg2mev / (gdat.kprc2cmet*1d3)**3 * (1. + gdat.meanreds)**7.68 * exp(-0.28 * gdat.meanreds) / (exp(1.77 * gdat.meanreds) + 26.3) # [1/cm^3/s/Hz/MeV]
    #  phqsemis = ((gdat.meanenph / gdat.enerrydb)**(-0.44) * frph) # phqsemis1 # [/cm^3/s/MeV]
    #  ienph = where(gdat.meanenph gt gdat.plnkcons * gdat.velolght / 1.3e-5)
    #  phqsemis[ienph,:] = ((gdat.meanenph[ienph] / gdat.enerrydb)**(-1.57) * frph) # phqsemis1 # [1/cm^3/s/MeV]
    #
    #
    # Haardt & Madau 2012 quasar + galaxy UV/X-ray background emissivity

    #path = gdat.pathdata + 'qgemis.dat'
    #hmemis = read_ascii(path)
    #wlemishm = phqgemisdat.field01[0,1:378] # [Angstrom]
    #fremishm = reverse(reform(gdat.velolght / (wlemishm * 1e-8))) # [Hz]
    #enemishm = gdat.plnkcons * fremishm # [MeV]
    #phqgemis = phqgemisdat.field01[1:59,1:378] # [erg/Mpc^3/s/Hz]
    #phqgemis = transpose(phqgemis)
    #phqgemis = reverse(phqgemis,1)
    #phqgemis *=  2.12e-68 * fremishm * (1. + zeros(60)) # [1/cm^3/s/MeV]
    #phemis[:,:,1] = interpolate(phqgemis, interpol(indgen(378), enemishm, gdat.meanenph), interpol(indgen(59), redhm, red)) 
    #phemisryd[:,1] = interpolate(phqgemis, interpol(indgen(378), enemishm, 1. * gdat.enerrydb), interpol(indgen(59), redhm, red)) 


def retr_fluxphotdmat(gdat):
    
    # DM annihilation spectrum
    # temp
    timeinit = time.time()
    multintptemp = interp1d(gdat.masspartp4dm, gdat.multp4dm, axis=1)(gdat.masspart) / gdat.masspart / gdat.enelscalp4dm
    indxeneltemp = where((gdat.meanenel > amin(gdat.enelscalp4dm * gdat.masspart)) & (gdat.meanenel < amax(gdat.enelscalp4dm * gdat.masspart)))[0]
    multintp = zeros(gdat.numbenel)
    multintp[indxeneltemp] = interp1d(gdat.enelscalp4dm * gdat.masspart, multintptemp)(gdat.meanenel[indxeneltemp])
    timefinl = time.time()
    gdat.listchro[0, gdat.indxswepcntr] = timefinl - timeinit

    # energy density and velocity variance in DM halos
    timeinit = time.time()
    velovarihalo = zeros((gdat.numbmass, gdat.numbreds, gdat.numbrsph))
    gdat.edendmathalo = zeros((gdat.numbmass, gdat.numbreds, gdat.numbrsph))
    for d in range(gdat.numbmass):
        for c in gdat.indxreds:
            gdat.edendmathalonorm = gdat.odenviri * gdat.edendmat[c] * \
                                        gdat.conc[d, c]**3 / 3. / (log(1. + gdat.conc[d, c]) - gdat.conc[d, c] / (1. + gdat.conc[d, c])) # [MeV/cm^3]
            gdat.rsphscal = gdat.rsphviri[d, c] / gdat.conc[d, c] # [kpc]
            gdat.rsphnorm = gdat.rsph[d, c, :] / gdat.rsphscal # [1]
            gdat.edendmathalo[d, c, :] = gdat.edendmathalonorm / gdat.rsphnorm / (1. + gdat.rsphnorm)**2 # [MeV/cm^3]
            edendemcnorm = (2. * gdat.demccons - 3.)**3 / 4. / (5. - 2. * gdat.demccons) * gdat.edendmathalonorm # [MeV/cm^3]
            gdat.rsphdemc = (5. - 2. * gdat.demccons) / (2. * gdat.demccons - 3.) * gdat.rsphscal # [kpc]
            velovarihalo[d, c, :] = 4. * pi * gdat.gravconsredu / 200. * 81. * edendemcnorm**(1. / 3.) * gdat.rsphdemc**2 * (gdat.edendmathalo[d, c, :] \
                * (gdat.rsph[d, c, :] / gdat.rsphdemc)**gdat.demccons)**(2. / 3.) / gdat.velolght**2 * gdat.kprc2cmet**2 # [1] 
    gdat.rvelvarihalo = gdat.demcsigm * velovarihalo
    timefinl = time.time()
    gdat.listchro[1, gdat.indxswepcntr] = timefinl - timeinit
          
    # DM annihilation cross section
    timeinit = time.time()
    gdat.csecvelohalo = zeros((gdat.numbmass, gdat.numbreds, gdat.numbrsph, gdat.numbmpol))
    gdat.csecvelohalo[:, :, :, 0] = gdat.csecvelo
    gdat.csecvelohalo[:, :, :, 1] = gdat.csecvelo * gdat.csecfrac * gdat.rvelvarihalo
    timefinl = time.time()
    gdat.listchro[2, gdat.indxswepcntr] = timefinl - timeinit

    # electron emissivity in DM halo
    timeinit = time.time()
    gdat.emiselechalohost = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbrsph, gdat.numbmpol))
    gdat.emiselechalosubh = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbrsph, gdat.numbmpol))
    gdat.emiselechalotemp = zeros((gdat.numbenel, gdat.numbmass, gdat.numbmass, gdat.numbreds, gdat.numbrsph, gdat.numbmpol))
    gdat.emiselechalo = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbrsph, gdat.numbmpol))
    diffemiselechalodiffrsph = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbrsph, gdat.numbmpol))
    lumielechalo = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbmpol))
    for a in range(gdat.numbenel):
        for d in range(gdat.numbmass):
            for c in gdat.indxreds:
                for m in range(gdat.numbmpol):
                    gdat.emiselechalo[a, d, c, :, m] = gdat.csecvelohalo[d, c, :, m] / gdat.masspart**2 * multintp[a] * gdat.fesc[a, d, c, :] * \
                                                                                                                            gdat.edendmathalo[d, c, :]**2 # [1/s/cm^3/MeV]
                    diffemiselechalodiffrsph[a, d, c, :, m] = 4. * pi * gdat.rsph[d, c, :]**2 * gdat.emiselechalo[a, d, c, :, m] * gdat.kprc2cmet**2 # [1/s/cm/MeV]
                    lumielechalo[a, d, c, m] = trapz(diffemiselechalodiffrsph[a, d, c, :, m] * gdat.kprc2cmet, gdat.rsph[d, c, :]) # [1/s/MeV]
    
    lumielechalo *= gdat.bostfact[None, :, None, None]
    
    gdat.emiselechaloenel = trapz(gdat.emiselechalo, gdat.meanenel, axis=0) # [1/s/cm^3]
    gdat.lumielechaloenel = trapz(lumielechalo, gdat.meanenel, axis=0) # [1/s]
    gdat.diffemiselechalodiffrsphenel = trapz(diffemiselechalodiffrsph, gdat.meanenel, axis=0) # [1/s/cm^3]
    timefinl = time.time()
    gdat.listchro[3, gdat.indxswepcntr] = timefinl - timeinit
    
    # mean DM relative velocity
    timeinit = time.time()
    gdat.diffrvelvarihavgdiffrsph = zeros((gdat.numbmass, gdat.numbreds, gdat.numbrsph))
    gdat.rvelvarihavg = zeros((gdat.numbmass, gdat.numbreds))
    for d in range(gdat.numbmass):
        for c in gdat.indxreds:
            gdat.diffrvelvarihavgdiffrsph[d, c, :] = 3. * gdat.rsph[d, c, :]**2 * gdat.rvelvarihalo[d, c, :] / amax(gdat.rsph[d, c, :])**3
            gdat.rvelvarihavg[d, c] = trapz(gdat.diffrvelvarihavgdiffrsph[d, c, :], gdat.rsph[d, c, :]) # [1]
    timefinl = time.time()
    gdat.listchro[4, gdat.indxswepcntr] = timefinl - timeinit

    # spatially averaged electron emissivity    
    ## smooth component
    # temp -- add p-wave to the smth component
    timeinit = time.time()
    gdat.emiselecsmth = zeros((gdat.numbenel, gdat.numbreds, gdat.numbmpol))
    for c in gdat.indxreds:
        for m in range(gdat.numbmpol):
            gdat.emiselecsmth[:, c, m] = gdat.csecvelo * gdat.edendmat[c]**2 / gdat.masspart**2 * multintp # [1/cm^3/s/MeV]
    gdat.emiselecsmthenel = trapz(gdat.emiselecsmth, gdat.meanenel, axis=0) # [1/cm^3/s]
    ## clumpy component
    diffemiselecclmpdiffmass = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbmpol))
    for a in range(gdat.numbenel):
        for d in range(gdat.numbmass):
            for c in gdat.indxreds:
                for m in range(gdat.numbmpol):
                    diffemiselecclmpdiffmass[a, :, c, m] = gdat.halodist[:, c] * lumielechalo[a, :, c, m] / gdat.kprc2cmet**3 * \
                                                                        (1. + gdat.meanreds[c])**3 # [1/cm^3/s/MeV/Msun]
    gdat.emiselecclmp = trapz(diffemiselecclmpdiffmass, gdat.meanmass, axis=1) # [1/cm^3/s/MeV]
    gdat.emiselecclmpenel = trapz(gdat.emiselecclmp, gdat.meanenel, axis=0) # [1/cm^3/s]
    gdat.diffemiselecclmpdiffmassenel = trapz(diffemiselecclmpdiffmass, gdat.meanenel, axis=0) # [1/cm^3/s/Msun]
    
    ## total
    gdat.emiselec = gdat.emiselecsmth + gdat.emiselecclmp
    gdat.emiselecenel = gdat.emiselecsmthenel + gdat.emiselecclmpenel
    timefinl = time.time()
    gdat.listchro[5, gdat.indxswepcntr] = timefinl - timeinit

    # dark matter velocity variance
    timeinit = time.time()
    diffrvelvariiavgdiffmass = gdat.halodist * gdat.rvelvarihavg * 4. * pi * amax(gdat.rsph, axis=2)**3 / 3. # [1/Msun]
    gdat.rvelvariiavgclmp = trapz(diffrvelvariiavgdiffmass, gdat.meanmass, axis=0) # [1]
    tempiavgsmth = gdat.tempcmbrnunc * (1. + gdat.meanreds) / (1. + 1e9  / (1. + gdat.meanreds) / (1. + ((1. + gdat.meanreds) / 1e9)**2.5))
    gdat.rvelvariiavgsmth = 3. * gdat.boltcons * tempiavgsmth / gdat.masspart
    gdat.rvelvariiavg = gdat.rvelvariiavgsmth + gdat.rvelvariiavgclmp
    timefinl = time.time()
    gdat.listchro[6, gdat.indxswepcntr] = timefinl - timeinit

    # mean electron number density in the IGM
    timeinit = time.time()
    gdat.numbelecenelreds = zeros((gdat.numbenel, gdat.numbreds, gdat.numbmpol))
    for c in gdat.indxreds:
        for m in range(gdat.numbmpol):
            for i in range(gdat.numbenel-1):
                gdat.numbelecenelreds[i, c, m] = trapz(gdat.emiselec[i:gdat.numbenel, c, m], gdat.meanenel[i:], axis=0) / gdat.edotegbl[i, c] # [1/cm^3/MeV]
    gdat.numbelecreds = trapz(gdat.numbelecenelreds, gdat.meanenel, axis=0) # [1/cm^3]
    timefinl = time.time()
    gdat.listchro[7, gdat.indxswepcntr] = timefinl - timeinit

    # photon emissivity
    timeinit = time.time()
    gdat.emisphot = trapz(gdat.specinvc[:, :, :, None] * gdat.numbelecenelreds[None, :, :, :], gdat.meanenel, axis=1) # [1/cm^3/s/MeV]
    gdat.emisphotrydb = interp1d(gdat.meanenph, gdat.emisphot, axis=0)(gdat.enerrydb) # [1/cm^3/s/MeV]
    gdat.emisphotenph = trapz(gdat.emisphot, gdat.meanenph, axis=0) # [1/cm^3/s]
    timefinl = time.time()
    gdat.listchro[8, gdat.indxswepcntr] = timefinl - timeinit

    # photon flux 
    timeinit = time.time()
    difffluxphotdiffreds = zeros((gdat.numbenph, gdat.numbreds, gdat.numbreds, gdat.numbmpol))
    gdat.fluxphotdmat = zeros((gdat.numbenph, gdat.numbreds, gdat.numbmpol))
    for a in range(gdat.numbenph):
        for c in gdat.indxreds:
            gdat.meanenphshft = gdat.meanenph[a] * (1. + gdat.meanreds[c]) / (1. + gdat.meanreds)
            gdat.indxredsshft = where((gdat.minmenph < gdat.meanenphshft) & (gdat.meanenphshft < gdat.maxmenph))[0]
            if gdat.indxredsshft.size == 0:
                continue
            for m in range(gdat.numbmpol):
                gdat.emisphotintp = interp1d(gdat.meanenph, gdat.emisphot[:, c, m])(gdat.meanenphshft[gdat.indxredsshft]) # [1/cm^3/s/MeV]
                difffluxphotdiffreds[a, c, gdat.indxredsshft, m] = gdat.velolght * gdat.timehubb[gdat.indxredsshft] / 4. / pi * \
                   gdat.emisphotintp * exp(-gdat.odep[a, c, gdat.indxredsshft]) * (1. + gdat.meanreds[c])**3 / (1. + gdat.meanreds[gdat.indxredsshft])**4 # [1/cm^2/s/sr/MeV]
                gdat.fluxphotdmat[a, c, m] = trapz(difffluxphotdiffreds[a, c, :, m], gdat.meanreds) # [1/cm^2/s/sr/MeV]  
    gdat.fluxphotrydb = interp1d(gdat.meanenph, gdat.fluxphotdmat, axis=0)(gdat.enerrydb) # [1/cm^2/s/sr/MeV]
    gdat.fluxphotenph = trapz(gdat.fluxphotdmat, gdat.meanenph, axis=0) # [1/cm^2/s/sr]
    timefinl = time.time()
    gdat.listchro[9, gdat.indxswepcntr] = timefinl - timeinit

    # photo-ionization rate
    timeinit = time.time()
    gdat.diffionrdiffenph = 4. * pi * gdat.csecionz[:, None, None] * gdat.fluxphotdmat # [1/s/MeV]
    gdat.ionrdmat = trapz(gdat.diffionrdiffenph, gdat.meanenph, axis=0) # [1/s]
    timefinl = time.time()
    gdat.listchro[10, gdat.indxswepcntr] = timefinl - timeinit
    gdat.indxswepcntr += 1


def plot_halo(gdat):

    listlinestyl = [':', '--', '-']
    listcolr = ['b', 'g', 'r']
    
    # concentration parameter
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for m in range(gdat.numbconcmodltype):
        for c in range(gdat.numbredsplotlate):
            axis.loglog(gdat.meanmass, gdat.conccatl[:, gdat.indxredsplotlate[c], m], ls=listlinestyl[m], color=listcolr[c])
    axis.set_xlabel('$M [M_\odot]$')
    axis.set_ylabel('$c(M,z)$')
    axis.set_ylim([1., None])

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
    listlabl.append(gdat.strgredsplotlate[0])
    listlabl.append(gdat.strgredsplotlate[1])
    listlabl.append(gdat.strgredsplotlate[2])
    
    axis.legend(listline, listlabl)

    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'halo.pdf')
    plt.close()

    # substructure
    #figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    #axis.loglog(rsph, subfs)
    #axis.loglog(gdat.rsph, subbfrsp)           
    #axis.set_xlabel('$r$ [kpc]')
    #axis.set_ylabel('$f_s(r)$')

    # DM energy density and relative velocity variance in halos
    temp = meshgrid(gdat.indxmassplot, gdat.indxredsplotlate, indexing='ij')
    labl = gdat.strgmass + gdat.strgredsplotlate

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    tdpy.util.plot_matr(axis, gdat.rsph[temp[0], temp[1], :], gdat.edendmathalo[temp[0], temp[1], :], labl=labl, loc=3)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlabel('$r$ [kpc]')
    axis.set_ylabel(r'$\rho(r,M,z)$ [MeV/cm$^3$]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'edendmathalo.pdf')
    plt.close()

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    tdpy.util.plot_matr(axis, gdat.rsph[temp[0], temp[1], :], gdat.rvelvarihalo[temp[0], temp[1], :], labl=labl, loc=2)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlabel('$r$ [kpc]')
    axis.set_ylabel('$v^2 (M,z,r)$')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'rvelvarihalo.pdf')
    plt.close()

    # virial and scale radii
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for c in range(gdat.numbredsplotlate):
        axis.loglog(gdat.meanmass, gdat.rsphviri[:, gdat.indxredsplotlate[c]], label=gdat.strgredsplotlate[c])
    axis.set_xlabel('$M [M_\odot]$')
    axis.set_ylabel('$r_{vir}(M,z)$ [kpc]')
    axis.legend(loc=2)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'rsphviri.pdf')
    plt.close()

    # e^-/e^+ emissivity in a halo
    figr, axisgrid = plt.subplots(gdat.numbmassplot, gdat.numbmpol, sharey='row', sharex='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.numbmassplot * gdat.plotsize))
    for d, axisrows in enumerate(axisgrid):
        for p, axis in enumerate(axisrows):
            for c in range(gdat.numbredsplotlate):
                axis.loglog(gdat.rsph[gdat.indxmassplot[d], gdat.indxredsplotlate[c], :], gdat.emiselechaloenel[gdat.indxmassplot[d], gdat.indxredsplotlate[c], :, p], \
                                                                                                                                                label=gdat.strgredsplotlate[c])
            axis.set_xlabel('$r$ [kpc]')
            if d == 0:
                axis.set_title(gdat.strgmpol[p])
            axis.text(0.23, 0.15, gdat.strgmass[d], ha='center', va='center', transform = axis.transAxes)
            if d == gdat.numbmassplot / 2:
                axis.set_ylabel(r'$dN_e/dtdV (M,z,r)$ [1/s/cm$^3$]')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'emiselechaloenel.pdf')
    plt.close()

    # e^-/e^+ differential luminosity of a halo
    figr, axisgrid = plt.subplots(gdat.numbmassplot, gdat.numbmpol, sharey='row', sharex='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.numbmassplot * gdat.plotsize))
    for d, axisrows in enumerate(axisgrid):
        for p, axis in enumerate(axisrows):
            for c in range(gdat.numbredsplotlate):
                axis.loglog(gdat.rsph[gdat.indxmassplot[d], gdat.indxredsplotlate[c], :], gdat.diffemiselechalodiffrsphenel[gdat.indxmassplot[d], \
                                                                                                        gdat.indxredsplotlate[c], :, p], label=gdat.strgredsplotlate[c])
            if d == gdat.numbmassplot - 1:
                axis.set_xlabel('$r$ [kpc]')
            if d == 0:
                axis.set_title(gdat.strgmpol[p])
            axis.text(0.23, 0.12, gdat.strgmass[d], ha='center', va='center', transform = axis.transAxes)
            if d == gdat.numbmassplot / 2 and p == 0:
                axis.set_ylabel(r'$dN_e/dt/dVd\log r$ [1/s/cm$^3$]')
    axis.legend()
    figr.subplots_adjust(hspace=0.1, wspace=0.1)
    plt.savefig(gdat.pathimag + 'diffemiselechalodiffrsphenel.pdf')
    plt.close()

    # e^-/e^+ luminosity of a halo
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='row', sharex='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for p, axis in enumerate(axisrows):
        for c in range(gdat.numbredsplotlate):
            axis.loglog(gdat.meanmass, gdat.lumielechaloenel[:, gdat.indxredsplotlate[c],p], label=gdat.strgredsplotlate[c])
        axis.set_xlabel(r'$M$ $[M_\odot]$')
        axis.set_title(gdat.strgmpol[p])
        if p == 0:
            axis.set_ylabel('$dN_e/dt (M,z)$ [1/s]')
    axis.legend(loc=2)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'lumielechaloenel.pdf')
    plt.close()

    # spatially averaged e^-/e^+ differential emissivity
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='row', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for p, axis in enumerate(axisrows):
        for c in range(gdat.numbredsplotlate):
            if False:
                print 'p'
                print p
                print 'c'
                print c
                print 'reds'
                print gdat.meanreds[gdat.indxredsplotlate[c]]
                print 'gdat.diffemiselecclmpdiffmassenel'
                print gdat.diffemiselecclmpdiffmassenel[:, gdat.indxredsplotlate[c], p] * gdat.meanmass
                print 
            axis.loglog(gdat.meanmass, gdat.meanmass * gdat.diffemiselecclmpdiffmassenel[:, gdat.indxredsplotlate[c], p], label=gdat.strgredsplotlate[c])
        axis.set_xlabel(r'$M$ $[M_\odot]$')
        axis.set_title(gdat.strgmpol[p])
        if p == 0:
            axis.set_ylabel(r'$dN_e/dtdVd\log M$ [1/s/cm$^3$]')
            axis.legend()
        axis.set_ylim([1e-36, 1e-28])
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'diffemiselecclmpdiffmassenel.pdf')
    plt.close()

    ylabel = '$dN_e/dVdt (z)$ [1/cm$^3$/s]'
    listvarb = [gdat.emiselecenel, gdat.emiselecclmpenel, gdat.emiselecsmthenel]
    listlabl = ['Total', 'Clumpy', 'Smooth']
    figr, axisrows = plt.subplots(1, gdat.numbmpol, figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for p, axis in enumerate(axisrows):
        axis.set_xlabel('$z$')
        for g in range(3):
            axis.loglog(gdat.meanreds, listvarb[g].squeeze()[:, p], label=listlabl[g])
        if p == 0:
            axis.set_ylabel(ylabel)
        axis.set_title(gdat.strgmpol[p])
        axis.set_ylim([1e-36, 1e-24])
        axis.loglog(gdat.meanredssfrd, gdat.emiselecsfrd, markersize=10, label='Star Formation', color='black', ls='none', marker='o')
        if p == 0:
            axis.legend(loc=1)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'emiselecenel.pdf')
    plt.close()

    ylabel = '$dN_e/dVdt (z)$ [1/cm$^3$/s]'
    listvarb = [gdat.emiselecenel, gdat.emiselecclmpenel, gdat.emiselecsmthenel]
    figr, axisrows = plt.subplots(1, gdat.numbmpol, figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for p, axis in enumerate(axisrows):
        axis.set_xlabel('$z$')
        for g in range(3):
            axis.loglog(gdat.meanreds, listvarb[g].squeeze()[:, p] / (1. + gdat.meanreds)**3, label=listlabl[g])
        if p == 0:
            axis.set_ylabel(ylabel)
        axis.set_title(gdat.strgmpol[p])
        axis.set_ylim([1e-36, 1e-24])
        axis.loglog(gdat.meanredssfrd, gdat.emiselecsfrd / (1. + gdat.meanredssfrd)**3, markersize=10, label='Star Formation', color='black', ls='none', marker='o')
        if p == 0:
            axis.legend(loc=1)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'emiselecenelcomo.pdf')
    plt.close()

    # Mean DM relative velocity variance in a halo
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for c in range(gdat.numbredsplotlate):
        axis.loglog(gdat.meanmass, sqrt(gdat.rvelvarihavg[:, gdat.indxredsplotlate[c]]), label=gdat.strgredsplotlate[c])
    axis.set_ylabel('$\sigma_v$')
    axis.set_xlabel(r'$M$ $[M_\odot]$')
    axis.legend(loc=2)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'rvelvarihavg.pdf')
    plt.close()

    # Mean DM relative velocity variance
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.set_xlabel('$z$')
    axis.set_ylim([1e-12, 1.])
    axis.loglog(gdat.meanreds, sqrt(gdat.rvelvariiavgclmp))
    axis.loglog(gdat.meanreds, sqrt(gdat.rvelvariiavgsmth))
    axis.loglog(gdat.meanreds, sqrt(gdat.rvelvariiavg))
    axis.set_ylabel('$\sigma_v$')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'rvelvariiavg.pdf')
    plt.close()

    
def plot_elec_flux(gdat):
    
    # e^-/e^+ number density in the IGM
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for p, axis in enumerate(axisrows):
        for c in range(gdat.numbredsplotlate):
            axis.loglog(gdat.meanenel, gdat.meanenel * gdat.numbelecenelreds[:, gdat.indxredsplotlate[c], p], label=gdat.strgredsplotlate[c])
        axis.set_xlabel(r'$E_e$ [GeV]')
        axis.set_title(gdat.strgmpol[p])
        if p == 0:
            axis.set_ylabel(r'$dN_e/dV/d\log E$ [1/cm$^3$]')
    axis.legend(loc=2)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'numbelecenelreds.pdf')
    plt.close()

    # e^-/e^+ number density in the IGM
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for p, axis in enumerate(axisrows):
        axis.loglog(gdat.meanreds, gdat.numbelecreds[:, p])
        axis.set_xlabel('$z$')
        axis.set_title(gdat.strgmpol[p])
        if p == 0:
            axis.set_ylabel('$dN_e/dV (z)$ [1/cm$^3$]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'numbelecreds.pdf')
    plt.close() 
    
        
def plot_invrcomp(gdat):

    listlinestyl = ['-', '--', '-.']
    listcolr = ['b', 'g', 'r']
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for b in range(gdat.numbenelplot):
        for c in range(gdat.numbredsplotlate):
            axis.loglog(gdat.meanenph * 1e6, gdat.meanenph * sum(gdat.specinvccatl[:, gdat.indxenelplot[b], gdat.indxredsplotlate[c], :], 1), ls=listlinestyl[b], c=listcolr[c])
    axis.set_xlabel('$E_\gamma$ [eV]')
    axis.set_ylabel(r'$dN_\gamma/dtd\log E_\gamma$ [1/s]')
    axis.set_ylim([1e-15, 1e-7])

    listlabl = []
    listlabl.append('$E_e$ = ' + tdpy.util.mexp(gdat.meanenel[gdat.indxenelplot[0]] * 1e-3) + ' GeV')
    listlabl.append('$E_e$ = ' + tdpy.util.mexp(gdat.meanenel[gdat.indxenelplot[1]] * 1e-3) + ' GeV')
    listlabl.append('$E_e$ = ' + tdpy.util.mexp(gdat.meanenel[gdat.indxenelplot[2]] * 1e-3) + ' GeV')
    listlabl.extend(gdat.strgredsplotlate)
    
    listline = []
    for k in arange(3):
        listline.append(plt.Line2D((0,1),(0,0), color='black', ls=listlinestyl[k]))
    for k in range(3):
        listline.append(plt.Line2D((0,1),(0,0), color=listcolr[k]))
    axis.legend(listline, listlabl, loc=9, ncol=2) 
    
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'specinvccatl.pdf')
    plt.close() 
    

def plot_igma(gdat):
      
    extt = [gdat.minmreds, gdat.maxmreds, gdat.minmcden, gdat.maxmcden]
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    zdat = gdat.meanreds[None, :] * gdat.meancden[:, None] * gdat.diffnabsdiffcdendiffreds
    imag = axis.imshow(zdat, extent=extt, cmap='Reds', norm=mpl.colors.LogNorm(), aspect='auto')
    for k in range(gdat.cdenbrek.size):
        axis.axhline(gdat.cdenbrek[k], color='grey', ls='--')
    for k in range(gdat.meanredsbrek.size):
        axis.axvline(gdat.meanredsbrek[k], color='grey', ls='--')
    axis.set_ylabel('$N_{HI}$ [cm$^{-2}$]')
    axis.set_xlabel('$z$')
    axis.set_title('$d^2N_{abs}/d\log N_{HI}/d\log z$')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.colorbar(imag, ax=axis, fraction=0.04)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'diffnabsdiffcdendiffreds.pdf')
    plt.close() 
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    zdat = gdat.meanreds[None, :] * gdat.meancden[:, None] * gdat.diffoptddiffcdendiffredsrydb
    zdat[where(zdat > 1.)] = 1.
    imag = axis.imshow(zdat, extent=extt, cmap='Reds', norm=mpl.colors.LogNorm(), aspect='auto')
    for k in range(gdat.cdenbrek.size):
        axis.axhline(gdat.cdenbrek[k], color='grey', ls='--')
    for k in range(gdat.meanredsbrek.size):
        axis.axvline(gdat.meanredsbrek[k], color='grey', ls='--')
    axis.set_ylabel('$N_{HI}$ [cm$^{-2}$]')
    axis.set_xlabel('$z$')
    axis.set_title(r'$d^2\tau_{abs}/d\log N_{HI}/d\log z$')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.colorbar(imag, ax=axis, fraction=0.04)
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'diffoptddiffcdendiffredsrydb.pdf')
    plt.close() 
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.loglog(gdat.meanreds, gdat.meanreds * gdat.diffoptddiffredsrydb)
    for k in range(gdat.meanredsbrek.size):
        axis.axvline(gdat.meanredsbrek[k], color='grey', ls='--')
    axis.set_xlabel('$z$')
    axis.set_ylabel(r'$d\tau / d\log z$')
    axis.set_ylim([None, 1.])
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'diffoptddiffredsrydb.pdf')
    plt.close() 
        
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for c in range(gdat.numbredsplotlate):
        axis.loglog(gdat.meanreds[gdat.indxredsplotlate[c]:], gdat.odeprydb[gdat.indxredsplotlate[c], gdat.indxredsplotlate[c]:], label=gdat.strgredsplotlate[c])
    for k in range(gdat.meanredsbrek.size):
        axis.axvline(gdat.meanredsbrek[k], color='grey', ls='--')
    axis.set_xlabel('$z$')
    axis.set_ylabel(r'$\tau$')
    axis.set_ylim([None, 1.])
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'odep.pdf')
    plt.close() 
    
    
def plot_fluxphot(gdat):
    
    # differential photon emissivity 
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for i, axis in enumerate(axisrows):
        for c in range(gdat.numbredsplotlate):
            plot = axis.loglog(gdat.meanenph * 1e6, gdat.meanenph * gdat.emisphot[:, gdat.indxredsplotlate[c], i], label=gdat.strgredsplotlate[c])
        axis.set_xlabel('$E_\gamma$ [eV]')
        axis.set_ylabel(r'$dN_\gamma/dtdVd\log E_\gamma$ [1/s]')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'emisphot.pdf')
    plt.close() 
    
    # photon emissivity 
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for i, axis in enumerate(axisrows):
        plot = axis.loglog(gdat.meanreds, gdat.enerrydb * gdat.emisphotrydb[:, i])
        axis.set_xlabel('$z$')
        axis.set_ylabel(r'$dN_\gamma/dtdV$ [1/cm$^3$/s]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'emisphotrydb.pdf')
    plt.close() 
    
    # differential photon flux 
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for m, axis in enumerate(axisrows):
        for c in range(gdat.numbredsplotlate):
            plot = axis.loglog(gdat.meanenph * 1e6, gdat.meanenph * gdat.fluxphotdmat[:, gdat.indxredsplotlate[c], m], label=gdat.strgredsplotlate[c])
        axis.set_xlabel('$E_\gamma$ [eV]')
        axis.set_ylabel(r'$dN_\gamma/dtdAd\log E_\gamma$ [1/cm$^2$/s/sr]')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'fluxphot.pdf')
    plt.close() 
    
    # differential photon flux at 1 Ryd
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for m, axis in enumerate(axisrows):
        plot = axis.loglog(gdat.meanreds, gdat.enerrydb * gdat.fluxphotrydb[:, m])
        axis.set_xlabel('$z$')
        axis.set_ylabel(r'$dN_\gamma/dtdAd\log E_\gamma$ [1/cm$^2$/s/sr]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'fluxphotrydb.pdf')
    plt.close() 
    
    # integrated photon flux 
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for m, axis in enumerate(axisrows):
        plot = axis.loglog(gdat.meanreds, gdat.fluxphotenph[:, m])
        axis.set_xlabel('$z$')
        axis.set_ylabel(r'$dN_\gamma/dtdA (z)$ [1/cm$^2$/s/sr]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'fluxphotenph.pdf')
    plt.close() 

        
def plot_ionr(gdat):

    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for i, axis in enumerate(axisrows):
        for c in arange(len(gdat.indxredsplotlate)):
            plot = axis.loglog(gdat.meanenph * 1e6, gdat.diffionrdiffenph[:, gdat.indxredsplotlate[c], i], label=gdat.strgredsplotlate[c])
        axis.set_xlabel('$E_\gamma [eV]$')
        axis.set_ylabel(r'$d\Gamma/d\log E_\gamma$ [1/s/H]')
        axis.axvline(gdat.enerrydb * 1e6, ls='--', color='grey')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'diffionrdiffenph.pdf')
    plt.close() 
        
    liststrg = ['ionrbolt.csv', 'ionrbeck.csv', 'ionrfauc.csv']
    listlabl = ['Bolton & Haehnelt (2007)', 'Becker et al. (2007)', 'Faucher-Giguere (2008)']
    listmrkr = ['o', 'D', 'x']
    
    figr, axisrows = plt.subplots(1, gdat.numbmpol, sharey='all', figsize=(gdat.numbmpol * gdat.plotsize, gdat.plotsize))
    for i, axis in enumerate(axisrows):
        plot = axis.loglog(gdat.meanreds, gdat.ionrdmat[:,i])
        axis.set_xlabel('$z$')
        axis.set_ylabel(r'$\Gamma$ [1/s/H]')
        
        for k, strg in enumerate(liststrg):
            path = gdat.pathdata + strg
            data = loadtxt(path)
            ndata = data.shape[0] / 3
            yerr = zeros((2, ndata))
            xdat = data[:ndata, 0]
            ydat = data[:ndata, 1] * 1e-12
            yerr[0, :] = data[ndata:2*ndata, 1]
            yerr[1, :] = data[2*ndata:3*ndata, 1]
            yerr = abs(yerr - ydat)
            axis.errorbar(xdat, ydat, yerr=yerr, ls='', marker=listmrkr[k])

    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'ionr.pdf')
    plt.close() 
        

def retr_fluxphotexpr(gdat):
 
    path = gdat.pathdata + 'xray_background.dat'
    tabl = loadtxt(path)
    
    gdat.meanenphexpr = tabl[:,0]
    gdat.meanenphexpr *= 1e-6 # [MeV]
    gdat.fluxphotexpr = tabl[:, 1]
    gdat.fluxphotexpr *= 1e-3 / gdat.meanenphexpr**2 # [1/cm^2/s/sr/MeV]
    
    gdat.indxenphexpr = where((amin(gdat.meanenphexpr) < gdat.meanenph) & (gdat.meanenph < amax(gdat.meanenphexpr)))[0]
    
    gdat.fluxphotexpr = interp1d(gdat.meanenphexpr, gdat.fluxphotexpr)(gdat.meanenph[gdat.indxenphexpr])
    gdat.fluxphotexprvari = (gdat.fluxphotexpr * 1.)**2
    

def retr_datapara(gdat):
    
    gdat.numbpara = 3

    datapara = tdpy.util.gdatstrt()

    datapara.indx = dict()
    datapara.minm = zeros(gdat.numbpara)
    datapara.maxm = zeros(gdat.numbpara)
    datapara.name = empty(gdat.numbpara, dtype=object)
    datapara.scal = empty(gdat.numbpara, dtype=object)
    datapara.labl = empty(gdat.numbpara, dtype=object)
    datapara.unit = empty(gdat.numbpara, dtype=object)
    datapara.true = empty(gdat.numbpara, dtype=object)
    datapara.vari = zeros(gdat.numbpara)

    datapara.indx['csecvelo'] = 0
    datapara.name[0] = 'csecvelo'
    datapara.minm[0] = 3e-26
    datapara.maxm[0] = 3e-20
    datapara.scal[0] = 'logt'
    datapara.labl[0] = '$a$'
    datapara.unit[0] = '[cm$^3$/s]'
    datapara.vari[0] = 3e-1
    datapara.true[0] = None
    
    datapara.indx['csecfrac'] = 1
    datapara.name[1] = 'csecfrac'
    datapara.minm[1] = 1e-2
    datapara.maxm[1] = 1e10
    datapara.scal[1] = 'logt'
    datapara.labl[1] = '$b/a$'
    datapara.unit[1] = ''
    datapara.vari[1] = 3e-1
    datapara.true[1] = None
    
    datapara.indx['masspart'] = 2
    datapara.name[2] = 'masspart'
    datapara.minm[2] = 1e4
    datapara.maxm[2] = 1e6
    datapara.scal[2] = 'logt'
    datapara.labl[2] = '$M$'
    datapara.unit[2] = '[MeV]'
    datapara.vari[2] = 3e-1
    datapara.true[2] = None

    if False:
        datapara.indx['dmatslop'] = 3
        datapara.name[3] = 'dmatslop'
        datapara.minm[3] = 0.8
        datapara.maxm[3] = 1.5
        datapara.scal[3] = 'self'
        datapara.labl[3] = r'$\gamma$'
        datapara.unit[3] = ''
        datapara.vari[3] = 3e-1

    datapara.strg = datapara.labl + ' ' + datapara.unit
    
    return datapara


def retr_fluxphotdmatintp(gdat):
    
    fluxphotdmatintp = gdat.fluxphotdmat[:, :, 0] * gdat.csecvelo / gdat.csecvelopivt + gdat.fluxphotdmat[:, :, 1] * gdat.csecfrac / gdat.csecfracpivt
    
    if False:
        
        print 'retr_fluxphotdmatintp'
        print 'csecvelo'
        print gdat.csecvelo
        print 'csecvelopivt'
        print gdat.csecvelopivt
        print 'csecfrac'
        print gdat.csecfrac
        print 'csecfracpivt'
        print gdat.csecfracpivt
        print 'fluxphotdmat[:, :, 0]'
        print gdat.fluxphotdmat[:, :, 0]
        print 'fluxphotdmat[:, :, 1]'
        print gdat.fluxphotdmat[:, :, 1]
        print 'fluxphotdmatintp[:, 0]'
        print fluxphotdmatintp[:, 0]
        
    return fluxphotdmatintp


def plot_sfrd(gdat):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    
    axis.plot(gdat.meanredssfrd, gdat.sfrd)
    axis.set_yscale('log')
    
    axis.set_xlabel('$z$')
    axis.set_ylabel(r'SFRD [erg/Mpc$^3$/yr]')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'sfrd.pdf')
    plt.close()


def plot_hm12(gdat, fluxphotdmat=None, listfluxphotdmat=None):
    
    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    
    xdat = gdat.meanenph[gdat.indxenphexpr] * 1e6
    ydatlowr = gdat.meanenph[gdat.indxenphexpr] * (gdat.fluxphotexpr - sqrt(gdat.fluxphotexprvari))
    ydatuppr = gdat.meanenph[gdat.indxenphexpr] * (gdat.fluxphotexpr + sqrt(gdat.fluxphotexprvari))
    axis.fill_between(xdat, ydatuppr, ydatlowr, color='lightblue')
    
    axis.loglog(gdat.meanenph[gdat.indxenphexpr] * 1e6, gdat.meanenph[gdat.indxenphexpr] * gdat.fluxphotexpr, label='ROSAT')
    
    liststrg = ['moretotl', 'moreothr', 'more']
    listlabl = ['Total XRB, Moretti et al. (2009)', 'Unresolved XRB, Worsley et al. (2006)', 'Unresolved XRB, Moretti et al. (2012)']
    listcolr = ['black', 'yellow', 'green']
    for k, strg in enumerate(liststrg):
        path = gdat.pathdata + strg + '.csv'
        datamore = loadtxt(path)
        gdat.meanenphmore = datamore[:, 0] * 1e-3 # [MeV]
        fluxphotmore = gdat.ergs2mgev * (180. / pi)**2 * datamore[:, 1] / gdat.meanenphmore**2 # [1/cm^2/s/sr/MeV]
        #axis.fill_between(gdat.meanenphmore * 1e6, 2. * gdat.meanenphmore * fluxphotmore, 0.5 * gdat.meanenphmore * fluxphotmore, color='lightgreen')
        axis.loglog(gdat.meanenphmore * 1e6, gdat.meanenphmore * fluxphotmore, label=listlabl[k], color=listcolr[k])

    listcolr = ['b', 'g', 'r']
    for c in range(gdat.numbredsplotlate):
        axis.loglog(gdat.meanenph * 1e6, gdat.meanenph * gdat.fluxphothm12[:, gdat.indxredsplotlate[c]], label='Haardt & Madau (2012), ' + gdat.strgredsplotlate[c])
        if fluxphotdmat != None:
            axis.loglog(gdat.meanenph * 1e6, gdat.meanenph * gdat.fluxphotdmat[:, gdat.indxredsplotlate[c]], label='DM, ' + gdat.strgredsplotlate[c])
        if listfluxphotdmat != None:
            tdpy.util.plot_braz(axis, gdat.meanenph * 1e6, gdat.meanenph[None, :] * listfluxphotdmat[:, :, gdat.indxredsplotlate[c]], \
                    lcol=listcolr[c], alpha=0.5, dcol=listcolr[c], mcol='black')
    axis.set_xlabel(r'$E_\gamma$ [eV]')
    axis.set_ylabel(r'$EdN/dE$ [1/cm$^2$/s/sr]')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'fluxphothm12.pdf')
    plt.close()


def init( \
         datatype='inpt', \
         datalabl='XMM-Newton', \
         numbswep=100000, \
         verbtype=1, \
         saveflux=True, \
         makeplot=True, \
         propmodltype='effi', \
         concmodltype='duff', \
         subsmodltype='smth', \
         igmamodltype='clum'
        ):

    # global object
    gdat = tdpy.util.gdatstrt()
    gdat.saveflux = saveflux
    gdat.makeplot = makeplot
    gdat.propmodltype = propmodltype
    gdat.concmodltype = concmodltype
    gdat.subsmodltype = subsmodltype
    gdat.igmamodltype = igmamodltype
    
    datapara = retr_datapara(gdat)
    gdat.numbpara = len(datapara.name)
    gdat.indxpara = arange(gdat.numbpara)

    # plotting
    gdat.plotsize = 6
    
    gdat.strgmpol = ['s-wave', 'p-wave']
    
    gdat.anch = 'b'
    gdat.demcsigm = 6.
    
    # constants
    ## cosmological constants

    gravlght = 1.19e-34 # [cm^3/MeV/s^2]
    gdat.omegbmat = 0.049 # baryonic matter abundance today
    gdat.omegdmat = 0.26 # dark matter abundance today
    gdat.omegradi = 4.8e-5 # radiation abundance today
    gdat.omegdene = 0.69 # dark energy abundance today
    gdat.stdvodennorm = 0.83 # rms density fluctuation in spheres of radius 8/h Mpc
    gdat.psecindx = 0.96 # spectral index of the primordial power spectrum
    gdat.hubbcons = 0.704 # reduced Hubble constant
    gdat.tempcmbrnunc = 2.725 # CMB temperature today [K]
    gdat.shtrwgth = 0.707 # Sheth-Tormen
    gdat.shtrnorm = 0.3222 # Sheth-Tormen
    gdat.shtrindx = 0.3 # Sheth-Tormen
    
    gdat.omegmatt = gdat.omegbmat + gdat.omegdmat
    
    gdat.masselec = 0.511 # electron gdat.meanmass [MeV]
    gdat.velolght = 2.998e10 # speed of light [cm/s]
    gdat.strtcons = 7.297e-3 # fine structure constant
    
    gdat.odencoll = 1.686 # linear overdensity at collapse
    gdat.odenviri = 18. * pi**2 # overdensity at virialization
    
    gdat.radisolr = 8.5 # radial distance from the GC to the Sun [kpc]
    gdat.edencritnunc = 5.3e-3 # critical density of the Universe [MeV/cm^3]

    gdat.csecthom = 6.65e-25 # [cm^2] Thompson cross section 
    gdat.csecionzrydb = 6.3e-18 # [cm^2] neutral Hydrogen photon-ionization cross section at 13.6 eV 
    gdat.enerrydb = 13.5984e-6 # Rydberg energy [MeV]

    gdat.plnkcons = 4.136e-21 # Planck constant [MeV s]
    gdat.boltcons = 8.6173e-11 # Boltzmann constant [MeV/K]

    gdat.massmilk = 1e12 # gdat.meanmass of the Milky Way [Msun]

    gdat.ergs2mgev = 6.241509e5 # conversion factor from erg to MeV
    gdat.myrs2secd = 3.154e13 # million year per second
    gdat.solm2mgev = 1.115e60 # Solar gdat.meanmass in MeVc^2
    gdat.kprc2cmet = 3.086e21 # conversion factor from kpc to cm
    gdat.magffac = 4.966835e-8 # [(MeV/cm^3)/(muG^2/mu0)] 
    
    gdat.demccons = 35. / 18. # Dehnen McLaughlin profile index
    gdat.sigm = 6. # ratio of mean squared relative velocity to 1-particle velocity variance
    gdat.gravconsredu = 1.19e-34 # [cm^3/MeV/s^2]
    gdat.cmet2angs = 1e8
    
    # axis
    # temp
    gdat.numbradi = 50
    gdat.numbrsph = 50
    gdat.numbzaxi = 50
    gdat.numbenph = 50
    gdat.numbenel = 50
    gdat.numbenpi = 50
    gdat.numbcden = 50
    gdat.numbreds = 50
    gdat.numbmass = 50
    gdat.numbwnum = 500
  
    gdat.indxradi = arange(gdat.numbradi)
    gdat.indxrsph = arange(gdat.numbrsph)
    gdat.indxzaxi = arange(gdat.numbzaxi)
    gdat.indxenph = arange(gdat.numbenph)
    gdat.indxenel = arange(gdat.numbenel)
    gdat.indxenpi = arange(gdat.numbenpi)
    gdat.indxcden = arange(gdat.numbcden)
    gdat.indxreds = arange(gdat.numbreds)
    gdat.indxmass = arange(gdat.numbmass)
    gdat.indxwnum = arange(gdat.numbwnum)

    gdat.minmenph = 1e-6 # [MeV]
    gdat.maxmenph = 1e-1 # [MeV]
    gdat.meanenph = logspace(log10(gdat.minmenph), log10(gdat.maxmenph), gdat.numbenph)
    gdat.indxenphplot = array([0, gdat.numbenph / 2, gdat.numbenph - 1])

    gdat.minmenel = 5e1 # [MeV]
    gdat.maxmenel = 1e5 # [MeV]
    gdat.meanenel = logspace(log10(gdat.minmenel), log10(gdat.maxmenel), gdat.numbenel)
    gdat.diffenel = gdat.meanenel[1:] - gdat.meanenel[:-1]
    gdat.indxenelplot = array([0, gdat.numbenel / 4, gdat.numbenel / 2])
    gdat.numbenelplot = gdat.indxenelplot.size

    gdat.minmreds = 1e-3
    gdat.maxmreds = 1e2
    gdat.meanreds = logspace(log10(gdat.minmreds), log10(gdat.maxmreds), gdat.numbreds)
    
    gdat.redsplotprox = array([0., 1e3, 1e6])
    gdat.numbredsplot = gdat.redsplotprox.size
    gdat.indxredsplot = empty(gdat.numbredsplot, dtype=int)
    for n, reds in enumerate(gdat.redsplotprox):
        gdat.indxredsplot[n] = argmin(fabs(gdat.meanreds - reds))
        
    gdat.meanredsplotlateprox = array([0.1, 2., 6.])
    gdat.numbredsplotlate = len(gdat.meanredsplotlateprox)
    gdat.indxredsplotlate = []
    gdat.strgredsplotlate = []
    for k in range(gdat.numbredsplotlate):
        gdat.indxredsplotlate.append(argmin(abs(gdat.meanreds - gdat.meanredsplotlateprox[k])))
        gdat.strgredsplotlate.append('$z = %.3g$' % gdat.meanredsplotlateprox[k])
        
    minmmass = 1e8 # [Solar Mass]
    maxmmass = 1e17 # [Solar Mass]
    gdat.meanmassprim = logspace(log10(minmmass), log10(maxmmass), gdat.numbmass + 1)
    gdat.meanmass = gdat.meanmassprim[:-1]
    
    gdat.meanmassprox = [1e10, 1e12, 1e15]
    gdat.numbmassplot = len(gdat.meanmassprox)
    gdat.indxmassplot = []
    gdat.strgmass = []
    for d in range(gdat.numbmassplot):
        gdat.indxmassplot.append(argmin(fabs(gdat.meanmass - gdat.meanmassprox[d])))
        gdat.strgmass.append('$M$ = ' + tdpy.util.mexp(gdat.meanmassprox[d]) + r' $M_\odot$')
       
    gdat.minmcden = 10**11. # [1/cm^2]
    gdat.maxmcden = 10**22. # [1/cm^2]
    gdat.meancden = logspace(log10(gdat.minmcden), log10(gdat.maxmcden), gdat.numbcden)

    gdat.minmenpi = 1e-12 # [MeV]
    gdat.maxmenpi = 1e-5 # [MeV]
    gdat.meanenpi = logspace(log10(gdat.minmenpi), log10(gdat.maxmenpi), gdat.numbenpi)

    # wavenumber axis.s
    gdat.minmwnum = 1e-4
    gdat.maxmwnum = 1e4
    gdat.meanwnum = logspace(log10(gdat.minmwnum), log10(gdat.maxmwnum), gdat.numbwnum)
    
    gdat.edenbmat = gdat.omegbmat * gdat.edencritnunc * (1. + gdat.meanreds)**3
    gdat.edendmat = gdat.omegdmat * gdat.edencritnunc * (1. + gdat.meanreds)**3
    gdat.edenmatt = gdat.omegmatt * gdat.edencritnunc * (1. + gdat.meanreds)**3
    gdat.edenradi = gdat.omegradi * gdat.edencritnunc * (1. + gdat.meanreds)**4
    gdat.edendene = gdat.omegdene * gdat.edencritnunc
    
    gdat.edendmatnunc = gdat.omegdmat * gdat.edencritnunc
    
    gdat.rsphviri = (3. * gdat.meanmass[:, None] * gdat.solm2mgev / 4. / pi / gdat.odenviri / gdat.edendmat / gdat.kprc2cmet**3)**(1. / 3.) # [kpc]
    minmrsph = gdat.rsphviri * 1e-3
    maxmrsph = gdat.rsphviri * 1e1
    gdat.rsph = zeros((gdat.numbmass, gdat.numbreds, gdat.numbrsph))
    for c in gdat.indxreds:
        for d in range(gdat.numbmass):             
            gdat.rsph[d, c, :] = logspace(log10(minmrsph[d,c]), log10(maxmrsph[d,c]), gdat.numbrsph) # [kpc]
    
    frph = gdat.meanenph / gdat.plnkcons # [Hz]
    gdat.csecionz = gdat.csecionzrydb * (gdat.enerrydb / gdat.meanenph)**3
    
    gdat.timehubbnunc = 1. / (100. * gdat.hubbcons / (gdat.kprc2cmet / 1e2))
    gdat.funchubb = sqrt(gdat.omegdene + gdat.omegdmat * (1. + gdat.meanreds**3) + gdat.omegradi * (1. + gdat.meanreds**4))
    gdat.timehubb = gdat.timehubbnunc / gdat.funchubb
    
    gdat.listchro = zeros((11, numbswep + 2))
    gdat.indxswepcntr = 0

    # get the time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # run tag
    gdat.rtag = '%s_%s_%s_%s_%s_%s' % (strgtimestmp, gdat.propmodltype, gdat.concmodltype, gdat.subsmodltype, gdat.igmamodltype, gdat.anch)
    
    # paths
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'phot_ionz/', 'phot_ionz/', gdat.rtag)
    
    # plot annihilation spectrum
    if gdat.makeplot:
        anchlabl = ['$e^-e^+$', r'$\mu\bar{\mu}$', r'$\tau^-\tau^+$', r'$b\bar{b}$']
        listmasspartplot = array([1e4, 1e5, 1e6])
        listlablmasspartplot = ['10 GeV', '100 GeV', '1 TeV']
        listanchplot = ['e', 'mu', 'tau', 'b']
        listcolr = ['b', 'g', 'r', 'm']
        listline = ['-', '--', '-.']

        figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
        for k, anch in enumerate(listanchplot):
            multp4dm, enelscalp4dm, gdat.masspartp4dm = tdpy.util.retr_p4dm_spec(anch)
            for n, masspart in enumerate(listmasspartplot):
                indxtemp = argmin(abs(gdat.masspartp4dm - listmasspartplot[n]))
                axis.loglog(enelscalp4dm * masspart * 1e-3, multp4dm[:, indxtemp], color=listcolr[k], ls=listline[n], label=anch + ', ' + listlablmasspartplot[n])
        axis.set_xlim([0.05, 1e3]) 
        axis.set_ylim([1e-1, 5e2])
        axis.set_xlabel('$E$ [GeV]')
        axis.set_ylabel(r'$dN/d\log E$')
        axis.legend(ncol=4, loc=2)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'multp4dm.pdf')
        plt.close()
    
    gdat.multp4dm, gdat.enelscalp4dm, gdat.masspartp4dm = tdpy.util.retr_p4dm_spec(gdat.anch)

    #if gdat.makeplot:
        #plot_edot
        #plot_galprop

    # Haardt Madau 2012 quasar and galaxy background model
    retr_fluxphothm12(gdat)
    
    # experimental background
    retr_fluxphotexpr(gdat)
    
    if gdat.makeplot:
        plot_hm12(gdat)
    
    # CMB energy density
    gdat.edencmbr = retr_edencmbr(gdat) # [1/cm^3/MeV]
    gdat.edencmbrenpi = trapz(gdat.meanenpi[:, None] * gdat.edencmbr, gdat.meanenpi, axis=0) # [MeV/cm^3]
    
    # EGBL energy density
    gdat.edenegbl = retr_edenegbl(gdat) # [1/cm^3/MeV]
    gdat.edenegblenpi = trapz(gdat.meanenpi[:, None] * gdat.edenegbl, gdat.meanenpi, axis=0) # [MeV/cm^3]
    
    # Energy loss rate on EGBL
    # temp
    gdat.edotegbl = 4. * gdat.csecthom * gdat.velolght / 3. / gdat.masselec**2 * gdat.edencmbrenpi[None, :] * gdat.meanenel[:, None]**2 # [MeV/s]
 
    if gdat.makeplot:
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        for c in range(gdat.numbredsplotlate):
            plot = axis.loglog(gdat.meanenpi * 1e6, gdat.meanenpi * gdat.edencmbr[:, gdat.indxredsplotlate[c]], label=gdat.strgredsplotlate[c])
            plot_ = axis.loglog(gdat.meanenpi * 1e6, gdat.meanenpi * gdat.edenegbl[:, gdat.indxredsplotlate[c]], '--', color=plot[0].get_color())
        axis.set_xlabel('$E_\gamma$ [eV]')
        axis.set_ylabel(r'$dN_\gamma/dVd\log E_\gamma$ [1/cm$^3$]')
        axis.set_ylim([1e-1, 1e6])
        leg = axis.legend(loc=2)
        figr.tight_layout()
        plt.savefig(gdat.pathimag + 'edencmbr.pdf')
        plt.close()

    # halo gdat.meanmass function
    retr_hmfn(gdat)
    
    gdat.numbmpol = 2

    # halo concentration model
    gdat.numbconcmodltype = 3
    gdat.conccatl = zeros((gdat.numbmass, gdat.numbreds, gdat.numbconcmodltype))

    ## Sanchez-Conde & Prada 2014
    cona = [37.5153, -1.5093, 1.636e-2, 3.66e-4, -2.89237e-5, 5.32e-7]
    for k in range(6):
        gdat.conccatl[:, :, 0] += cona[k] * log(gdat.meanmass[:, None] / gdat.hubbcons)**k * gdat.funchubb[None, :]**(-2./3.)
    if gdat.concmodltype == 'sanc':
        gdat.conc = gdat.conccatl[:, :, 0]

    ## Duffy et al. 2007
    gdat.concmassslop = -0.081
    gdat.concredsslop = -0.71
    gdat.concmasspivt = 2e12 / gdat.hubbcons
    gdat.concnorm = 7.8
    gdat.conccatl[:,:,1] = gdat.concnorm * outer((gdat.meanmass / gdat.concmasspivt)**gdat.concmassslop, (1. + gdat.meanreds)**gdat.concredsslop)
    if gdat.concmodltype == 'duff':
        gdat.conc = gdat.conccatl[:,:,1]

    ## Comerford & Natarajan 2007
    gdat.concmassslop = -0.15
    gdat.concredsslop = -1.
    gdat.concmasspivt = 1.3e13 / gdat.hubbcons
    gdat.concnorm = 14.5
    gdat.conccatl[:, :, 2] = gdat.concnorm * (gdat.meanmass[:, None] / gdat.concmasspivt)**gdat.concmassslop * (1. + gdat.meanreds[None, :])**gdat.concredsslop
    if gdat.concmodltype == 'come':
        gdat.conc = gdat.conccatl[:,:,2]

    # substructure model
    gdat.bostfact = empty(gdat.numbmass)
    if gdat.subsmodltype == 'smth':
        gdat.bostfact[:] = 1.
    if gdat.subsmodltype == 'clmp':
        gdat.bostfact[:] = 1.6e-3 * mass**.39
    
    # cosmic ray escape model
    uisrf = 6e-7 # [MeV/cm^3]
    magfd = 5. # [muG]
    gasdn = 1. # [1/cm^3]

    # cosmic raw propagation model
    gdat.fesc = zeros((gdat.numbenel, gdat.numbmass, gdat.numbreds, gdat.numbrsph))
    if gdat.propmodltype == 'effi':
        gdat.fesc[:] = 1.
    else:
        if gdat.propmodltype == 'minm':
            magfd *= 0.2 # [muG]
            gasdn *= 0.2 # [1/cm^3]
            uisrf *= 0.2 # [MeV/cm^3]
        if gdat.propmodltype == 'maxm':
            magfd *= 5. # [muG]
            gasdn *= 5. # [1/cm^3]
            uisrf *= 5. # [MeV/cm^3]

        e0 = 1. # [GeV]
        if gdat.propmodltype == 'minm':
            k0 = 5.074e-17 # [kpc^2/s]
            delta = 0.85
            lz = 1. # [kpc]
        if gdat.propmodltype == 'medi':
            k0 = 3.551e-16 # [kpc^2/s]
            delta = 0.7
            lz = 4. # [kpc]
        if gdat.propmodltype == 'maxm':
            k0 = 2.426e-15 # [kpc^2/s]
            delta = 0.46
            lz = 8. # [kpc]
        
        gdat.fesc = retr_fesc(gdat)
        
    # IGM effective optical depth model
    gdat.cdenbrek = array([10**15., 10**17.5, 10**19., 10**20.3])
    gdat.meanredsbrek = array([1.56, 5.5])
    gdat.odep = zeros((gdat.numbenph, gdat.numbreds, gdat.numbreds))
    gdat.odeprydb = zeros((gdat.numbreds, gdat.numbreds))
    if gdat.igmamodltype == 'clum':
        
        gdat.diffnabsdiffcdendiffreds = zeros((gdat.numbcden, gdat.numbreds))

        # lower Lyman - alpha forest
        indxcdentemp = where(gdat.meancden < 10**15.)[0]

        gdat.indxredstemp = where((gdat.meanreds > 1.56) & (gdat.meanreds < 5.5))[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**7.079 * gdat.meancden[indxcdentemp, None]**(-1.5) * (1. + gdat.meanreds[None, gdat.indxredstemp])**3

        gdat.indxredstemp = where(gdat.meanreds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**8.238 * gdat.meancden[indxcdentemp, None]**(-1.5) * (1. + gdat.meanreds[None, gdat.indxredstemp])**0.16

        gdat.indxredstemp = where(gdat.meanreds > 5.5)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**1.470 * gdat.meancden[indxcdentemp, None]**(-1.5) * (1. + gdat.meanreds[None, gdat.indxredstemp])**9.9

        # upper Lyman - alpha forest
        indxcdentemp = where((gdat.meancden > 10**15.) & (gdat.meancden < 10**17.5))[0]
        
        gdat.indxredstemp = where((gdat.meanreds > 1.56) & (gdat.meanreds < 5.5))[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**14.58 * gdat.meancden[indxcdentemp, None]**(-2.) * (1. + gdat.meanreds[None, gdat.indxredstemp])**3.

        gdat.indxredstemp = where(gdat.meanreds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**15.74 * gdat.meancden[indxcdentemp, None]**(-2.) * (1. + gdat.meanreds[None, gdat.indxredstemp])**0.16

        gdat.indxredstemp = where(gdat.meanreds > 5.5)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**8.97 * gdat.meancden[indxcdentemp, None]**(-2.) * (1. + gdat.meanreds[None, gdat.indxredstemp])**9.9

        # Super Lyman limit systems
        indxcdentemp = where((gdat.meancden > 10**19.) & (gdat.meancden < 10**20.3))[0]

        gdat.indxredstemp = where(gdat.meanreds > 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**(-0.347) * gdat.meancden[indxcdentemp, None]**(-1.05) * (1. + gdat.meanreds[None, gdat.indxredstemp])**1.27

        gdat.indxredstemp = where(gdat.meanreds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**(0.107) * gdat.meancden[indxcdentemp, None]**(-1.05) * (1. + gdat.meanreds[None, gdat.indxredstemp])**0.16

        # Damped Lyman alpha systems
        indxcdentemp = where(gdat.meancden > 10**20.3)[0]

        gdat.indxredstemp = where(gdat.meanreds > 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**18.94 * gdat.meancden[indxcdentemp, None]**(-2.) * (1. + gdat.meanreds[None, gdat.indxredstemp])**1.27

        gdat.indxredstemp = where(gdat.meanreds < 1.56)[0]
        indxtemp = meshgrid(indxcdentemp, gdat.indxredstemp, indexing='ij')
        gdat.diffnabsdiffcdendiffreds[indxtemp] = 10**19.393 * gdat.meancden[indxcdentemp, None]**(-2.) * (1. + gdat.meanreds[None, gdat.indxredstemp])**0.16

        diffoptddiffcdendiffreds = gdat.diffnabsdiffcdendiffreds[:, None, :] * (1. - exp(-gdat.meancden[:, None, None] * gdat.csecionz[None, :, None]))
        diffoptddiffreds = trapz(diffoptddiffcdendiffreds, gdat.meancden, axis=0)
        gdat.diffoptddiffredsrydb = interp1d(gdat.meanenph, diffoptddiffreds, axis=0)(gdat.enerrydb)
        gdat.diffoptddiffcdendiffredsrydb = interp1d(gdat.meanenph, diffoptddiffcdendiffreds, axis=1)(gdat.enerrydb)
        gdat.odep = zeros((gdat.numbenph, gdat.numbreds, gdat.numbreds))
        for c in gdat.indxreds:
            for g in range(c + 1, gdat.numbreds):
                gdat.odep[:, c, g] = trapz(diffoptddiffreds[:, c:g], gdat.meanreds[c:g], axis=1)
        gdat.odeprydb = interp1d(gdat.meanenph, gdat.odep, axis=0)(gdat.enerrydb)
        
        if gdat.makeplot:
            plot_igma(gdat)
                
    # compute the differential ICS power as a function of final photon and electron energies
    gdat.specinvccatl = zeros((gdat.numbenph, gdat.numbenel, gdat.numbreds, 2))
    gdat.numbqics = 100
    gdat.maxmqics = 1.
    for b in gdat.indxenel:
        for a in gdat.indxenph:
            elecrelefact = gdat.meanenel[b] / gdat.masselec
            eps = gdat.meanenph[a] / gdat.meanenel[b]
            if eps >= 1.:
                continue
            for c in gdat.indxreds:
                minmqics = 1. / 4. / elecrelefact**2
                meanqics = logspace(log10(minmqics), log10(gdat.maxmqics), gdat.numbqics)
                enpitemp = gdat.masselec**2 / 4. / gdat.meanenel[b] * eps / (1. - eps) / meanqics
                qfac = 2. * meanqics * log(meanqics) + meanqics + 1. - 2. * meanqics**2 + eps**2 * (1. - meanqics) / 2. / (1. - eps)
                gdat.specinvctemp = 3. * gdat.csecthom * gdat.velolght / 4. / elecrelefact**2 * (gdat.meanenph[a] - enpitemp) * qfac / meanqics / gdat.meanenph[a] # [cm^3/s]
                indxqicstemp = where((enpitemp < gdat.maxmenpi) & (enpitemp > gdat.minmenpi))[0]
                
                # interpolate the energy densities to the current redshift
                edencmbrtemp = interp1d(gdat.meanenpi, gdat.edencmbr[:, c])(enpitemp[indxqicstemp])
                edenegbltemp = interp1d(gdat.meanenpi, gdat.edenegbl[:, c])(enpitemp[indxqicstemp])

                # integrate contibutions from previous redshift
                gdat.specinvccatl[a, b, c, 0] = trapz(gdat.specinvctemp[indxqicstemp] * edencmbrtemp, meanqics[indxqicstemp]) # [1/s/MeV] 
                gdat.specinvccatl[a, b, c, 1] = trapz(gdat.specinvctemp[indxqicstemp] * edenegbltemp, meanqics[indxqicstemp]) # [1/s/MeV]
        
    # temp
    gdat.specinvc = gdat.specinvccatl[:, :, :, 0]
    gdat.specinvccatlenph = trapz(gdat.specinvccatl, gdat.meanenph, axis=0) # [1/s]
    
    # star formation rate density
    path = gdat.pathdata + '/sfrd.csv'
    sfrddata = loadtxt(path)
    gdat.meanredssfrd = sfrddata[:, 0]
    gdat.sfrd = sfrddata[:, 1]
    gdat.emiselecsfrd = 1e51 * gdat.ergs2mgev / (1e3 * gdat.kprc2cmet)**3 / (1e-6 * gdat.myrs2secd) * 0.15 * gdat.sfrd * 1e-2 * (1. + gdat.meanredssfrd)**3
         
    gdat.csecvelopivt = 3e-26
    gdat.csecfracpivt = 1e4
    gdat.masspartpivt = 1e4
    gdat.dmatsloppivt = 1.
    
    gdat.masspart = 5e3 # [MeV]

    path = gdat.pathdata + 'data_%s.fits' % gdat.rtag
    # temp
    if gdat.saveflux:
        if os.path.isfile(path):
            print 'Reading %s...' % path
            gdat.ionrdmat = pf.getdata(path, 0)
            gdat.fluxphotdmat = pf.getdata(path, 1)
        else:
            gdat.csecvelo = gdat.csecvelopivt
            gdat.csecfrac = gdat.csecfracpivt
            gdat.masspart = gdat.masspartpivt
            gdat.dmatslop = gdat.dmatsloppivt
            retr_fluxphotdmat(gdat)
            pf.append(path, gdat.ionrdmat)
            pf.append(path, gdat.fluxphotdmat)

            if gdat.makeplot:
                plot_halo(gdat)
                plot_elec_flux(gdat)
                plot_invrcomp(gdat)
                plot_fluxphot(gdat)
                plot_ionr(gdat)
                plot_sfrd(gdat)

    optiprop = False

    numbproc = 1
    thissamp = rand(numbproc * gdat.numbpara).reshape((numbproc, gdat.numbpara))

    numbplotside = gdat.numbpara
    sampbund = tdpy.mcmc.init(retr_llik, datapara, initsamp=thissamp, \
                              numbswep=numbswep, \
                              factthin=100, \
                              gdatextr=gdat, optiprop=optiprop, \
                                    verbtype=verbtype, pathdata=gdat.pathdata, pathimag=gdat.pathimag, numbplotside=numbplotside)
    
    listsampvarb = sampbund[0]
    listsamp = sampbund[1]
    listsampcalc = sampbund[2]
    listllik = sampbund[3]
    listaccp = sampbund[4]
    listindxsampvari = sampbund[5]
    
    listfluxphotdmat = listsampcalc[0]
    
    numbsamp = listsamp.shape[0]
            
    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    uppr = gdat.meanenph[gdat.indxenphexpr] * (gdat.fluxphotexpr + sqrt(gdat.fluxphotexprvari))
    lowr = gdat.meanenph[gdat.indxenphexpr] * (gdat.fluxphotexpr - sqrt(gdat.fluxphotexprvari))
    axis.fill_between(gdat.meanenph[gdat.indxenphexpr] * 1e6, uppr, lowr, color='lightblue')
    axis.loglog(gdat.meanenph[gdat.indxenphexpr] * 1e6, gdat.meanenph[gdat.indxenphexpr] * gdat.fluxphotexpr, label='ROSAT')
    axis.loglog(gdat.meanenph * 1e6, gdat.meanenph * gdat.fluxphothm12[:, 0], label='Haardt & Madau (2012)')
    tdpy.util.plot_braz(axis, gdat.meanenph * 1e6, gdat.meanenph[None, :] * listfluxphotdmat, alpha=0.5, mcol='black', labl='DM')
    axis.set_xlabel(r'$E_\gamma$ [eV]')
    axis.set_ylabel(r'$EdN/dE$ [1/cm$^2$/s/sr]')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'fluxphotpost.pdf')
    plt.close()

def cnfg_nomi():
    
    init()


if __name__ == '__main__':
    
    cnfg_nomi()

