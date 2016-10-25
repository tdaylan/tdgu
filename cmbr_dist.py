from __init__ import *

def plot_cmbrtemppsec():
    
    numbside = 256
    lghp, bghp, numbpixl, apix = tdpy.util.retr_healgrid(numbside) 

    path = gdat.pathdata + 'COM_PowerSpect_CMB-TT-loL-full_R2.01.txt'
    plnkdata = loadtxt(path)
    lmodloww = plnkdata[:, 0]
    clhploww = plnkdata[:, 1]
    path = gdat.pathdata + 'COM_PowerSpect_CMB-TT-hiL-full_R2.01.txt'
    plnkdata = loadtxt(path)
    lmodhigh = plnkdata[:, 0]
    clhphigh = plnkdata[:, 1]

    flux = hp.synfast(clhphigh, numbside)
    fluxcart = tdpy.util.retr_cart(flux)
    figr, axis = plt.subplots()
    imag = axis.imshow(fluxcart, origin='lower', cmap='Reds')
    plt.colorbar(imag, ax=axis, fraction=0.05)

    flux = hp.synfast(clhploww, numbside)
    fluxcart = tdpy.util.retr_cart(flux)
    figr, axis = plt.subplots()
    imag = axis.imshow(fluxcart, origin='lower', cmap='Reds')
    plt.colorbar(imag, ax=axis, fraction=0.05)

    figr, axis = plt.subplots()
    axis.plot(lmodhigh, clhphigh)
    axis.plot(lmodloww, clhploww)               


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

    figr, axis = plt.subplots()
    axis.loglog(gdat.meanreds, tempcmbr, label='CMB')
    axis.loglog(gdat.meanreds, tempdmatwarm, label='WDM')
    axis.loglog(gdat.meanreds, tempdmatcold, label='CDM')
    axis.loglog(gdat.meanreds, tempmatt, label='Baryonic matter')
    axis.set_yscale('log')
    axis.set_xscale('log')
    axis.set_xlabel(r'$z$')
    axis.set_ylabel(r'T(z) [K]')

    enerradilate = gdat.tempcmbrconc * (1. + gdat.meanreds) * gdat.consboltnatu

    axistwin = axis.twiny()
    axistwin.set_xscale('log')
    axistwin.set_xlim([amin(enerradilate), amax(enerradilate)])
    axistwin.set_xlabel(r'$E_\gamma$ [eV]')

    axis.legend(loc=2)
    plt.savefig(gdat.pathimag + 'tempevol.pdf')
    plt.close()
    

def plot_silkscal():
    
    wnumsilk = (5.92e10)**(-0.5) * (1. + gdat.meanreds)**1.5
    wlensilk = 2. * pi / wnumsilk

    wlenahor = (1. + gdat.meanreds) * gdat.time * yr2s / mp2m * gdat.velolght / sqrt(3. * (1. + edenbmat / gdat.edenradi))
    wnumahor = 2. * pi / wlenahor

    wlenphor = (1. + gdat.meanreds) * gdat.time * yr2s / mp2m * gdat.velolght
    wnumphor = 2. * pi / wlenphor

    figr, axis = plt.subplots()
    axis.loglog(gdat.meanreds, wnumsilk, label='Dissipation scale')
    axis.set_ylabel('$k$ [Mpc$^{-1}$]')
    axis.set_xlabel('$z$')
    axis.axvline(5e4, ls='--', color='black')
    axis.axvline(2e6, ls='--', color='black')
    axis.axhline(interp1d(gdat.meanreds, wnumsilk)(5e4), ls='--', color='black')
    axis.axhline(interp1d(gdat.meanreds, wnumsilk)(2e6), ls='--', color='black')

    plt.figtext(0.19, 0.2, 'y-type distortion', fontsize=20)
    plt.figtext(0.52, 0.2, r'$\mu$-type distortion', fontsize=20)
    plt.figtext(0.78, 0.2, 'Black body', fontsize=20)
    axistwin = axis.twinx()
    axis.set_ylim([amin(1e-3), amax(wnumsilk)])
    axistwin.set_ylim([2. * pi / 1e-3, amin(wlensilk)])
    axistwin.set_yscale('log')
    axistwin.set_ylabel(r'$\lambda$ [Mpc]')
    axis.set_title('Silk dissipation scale', fontsize=20)

    plt.savefig(gdat.pathimag + 'wlensilk.pdf')
    figr.subplots_adjust(top=0.9)
    plt.close()


def plot_resi_post(gdat, postflux):

    figr = plt.figure(figsize=(16, 12))
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 3]) 
    axis = [] 
    axis.append(plt.subplot(gs[0]))
    axis.append(plt.subplot(gs[1], sharex=axis[0]))

    axis[1].set_title('')
    axis[1].errorbar(gdat.freqexpr * 1e-9, gdat.dataflux * 1e-6, ls='none', yerr=gdat.datafluxstdv*1e-6, label=gdat.datalabl, marker='o', markersize=5, color='k')
    
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxtotl * 1e-6, lcol='salmon', alpha=0.5, dcol='red', mcol='black')
    if gdat.inclcmbrmono:
        tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxcmbr * 1e-6, lcol='lightblue', alpha=0.5, dcol='blue', mcol='black')
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxdustcold * 1e-6, lcol='lightgreen', alpha=0.5, dcol='green', mcol='black')
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxdustwarm * 1e-6, lcol='lightgreen', alpha=0.5, dcol='green', mcol='black')
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxsync * 1e-6, lcol='lightyellow', alpha=0.5, dcol='yellow', mcol='black')
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxfree * 1e-6, lcol='lightcoral', alpha=0.5, dcol='coral', mcol='black')
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxydis * 1e-6, lcol='lightyellow', alpha=0.5, dcol='yellow', mcol='black')
    tdpy.mcmc.plot_braz(axis[1], gdat.freqmodl * 1e-9, listfluxdeca * 1e-6, lcol='lightcyan', alpha=0.5, dcol='cyan', mcol='black')
    
    axis[1].set_xscale('log')
    axis[1].set_yscale('log')
    axis[1].set_xlabel(r'$\nu$ [GHz]')
    axis[1].set_ylabel(r'$I_\nu$ [MJy/sr]')
    axis[1].set_ylim([1e-7, 1e4])
    axis[1].legend(loc=9, ncol=4)

    tdpy.mcmc.plot_braz(axis[0], gdat.freqexpr * 1e-9, listresiflux, lcol='lightgrey', alpha=0.5, dcol='darkgrey', mcol='black')
    axis[0].set_ylabel(r'$I_\nu^{res}$ [Jy/sr]')
    axis[0].axhline(0., ls='--', color='black', alpha=0.1)

    plt.savefig(path)
    plt.close(figr)


def plot_flux(gdat, fluxtotl=None, fluxtotlintp=None, fluxcmbr=None, fluxdustcold=None, fluxdustwarm=None, \
                                            fluxsync=None, fluxfree=None, fluxydis=None, fluxdeca=None, plottype='mock'):
    
    if plottype == 'data':
        frac = 1.
    else:
        frac = 1.25
    figr = plt.figure(figsize=(2 * gdat.plotsize, frac * gdat.plotsize))
    
    if plottype == 'data':
        axis = [plt.gca()]
    else:
        axisgridtemp = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 3])
        axisgridtemp.update(hspace=0.1)
        axis = [] 
        axis.append(plt.subplot(axisgridtemp[1]))
        axis.append(plt.subplot(axisgridtemp[0], sharex=axis[0]))

    # plot data
    yerr = gdat.datafluxstdv * 1e-6

    axis[0].errorbar(gdat.freqexpr * 1e-9, gdat.dataflux * 1e-6, ls='', yerr=yerr, label=gdat.datalabl, marker='o', color='black')
    
    if plottype != 'data':
        axis[0].plot(gdat.freqmodlscal, gdat.mockfluxtotl * 1e-6, label='Total model', color='b')
        if gdat.inclcmbrmono:
            axis[0].plot(gdat.freqmodlscal, gdat.mockfluxcmbr * 1e-6, label='CMB', color='g')
        axis[0].plot(gdat.freqmodlscal, gdat.mockfluxdustcold * 1e-6, label='Cold Dust', color='r', ls='--')
        axis[0].plot(gdat.freqmodlscal, gdat.mockfluxdustwarm * 1e-6, label='Warm Dust', color='r', ls='-.')
        axis[0].plot(gdat.freqmodlscal, gdat.mockfluxsync * 1e-6, label='Synchrotron', color='orange')
        axis[0].plot(gdat.freqmodlscal, gdat.mockfluxfree * 1e-6, label='Brem', color='purple')
        axis[0].plot(gdat.freqmodlscal, fabs(gdat.mockfluxydis) * 1e-6, label='Reionization', color='yellow')
        axis[0].plot(gdat.freqmodlscal, fabs(gdat.mockfluxdeca) * 1e-6, label='Particle decay', color='cyan')
        if plottype == 'samp':
            axis[0].plot(gdat.freqmodlscal, fluxtotl * 1e-6, color='b', alpha=gdat.alphsamp)
            if gdat.inclcmbrmono:
                axis[0].plot(gdat.freqmodlscal, fluxcmbr * 1e-6, color='g', alpha=gdat.alphsamp)
            axis[0].plot(gdat.freqmodlscal, fluxdustcold * 1e-6, color='r', ls='--', alpha=gdat.alphsamp)
            axis[0].plot(gdat.freqmodlscal, fluxdustwarm * 1e-6, color='r', ls='-.', alpha=gdat.alphsamp)
            axis[0].plot(gdat.freqmodlscal, fluxsync * 1e-6, color='orange', alpha=gdat.alphsamp)
            axis[0].plot(gdat.freqmodlscal, fluxfree * 1e-6, color='purple', alpha=gdat.alphsamp)
            axis[0].plot(gdat.freqmodlscal, fabs(fluxydis) * 1e-6, color='yellow', alpha=gdat.alphsamp)
            axis[0].plot(gdat.freqmodlscal, fabs(fluxdeca) * 1e-6, color='cyan', alpha=gdat.alphsamp)
            
    axis[0].set_title('')
    axis[0].set_xscale('log')
    axis[0].set_yscale('log')
    axis[0].set_xlabel(r'$\nu$ [GHz]')
    axis[0].set_ylabel(r'$I_\nu$ [MJy/sr]')
    axis[0].set_ylim([1e-7, 1e4])
    axis[0].legend(loc=2)
    xlim = [gdat.minmfreqmodl * 1e-9, gdat.maxmfreqmodl * 1e-9]
    axis[0].set_xlim(xlim)

    if plottype != 'data':
        if plottype == 'mock':
            resiflux = gdat.dataflux - gdat.mockfluxtotlintp
        else:
            resiflux = gdat.dataflux - fluxtotlintp

        axis[1].errorbar(gdat.freqexpr * 1e-9, resiflux, yerr=gdat.datafluxstdv, marker='o', lw=1, ls='none', markersize=5, color='k')
        axis[1].set_ylabel(r'$I_\nu^{res}$ [Jy/sr]')
        axis[1].axhline(0., ls='--', color='black', alpha=0.1)
        axis[1].set_xlim(xlim)
        maxmresi = max(amax(resiflux + gdat.datafluxstdv), amax(fabs(resiflux - gdat.datafluxstdv)))
        axis[1].set_yticks(linspace(-maxmresi, maxmresi, 5))
        plt.setp(axis[1].get_xticklabels(), visible=False)
    
    if plottype == 'samp':
        strg = '_%d' % gdat.cntrswep
    else:
        strg = ''
    plt.savefig(gdat.pathimag + 'postflux_%s%s.pdf' % (plottype, strg))
    plt.close(figr)


def plot_plnk(gdat):
    
    figr, axis = plt.subplots(figsize=(10, 8))
    for k in arange(gdat.numbfreqplnk):
        axis.loglog(gdat.freqplnk[k] * 1e-9, gdat.tran[k])
    axis.set_ylim([1e-13, 1e-9])
    axis.set_xlim([1e1, 1e3])
    axis.set_ylabel(r'$T(\nu)$')
    axis.set_xlabel(r'$\nu$ [GHz]')
    path = gdat.pathdata + 'plnktran.pdf'
    plt.savefig(path)
    plt.close()
    
    indxpixltemp = random_integers(0, 12*256**2, size=100)
    path = gdat.pathdata + 'plnkflux.dat'
    flux = loadtxt(path)[indxpixltemp, :]
    path = gdat.pathdata + 'plnkfluxstdv.dat'
    fluxstdv = loadtxt(path)[indxpixltemp, :]

    bolo = empty(gdat.numbfreqexpr)
    for i in range(gdat.numbfreqexpr):
        bolo[i] = trapz(tran[i][1:] * flux[0, i] * gdat.freqexpr[i] / freqtran[i][1:], freqtran[i][1:])
        plt.loglog(freqtran[i] * 1e-9, tran[i])  
    plt.loglog(gdat.freqexpr * 1e-9, fluxstdv)

    for k in range(10):
        plt.loglog(gdat.freqexpr * 1e-9, flux[k, :])
    plt.loglog(gdat.freqexpr * 1e-9, bolo)
    plt.loglog(gdat.freqexpr * 1e-9, 100. * (bolo - flux) / flux)
    
    strgpara = []
    for k in range(gdat.numbfreqplnk):
        strgpara.append(r'$\mathcal{I}_{%d}$' % (gdat.freqplnk[k] / 1e9))
    scalpara = ['self'] * gdat.numbfreqplnk
    path = gdat.pathimag + 'plnkcros'
    tdpy.mcmc.plot_grid(path, gdat.plnkflux * 1e-6, strgpara, scalpara=scalpara)
    
    for k in range(gdat.numbfreqplnk):
        hmapcart = tdpy.util.retr_cart(gdat.plnkflux[:, k])


def plot_intr(eventype='norm'):
    
    with plt.xkcd():

        def sigd(varb, varbcntr):
            sigd = 1. / (1 + exp(-0.1 * (varb - varbcntr)))
            return sigd / amax(sigd)

        figr, axis = plt.subplots(figsize=(14, 6))

        redsxkcd = arange(1000)
        rele = sigd(redsxkcd, 300.) - sigd(redsxkcd, 500.) + sigd(redsxkcd, 950.)

        axis.plot(redsxkcd, rele)

        axis.text(200, -0.1, "z=One Million")
        axis.text(450, -0.1, "z=One Thousand")
        axis.text(650, -0.1, "z=One Hundred")
        axis.text(900, -0.1, "04/01/2016")
        axis.set_xticks([])
        if eventype == 'scra':
            axis.annotate('Nucleosynthesis', xy=[950, 0.3], xytext=[20, 0.7], rotation=45)
            axis.annotate('Reheating', xy=[950, 0.3], xytext=[120, 0.7], rotation=45)
            axis.annotate('Inflation', xy=[950, 0.3], xytext=[200, 0.7], rotation=45)
            axis.annotate('CMB Spectral Distortions', xy=[950, 0.3], xytext=[300, 0.45], rotation=45)
            axis.annotate('Reionization', xy=[550, 0.3], xytext=[550, 0.7], rotation=45)
            axis.annotate('Structure formation', xy=[750, 0.3], xytext=[700, 0.7], rotation=45)
            axis.annotate('Recombination', xy=[750, 0.3], xytext=[650, 0.7], rotation=45)
            axis.set_title("History of the Universe")
        else:
            axis.annotate('Inflaton', xy=[950, 0.3], xytext=[20, 0.7], rotation=45)
            axis.annotate('Reheating', xy=[950, 0.3], xytext=[70, 0.7], rotation=45)
            axis.annotate('Nucleosynthesis', xy=[950, 0.3], xytext=[150, 0.7], rotation=45)
            axis.annotate('CMB spectral distortions', xy=[950, 0.3], xytext=[300, 0.45], rotation=45)
            axis.annotate('Dark ages', xy=[550, 0.3], xytext=[550, 0.7], rotation=45)
            axis.annotate('Structure formation', xy=[750, 0.3], xytext=[700, 0.7], rotation=45)
            axis.annotate('Reionization', xy=[850, 0.3], xytext=[850, 0.7], rotation=45)
            axis.set_title("History of the Universe")

        axis.set_xlabel('Cosmic Time')
        axis.set_ylabel("Relevance to Today's Talk")

        path = gdat.pathdata + 'talkintr_%s.pdf' % eventype
        plt.savefig(path)
        plt.close()


def plot_fluxdust(gdat, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx):
    
    logtdustodep = -4.
    logtdustemisrati  = 1.
    dustpowrfrac = 0.05
    dustwarmindx = 2.5
    dustwarmtemp = 20.
    
    fluxdust, fluxdustcold, fluxdustwarm = retr_fluxdust(gdat, gdat.freqmodl, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx)
    
    figr, axis = plt.subplots(figsize=(12, 6))
    axis.loglog(gdat.freqmodl * 1e-9, 1e-6 * fluxdust, label='Total')
    axis.loglog(gdat.freqmodl * 1e-9, 1e-6 * fluxdustcold, label='Cold')
    axis.loglog(gdat.freqmodl * 1e-9, 1e-6 * fluxdustwarm, label='Warm')
    axis.set_title('Two-component thermal dust SED')
    axis.set_xlabel(r'$\nu$ [GHz]')
    axis.set_ylabel(r'$I_\nu$ [MJy/sr]')
    axis.legend(loc=2)
    axis.set_ylim([1e-6, 1e5])
    plt.savefig(gdat.pathimag + 'visifunc.pdf')
    plt.close()


def plot_grnf(gdat):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.meanreds, gdat.vify, label=r'$\mathcal{J}_y(z_i)$')
    axis.plot(gdat.meanreds, gdat.vift, label=r'$\mathcal{J_B}(z_i)$')
    axis.plot(gdat.meanreds, gdat.vifm, label=r'$\mathcal{J}_\mu(z_i)$')
    axis.plot(gdat.meanreds, gdat.vifm * gdat.vift, label=r'$\mathcal{J}_B(z_i)\mathcal{J}_\mu(z_i)$')
    axis.set_xscale('log')
    axis.set_xlabel('$z$')
    axis.set_ylabel(r'$\mathcal{J}(z)$')
    axis.legend()
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'visifunc.pdf')
    plt.close()
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.meanreds, gdat.vifm * gdat.vift + gdat.vify - gdat.vift)
    axis.set_xscale('log')
    axis.set_xlabel('$z$')
    axis.set_ylabel(r'$\Delta\mathcal{J}(z)$')
    figr.tight_layout()
    plt.savefig(gdat.pathimag + 'errrgren.pdf')
    plt.close()

    with sns.color_palette("Blues", gdat.numbredsplot):
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        for c in range(gdat.numbredsplot):
            axis.plot(gdat.freqmodl * 1e-9, 1e-6 * gdat.fluxfdisgrenconc[:,gdat.indxredsplot[c]], label='$z$ = %3.1g' % gdat.meanreds[gdat.indxredsplot[c]])
        axis.set_xscale('log')
        axis.set_xlabel(r'$\nu$ [GHz]')
        axis.set_ylabel('$\Delta I_\nu^G$')
        plt.savefig(gdat.pathimag + 'distfull.pdf')
        plt.close()
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        for c in range(gdat.numbredsplot):
            axis.plot(gdat.freqmodl * 1e-9, 1e-6 * gdat.fluxodisgrenconc[:,gdat.indxredsplot[c]], label='$z$ = %3.1g' % gdat.meanreds[gdat.indxredsplot[c]])
        axis.set_xscale('log')
        axis.set_xlabel(r'$\nu$ [GHz]')
        axis.set_ylabel('$\Delta I_\nu^{O}$')
        plt.savefig(gdat.pathimag + 'distobsv.pdf')
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

    figr, axis = plt.subplots()
    figr.suptitle('Time scales relevant to Comptonization in the early Universe', fontsize=20)

    #axis.loglog(gdat.meanreds, timebrem, label='Bremsstrahlung')
    #axis.loglog(gdat.meanreds, timedcom, label='Double Compton')
    #axis.loglog(gdat.meanreds, timecomp, label='Compton')
    axis.loglog(gdat.meanreds, timecoul, label='Coulomb')
    axis.loglog(gdat.meanreds, timebose, label='Compton')
    axis.loglog(gdat.meanreds, timether, label='Thermalization')
    axis.loglog(gdat.meanreds, gdat.timehubb / yr2s, color='black', label='Hubble')

    line0 = axis.axvline(5e4, ls='--', color='grey')
    line1 = axis.axvline(2e6, ls='-.', color='grey')

    leg0 = plt.legend([line0, line1], ['Compton Surface', 'Blackbody Surface'], loc='center', bbox_to_anchor=[0.2, 0.1])
    leg1 = axis.legend(loc='center', bbox_to_anchor=[0.65, 0.83])
    axis.add_artist(leg0)

    axis.set_xlabel('$z$')
    axis.set_ylabel(r'$\tau$ [yr]')

    axistwin = axis.twiny()
    axistwin.set_xscale('log')
    axistwin.set_xlim([amax(gdat.time[:-1]), amin(gdat.time[:-1])])
    axistwin.set_xlabel('$t$ [yr]')

    figr.subplots_adjust(top=0.85)
    plt.savefig(gdat.pathimag + 'timescal.pdf')
    plt.close()


def plot_dist(dist):
    figr, axis = plt.subplots()
    axis.set_xscale('log')
    axis.set_xlabel(r'$\nu$ [GHz]')
    axis.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')


def plot_pure_dist(gdat):
    figr, axis = plt.subplots()
    axis.plot(gdat.freqmodl * 1e-9, gdat.fluxcmbrconc * 1e-6, color='black', label=r'$I^B_\nu(\nu)$')
    axis.plot(gdat.freqmodl * 1e-9, gdat.fluxtdisgrenconc * 1e-6, label=r'$I^T_\nu(\nu)$')
    axis.plot(gdat.freqmodl * 1e-9, gdat.fluxydisgrenconc * 1e-6, label=r'$I^Y_\nu(\nu)$')
    axis.plot(gdat.freqmodl * 1e-9, gdat.fluxmdisgrenconc * 1e-6, label=r'$I^M_\nu(\nu)$')
    axis.set_xscale('log')
    axis.set_xlabel(r'$\nu$ [GHz]')
    axis.set_ylabel(r'$I_\nu$ [MJy/sr]')
    axistwin = axis.twinx()
    axistwin.set_ylim(array(axis.get_ylim()) * 1e-20)
    axistwin.set_ylabel('[J/m$^2$/s/sr/Hz]')
    axistwin = axis.twiny()
    axistwin.set_xlim([amin(gdat.sfrqconc), amax(gdat.sfrqconc)])
    axistwin.set_xlabel('$x$')
    axistwin.set_xscale('log')
    plt.legend(loc=2)
    plt.savefig(gdat.pathimag + 'puredist.pdf')
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

    figr, axis = plt.subplots()
    for k in range(1, ninjetype):
        axis.loglog(gdat.meanreds, gdat.meanreds * heat[:, k], label=heatstrg[k])
    axis.set_xlabel(r'$z$')
    axis.set_ylabel(r'$d\kappa/d\ln z$')
    axis.set_ylim([1e-12, 1e-6])
    axis.legend(loc=2)

    axistwin = axis.twiny()
    axistwin.set_xscale('log')
    axistwin.set_xlim([amax(gdat.time[:-1]), amin(gdat.time[:-1])])
    axistwin.set_xlabel('$t$ [yr]')

    plt.savefig(gdat.pathimag + 'heatrate.pdf')
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

    figr, axis = plt.subplots()
    for k in range(1, 4):
        axis.loglog(gdat.meanredsextd, gdat.meanredsextd * heatextd[:, k], label=heatstrg[k])
    axis.set_xlabel(r'$z$')
    axis.set_ylabel(r'$d\kappa/d\ln z$')
    #axis.set_ylim([1e-12, 1e-1])
    axis.legend(loc=2)

    axistwin = axis.twiny()
    axistwin.set_xscale('log')
    axistwin.set_xlim([amax(timeextd[:-1]), amin(timeextd[:-1])])
    axistwin.set_xlabel('$t$ [s]')
    axistwin.axvline(1., ls='--', color='grey')
    axistwin.annotate('BBN', xy=[1., 1e-12], xytext=[1e2, 1e-12], arrowprops=dict(arrowstyle="->"), fontsize=20)
    plt.savefig(gdat.pathimag + 'heatextd.pdf')
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
        disttemp, heattemp, redsdecatemp = retr_fluxdeca(gdat, gdat.fluxodisgrenconc, ampldeca, timedecalist[k])
        heatdeca[:, k] = heattemp
        heatdeca[:, k] = dist
        redsdeca[k] = redsdecatemp


    figr, axis = plt.subplots()
    for k in range(njdeca):
        axis.loglog(gdat.meanreds, gdat.meanreds * heatdeca[:, jdeca[k]], label=r'$\tau = %d$ yr' % timedecalist[jdeca[k]])
        axis.set_xlabel(r'$z$')
        axis.set_ylabel(r'$d\kappa/d\ln z$')
        axis.set_ylim([1e-10, 1e-6])
        axis.legend(loc=2, ncol=4)

        axistwin = axis.twiny()
        axistwin.set_xscale('log')
        axistwin.set_xlim([amax(gdat.time[:-1]), amin(gdat.time[:-1])])
        axistwin.set_xlabel('$t$ [yr]')

    plt.savefig(gdat.pathimag + 'heatdeca.pdf')
    plt.close()
    

    timedecaplot = logspace(-4., 6, 10)
    with sns.color_palette("Blues", gdat.numbredsplot): 
        figr, axis = plt.subplots()    
        for k in range(timedecaplot.size):
            axis.loglog(gdat.freqmodl * 1e-9, abs(retr_fluxdeca(gdat, fluxodisgren, 1., timedecaplot[k])), label='$t = %.3g$' % timedecaplot[k])
        axis.set_xscale('log')
        axis.set_xlabel(r'$\nu$ [GHz]')
        axis.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')
        axis.legend(loc=4)
        plt.close()


def plot_sampdist():
    diffdistdiffreds = heat[None, :, :] * fluxodisgren[:, :, None]
    dist = zeros((nfreq, ninjetype))
    for k in range(ninjetype):
        dist[:, k] = trapz(diffdistdiffreds[:, :, k], gdat.meanreds, axis=1)

    figr, axis = plt.subplots()
    for k in range(1, ninjetype):
        axis.plot(gdat.freqmodl * 1e-9, dist[:, k], label=heatstrg[k])
    axis.set_xscale('log')
    axis.set_xlabel(r'$\nu$ [GHz]')
    axis.set_ylabel(r'$\Delta I_\nu$ [Jy/sr]')
    axis.legend(loc=2)
    axis.axhline(5., ls='--', color='grey')
    axis.axhline(-5., ls='--', color='grey')
    axis.fill_between(gdat.freqmodl * 1e-9, ones_like(gdat.freqmodl) * 5., ones_like(gdat.freqmodl) * -5., color='grey')
    axis.text(2, 10, r'PIXIE 1-$\sigma$ sensitivity', color='grey', fontsize=20)
    plt.savefig(gdat.pathimag + 'totldist.pdf')
    plt.close()
    
    
    diffdistdifflred = heat[None, :, :] * fluxodisgren[:, :, None] * gdat.meanreds[None, :, None]
    for k in range(1, 2):
        ylim = [amin(diffdistdifflred[:, :, k]), amax(diffdistdifflred[:, :, k])]
        for c in range(gdat.numbredsplot):
            figr, axis = plt.subplots()
            axis.plot(gdat.freqmodl * 1e-9, diffdistdifflred[:, numbreds-1-gdat.indxredsplot[c], k])
            text = plt.figtext(0.2, 0.8, '$z_i$ = %.3g' % gdat.meanreds[nreds-gdat.indxredsplot[c]-1], fontsize=20)

            axis.set_xscale('log')
            axis.set_xlabel(r'$\nu$ [GHz]')
            axis.set_ylabel(r'$d\Delta I_\nu/ d\ln z$ [Jy/sr]')
            axis.set_ylim(ylim)
            axis.set_title(heatstrg[k])

            plt.savefig(gdat.pathimag + 'diffdistdifflred_' + heatrtag[k] + '_%d.pdf' % gdat.indxredsplot[c])
            plt.close(figr)
        

def retr_fluxdeca(gdat, fluxodisgren, ampldeca, timedeca):

    ratedeca = 1. / timedeca
    redsdeca = interp1d(gdat.time, gdat.meanreds)(timedeca)
    
    heatdeca = ampldeca * gdat.ndenbmat * ratedeca * exp(-timedeca / gdat.time) * gdat.timehubb / (1. + gdat.meanreds) / gdat.edenradi

    difffluxdecadiffreds = heatdeca[None, :] * fluxodisgren
    fluxdeca = trapz(difffluxdecadiffreds, gdat.meanreds, axis=1)

    return fluxdeca


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
    
    sfrq = gdat.consplnk * thisfreq / gdat.consbolt / tempcmbr
    
    occpplnk = retr_occp(sfrq)

    fluxcmbr = 2. * gdat.consplnk * thisfreq**3 / gdat.velolght**2 * occpplnk * 1e26
    
    return fluxcmbr


def retr_fluxgren(gdat, thisfreq, tempcmbr, disttype):
    
    sfrq = gdat.consplnk * thisfreq / gdat.consbolt / tempcmbr

    occpplnk, occptdis, occpydis, occpmdis = retr_occp(sfrq, dist=True)
    
    fluxtdisgren = 0.25 * 2. * gdat.consplnk * thisfreq**3 / gdat.velolght**2 * occptdis * 1e26
    fluxydisgren = 0.25 * 2. * gdat.consplnk * thisfreq**3 / gdat.velolght**2 * occpydis * 1e26
    fluxmdisgren = 1.41 * 2. * gdat.consplnk * thisfreq**3 / gdat.velolght**2 * occpmdis * 1e26
    
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
    gaunfact = log(4.955e2 * (thisfreq / 1e9)**(-1.)) + 1.5 * log(tempfree)
    odepfree = 3.014e-2 * (tempfree)**(-1.5) * (thisfreq / 1e9)**(-2.) * emmefree * gaunfact
    fluxfree = (1. - exp(-odepfree)) * thisplnkfunc
    
    return fluxfree


def retr_plnkfunc(gdat, thisfreq, temp):
    
    thissfrq = gdat.consplnk * thisfreq / gdat.consbolt / temp
    
    plnk = 2. * gdat.consplnk * thisfreq**3 / gdat.velolght**2 / (exp(thissfrq) - 1.)
    
    return plnk


def retr_ydis_trac():
    
    path = gdat.pathimag + 'electron_photonspectrum_results.fits'
    data, hdr = pf.getdata(path, 1, header=True)

    redsoutp = data['OUTPUT_REDSHIFT'].squeeze()
    enerpart = 10**data['LOG10ENERGY'].squeeze()
    redsinpt = data['INPUT_REDSHIFT'].squeeze()
    enerphot = data['PHOTON_ENERGY'].squeeze()
    freqphot = enerphot / gdat.consplnk

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

    specchek = 8. * pi * enerphot[:,:,None]**2 / (gdat.velolght * gdat.consplnk)**3 / (exp(enerphot[:,:,None] / \
                                                    tempcmbrnunc / gdat.consbolt / (1. + redsinpt[None,None,:])) - 1.) # [1/cm^3/eV]

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
                labl = '$z_{out} = %.3g, E_{inj} = %.3g, z_{inp} = %.3g$' % (redsoutp[jredsoutp[c]], enerpart[jenerpart[a]], redsinpt[jredsinpt[f]])
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

    figr, axis = plt.subplots(figsize=(10,10));
    line, = axis.plot(freqphot[:,0], spec[:,0,c])
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.close()


def retr_flux(gdat, thisfreq, thispara):

    # CMB
    if modlcmbrtype == 'full':
        gdatmodi.flux.cmbr = retr_fluxcmbr(gdat, thisfreq, thispara.tempcmbr)
    

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
        
    # dust
    if gdat.modldusttype == 'singtemp':
        fluxdust, fluxdustcold, fluxdustwarm = retr_fluxdust_singtemp(gdat, thisfreq, dustodep, dustwarm, dusttemp)
    if gdat.modldusttype == 'doubtemp':
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

    fluxtotl, fluxcmbr, fluxdustcold, fluxdustwarm, fluxsync, fluxfree, fluxydis, fluxdeca
    return flux


def retr_fluxdoubdust(gdat, thisfreq, dustodep, dustemisrati, dustpowrfrac, dustwarmindx, dustwarmtemp, dustcoldindx):

    factwarm = (sp.special.zetac(4. + dustwarmindx) + 1) * sp.special.gamma(4. + dustwarmindx)
    factcold = (sp.special.zetac(4. + dustcoldindx) + 1) * sp.special.gamma(4. + dustcoldindx)
    dustcoldtemp = (factwarm / factcold / dustemisrati * (gdat.consplnk * gdat.freqpivt / gdat.consbolt)**(dustcoldindx - dustwarmindx) * \
            dustwarmtemp**(4. + dustwarmindx))**(1. / (4. + dustcoldindx))
    
    fluxdustcoldfact = dustpowrfrac * dustemisrati * (thisfreq / 3e12)**dustcoldindx
    fluxdustwarmfact = (1. - dustpowrfrac) * (thisfreq / 3e12)**dustwarmindx
        
    fluxdustcold = dustodep * fluxdustcoldfact * retr_plnkfunc(gdat, thisfreq, dustcoldtemp) * 1e26 / (fluxdustcoldfact + fluxdustwarmfact)
        
    fluxdustwarm = dustodep * fluxdustwarmfact * retr_plnkfunc(gdat, thisfreq, dustwarmtemp) * 1e26 / (fluxdustcoldfact + fluxdustwarmfact)
        
    fluxdust = fluxdustcold + fluxdustwarm
    
    return fluxdust, fluxdustcold, fluxdustwarm


def retr_llik(sampvarb, gdat, gdatintr):
    
    modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, modlfluxfree, modlfluxydis, modlfluxdeca = retr_flux(gdat, gdat.freqmodl, sampvarb)
    

    modlfluxtotlintp = interp1d(gdat.freqmodl, modlfluxtotl)(gdat.freqexpr)
    resiflux = gdat.dataflux - modlfluxtotlintp

    #llik = sum(-log(sqrt(2. * pi) * gdat.datafluxstdv) - 0.5 * (modlfluxtotlintp - gdat.dataflux) / gdat.datafluxstdv**2 * (modlfluxtotlintp - gdat.dataflux))
    llik = -0.5 * sum((modlfluxtotlintp - gdat.dataflux) / gdat.datafluxstdv**2 * (modlfluxtotlintp - gdat.dataflux))
   
    if False:
        print 'modlfluxtotl'
        print modlfluxtotl
        print 'modlfluxcmbr'
        print modlfluxcmbr
        print 'modlfluxdustcold'
        print modlfluxdustcold
        print 'modlfluxdustwarm'
        print modlfluxdustwarm
        print 'modlfluxsync'
        print modlfluxsync
        print 'modlfluxfree'
        print modlfluxfree
        print 'modlfluxydis'
        print modlfluxydis
        print 'modlfluxdeca'
        print modlfluxdeca

        if gdatintr.cntrswep == 0:
            print 'gdat.dataflux'
            print gdat.dataflux
            print 'gdat.datafluxstdv'
            print gdat.datafluxstdv
            print 'modlfluxtotlintp'
            print modlfluxtotlintp
        print

    if gdat.plotsamp and gdatintr.cntrswep % gdat.numbswepplot == 0:
        gdat.cntrswep = gdatintr.cntrswep
        plot_flux(gdat, fluxtotl=modlfluxtotl, fluxtotlintp=modlfluxtotlintp, fluxcmbr=modlfluxcmbr, fluxdustcold=modlfluxdustcold, fluxdustwarm=modlfluxdustwarm, \
                        fluxsync=modlfluxsync, fluxfree=modlfluxfree, fluxydis=modlfluxydis, fluxdeca=modlfluxdeca, plottype='samp')
    
    if gdat.savepost:
        sampcalc = modlfluxtotl, modlfluxcmbr, modlfluxdustcold, modlfluxdustwarm, modlfluxsync, modlfluxfree, modlfluxydis, modlfluxdeca, resiflux
    else:
        sampcalc = []

    return llik, sampcalc


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
    datapara.true = zeros(numbpara)
    

    datapara.indx['tempcmbr'] = 0
    datapara.name[0] = 'tempcmbr'
    datapara.minm[0] = 2.72
    datapara.maxm[0] = 2.73
    datapara.scal[0] = 'self'
    datapara.labl[0] = '$T_{cmb}$'
    datapara.unit[0] = '[K]'
    datapara.vari[0] = 4.92e-7
    datapara.true[0] = 2.725

    datapara.indx['dustodep'] = 1
    datapara.name[1] = 'dustodep'
    datapara.minm[1] = 1e-4
    datapara.maxm[1] = 1e-2
    datapara.scal[1] = 'logt'
    datapara.labl[1] = r'$\tau_{dust}$'
    datapara.unit[1] = ''
    datapara.vari[1] = 1.33e-8
    datapara.true[1] = 1e-3

    datapara.indx['dustemisrati'] = 2
    datapara.name[2] = 'dustemisrati'
    datapara.minm[2] = 1e0
    datapara.maxm[2] = 1e2
    datapara.scal[2] = 'logt'
    datapara.labl[2] = '$q_1/q_2$'
    datapara.unit[2] = ''
    datapara.vari[2] = 4.34e-8
    datapara.true[2] = 1e1
    
    datapara.indx['dustpowrfrac'] = 3
    datapara.name[3] = 'dustpowrfrac'
    datapara.minm[3] = 0.
    datapara.maxm[3] = 0.05
    datapara.scal[3] = 'self'
    datapara.labl[3] = '$f_1$'
    datapara.unit[3] = ''
    datapara.vari[3] = 1.72e-7
    datapara.true[3] = 0.025
    
    datapara.indx['dustwarmindx'] = 4
    datapara.name[4] = 'dustwarmindx'
    datapara.minm[4] = 2.5
    datapara.maxm[4] = 3.
    datapara.scal[4] = 'self'
    datapara.labl[4] = r'$\beta_2$'
    datapara.unit[4] = ''
    datapara.vari[4] = 7.50e-7
    datapara.true[4] = 2.75
    
    datapara.indx['dustwarmtemp'] = 5
    datapara.name[5] = 'dustwarmtemp'
    datapara.minm[5] = 10.
    datapara.maxm[5] = 20.
    datapara.scal[5] = 'self'
    datapara.labl[5] = '$T_2$'
    datapara.unit[5] = '[K]'
    datapara.vari[5] = 2.58e-8
    datapara.true[5] = 15.
    
    datapara.indx['dustcoldindx'] = 6
    datapara.name[6] = 'dustcoldindx'
    datapara.minm[6] = 1.
    datapara.maxm[6] = 2.5
    datapara.scal[6] = 'self'
    datapara.labl[6] = r'$\beta_1$'
    datapara.unit[6] = ''
    datapara.vari[6] = 2.49e-7
    datapara.true[6] = 1.75
    
    datapara.indx['syncnorm'] = 7
    datapara.name[7] = 'syncnorm'
    datapara.minm[7] = 1e3
    datapara.maxm[7] = 1e7
    datapara.scal[7] = 'logt'
    datapara.labl[7] = '$A_{sync}$'
    datapara.unit[7] = ''
    datapara.vari[7] = 3.96e-5
    datapara.true[7] = 1e5
    
    datapara.indx['syncindx'] = 8
    datapara.name[8] = 'syncindx'
    datapara.minm[8] = -1.5
    datapara.maxm[8] = 0.
    datapara.scal[8] = 'self'
    datapara.labl[8] = r'$\alpha_{sync}$'
    datapara.unit[8] = ''
    datapara.vari[8] = 5.78e-5
    datapara.true[8] = -0.75
    
    datapara.indx['emmefree'] = 9
    datapara.name[9] = 'emmefree'
    datapara.minm[9] = 1e0
    datapara.maxm[9] = 1e4
    datapara.scal[9] = 'logt'
    datapara.labl[9] = 'EM'
    datapara.unit[9] = '[pc/cm$^6$]'
    datapara.vari[9] = 2.44e-6
    datapara.true[9] = 1e2
    
    datapara.indx['tempfree'] = 10
    datapara.name[10] = 'tempfree'
    datapara.minm[10] = 1e1
    datapara.maxm[10] = 1e3
    datapara.scal[10] = 'logt'
    datapara.labl[10] = r'$T_e$'
    datapara.unit[10] = '[K]'
    datapara.vari[10] = 2.55e-5
    datapara.true[10] = 1e2
    
    datapara.indx['ydisampl'] = 11
    datapara.name[11] = 'ydisampl'
    datapara.minm[11] = 1e-9
    datapara.maxm[11] = 1e-5
    datapara.scal[11] = 'logt'
    datapara.labl[11] = '$y_{ri}$'
    datapara.unit[11] = ''
    datapara.vari[11] = 6.11e-3
    datapara.true[11] = 1e-7
    
    datapara.indx['ampldeca'] = 12
    datapara.name[12] = 'ampldeca'
    datapara.minm[12] = 1e-7
    datapara.maxm[12] = 1e-3
    datapara.scal[12] = 'logt'
    datapara.labl[12] = '$f_X$'
    datapara.unit[12] = '[eV]'
    datapara.vari[12] = 7.40e-1
    datapara.true[12] = 1e-5
    
    datapara.indx['timedeca'] = 13
    datapara.name[13] = 'timedeca'
    datapara.minm[13] = 1e8
    datapara.maxm[13] = 1e12
    datapara.scal[13] = 'logt'
    datapara.labl[13] = r'$\tau_X$'
    datapara.unit[13] = '[s]'
    datapara.vari[13] = 7.64e-1
    datapara.true[13] = 1e10
    
    datapara.strg = datapara.labl + ' ' + datapara.unit
    
    return datapara


def retr_egbl(thisfreq):
    
    path = gdat.pathdata + 'egbl.csv'
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
         modldusttype='doubtemp', \
         mockdusttype='doubtemp', \
         datalabl='PIXIE', \
         numbswep=10000, \
         numbburn=100, \
         factthin=99, \
         numbproc=1, \
         #numbproc=20, \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmgang=1., \
         exprflux=None, \
         exprfluxstdv=None, \
         freqexpr=None, \
         numbfreqexpr=None, \
         minmfreqexpr=None, \
         maxmfreqexpr=None, \
         fluxstdvinst=None, \
         fluxstdvfrac=None, \
         inclcmbrmono=True, \
         verbtype=1, \
         makeplot=True, \
         plotsamp=False, \
         numbswepplot=1000, \
        ):
    
    # global object
    gdat = tdpy.util.gdatstrt()
   
    gdat.datatype = datatype
    gdat.exprtype = exprtype
    gdat.datalabl = datalabl
    gdat.verbtype = verbtype
    gdat.numbfreqexpr = numbfreqexpr
    gdat.minmfreqexpr = minmfreqexpr
    gdat.maxmfreqexpr = maxmfreqexpr
    gdat.fluxstdvinst = fluxstdvinst
    gdat.fluxstdvfrac = fluxstdvfrac
    gdat.inclcmbrmono = inclcmbrmono
    gdat.freqexpr = freqexpr
    gdat.plotsamp = plotsamp
    gdat.numbswepplot = numbswepplot
    
    # data path
    gdat.pathdata = tdpy.util.retr_path('tdgu', 'cmbr_dist/', onlydata=True)
   
    # Planck variables
    gdat.freqplnk = array([3e10, 4.4e10, 7e10, 1e11, 1.43e11, 2.17e11, 3.53e11, 5.45e11, 8.57e11]) # [Hz]
    gdat.numbfreqplnk = gdat.freqplnk.size 
    gdat.fluxstdvinst = 5e3 # [Jy/sr]
    gdat.fluxstdvfrac = 1e-3
    gdat.numbpixlplnk = 2048 

    # PIXIE variables
    gdat.numbfreqpixi = 400
    gdat.minmfreqpixi = 3e10 # [Hz]
    gdat.maxmfreqpixi = 6e12 # [Hz]
    gdat.freqpixi = logspace(log10(gdat.minmfreqpixi), log10(gdat.maxmfreqpixi), gdat.numbfreqpixi)
    gdat.fluxstdvinstpixi = 5e0 # [Jy/sr]
    gdat.fluxstdvfracpixi = 1e-20
    
    # defaults
    if gdat.exprtype == 'plnk':
        gdat.freqexpr = gdat.freqplnk
        gdat.numbfreqexpr = gdat.freqexpr.size
        
    if gdat.exprtype == 'pixi':
        if gdat.minmfreqexpr == None or gdat.maxmfreqexpr == None or gdat.numbfreqexpr == None:
            if gdat.minmfreqexpr == None:
                gdat.minmfreqexpr = gdat.minmfreqpixi
            if gdat.maxmfreqexpr == None:
                gdat.maxmfreqexpr = gdat.maxmfreqpixi
            if gdat.numbfreqexpr == None:
                gdat.numbfreqexpr = gdat.numbfreqpixi
            gdat.freqexpr = logspace(log10(gdat.minmfreqexpr), log10(gdat.maxmfreqexpr), gdat.numbfreqexpr)
        else:
            gdat.freqexpr = gdat.freqpixi
        if fluxstdvinst == None:
            gdat.fluxstdvinst = gdat.fluxstdvinstpixi
        if fluxstdvfrac == None:
            gdat.fluxstdvfrac = gdat.fluxstdvfracpixi
        
    # get expected instrumental uncertainty for PIXIE
    # temp
    if fluxstdvinst == None:
        path = gdat.pathdata + 'pixifluxstdv.csv'
        gdat.fluxstdvinstfreq = loadtxt(path)
        gdat.fluxstdvinstfreq = interp1d(gdat.fluxstdvinstfreq[:, 0] * 1e9, gdat.fluxstdvinstfreq[:, 1] * 1e26)(gdat.freqexpr)

    # get the time stamp
    gdat.strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # get the Healpix grid for Planck maps
    numbsideplnk = 2048
    lgalheal, bgalheal, numbpixlplnk, apixplnk = tdpy.util.retr_healgrid(numbsideplnk) 

    indxpixlrofi = where((lgalheal - lgalcntr < maxmgang) & (bgalheal - bgalcntr < maxmgang))[0]
    
    # run tag
    rtag = retr_rtag(gdat)

    # paths
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'cmbr_dist/', 'cmbr_dist/', gdat.strgtimestmp + '_' + rtag)
    
    # Boolean flag to indicate whether the posterior of intermediate variables will be stored
    gdat.savepost = False
    
    # plotting settings
    gdat.plotsize = 6

    # constants
    ## fundmental
    gdat.velolght = 3e8 # [m/s]
    mp2m = 3.1e22 # [Mpc/m]
    yr2s = 364. * 3600. * 24. # [year/s]
    gdat.consplnk = 6.63e-34 # [J s]
    gdat.consbolt = 1.38e-23 # [J/K]
    gdat.freqpivt = 1e12 # [Hz]
    gdat.tempcmbrconc = 2.725 # Planck concordance model temperature of the CMB today [K]
    gdat.consboltnatu = 8.6173e-5 # Boltzmann constant [eV/K]
    massprot = 9.38e8 # [eV]
    ## distortion
    alph = 1.401
    redsdism = 2e6
    ## cosmological constants
    edencrit = 4e9 # [eV/m^3]
    hubb = 0.68
    omegbmat = 0.049
    omegdmat = 0.26
    omegmatt = omegbmat + omegdmat
    omegradi = 4.8e-5
    omegdene = 0.69

    # axes
    ## frequencies at which the model will be evaluated
    numbfreqmodl = 100
    gdat.minmfreqmodl = 1e9
    gdat.maxmfreqmodl = 1e13
    gdat.freqmodl = logspace(log10(gdat.minmfreqmodl), log10(gdat.maxmfreqmodl), numbfreqmodl) # [Hz]
    gdat.freqmodlscal = gdat.freqmodl * 1e-9
    ## redshift
    numbreds = 100
    minmreds = 1e2
    maxmreds = 1e10
    gdat.meanreds = logspace(log10(minmreds), log10(maxmreds), numbreds)
    gdat.numbredsplot = 10
    gdat.indxredsplot = [k * numbreds / gdat.numbredsplot for k in range(gdat.numbredsplot)]
    diffreds = gdat.meanreds[1:] - gdat.meanreds[:-1]
    gdat.meanredsextd = logspace(2., 10., numbreds)
    ## time
    gdat.timehubb = 4.5e17 * (0.27 * (1. + gdat.meanreds)**3 + 9.2e-5 * (1. + gdat.meanreds)**4.)**(-0.5)
    gdat.time = zeros(numbreds)
    for c in range(numbreds-1):
        gdat.time[c] = trapz(gdat.timehubb[c:] / (1. + gdat.meanreds[c:]), gdat.meanreds[c:])
    thertempantntemp = (1.76e-11 * gdat.freqmodl)**2 * exp(1.76e-11 * gdat.freqmodl) / (exp(1.76e-11 * gdat.freqmodl) - 1.)**2
    antntempflux = 0.0307 * (gdat.freqmodl / 1e9)**2
    # scaled frequency
    gdat.sfrqconc = gdat.consplnk * gdat.freqmodl / gdat.consbolt / gdat.tempcmbrconc
   
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
    
    gdat.fluxcmbrconc = retr_fluxcmbr(gdat, gdat.freqmodl, gdat.tempcmbrconc)
    gdat.fluxtdisgrenconc, gdat.fluxydisgrenconc, gdat.fluxmdisgrenconc, gdat.fluxfdisgrenconc, \
                                        gdat.fluxodisgrenconc = retr_fluxgren(gdat, gdat.freqmodl, gdat.tempcmbrconc, disttype='full')
        
    if makeplot:
        plot_pure_dist(gdat)
        plot_grnf(gdat)
    
    # plot settings
    gdat.alphsamp = 0.3

    # sampler setup
    gdat.numbproc = numbproc
    gdat.indxproc = arange(gdat.numbproc)
    
    # temp
    optiprop = False
    #optiprop = True
    
    # parameters
    datapara = retr_datapara()
    mocksampvarb = datapara.true
    numbpara = len(datapara.name)
    indxpara = arange(numbpara)
    
    # initial state
    thissampvarb = copy(mocksampvarb)
    thissamptemp = tdpy.mcmc.cdfn_samp(thissampvarb, datapara)
    thissamp = empty((gdat.numbproc, numbpara))
    for k in gdat.indxproc:
        thissamp[k, :] = thissamptemp

    
    if gdat.datatype == 'mock':
         
        gdat.mockfluxtotl, gdat.mockfluxcmbr, gdat.mockfluxdustcold, gdat.mockfluxdustwarm, gdat.mockfluxsync, \
                gdat.mockfluxfree, gdat.mockfluxydis, gdat.mockfluxdeca = retr_flux(gdat, gdat.freqmodl, mocksampvarb)
        gdat.mockfluxtotlintp = interp1d(gdat.freqmodl, gdat.mockfluxtotl)(gdat.freqexpr)
                
        if gdat.exprtype == 'pixi':
            gdat.datafluxstdv = gdat.mockfluxtotlintp * gdat.fluxstdvfrac + gdat.fluxstdvinstfreq * gdat.fluxstdvinst
        if gdat.exprtype == 'plnk':
            gdat.datafluxstdv = gdat.mockfluxtotlintp * gdat.fluxstdvfrac + gdat.fluxstdvinst
            
        gdat.dataflux = gdat.mockfluxtotlintp + gdat.datafluxstdv * randn(gdat.datafluxstdv.size)
        
        print 'gdat.dataflux'
        print gdat.dataflux

        # temp
        #fluxegbl = retr_egbl(gdat.freqexpr)
        #gdat.dataflux += fluxegbl
        plot_flux(gdat, plottype='mock')
    else:
        if gdat.exprtype == 'pixi':
            raise Exception('PIXIE has not been launched!')
        if gdat.exprtype == 'plnk':
            gdat.datafluxplnk, gdat.datafluxstdvplnk = tdpy.util.retr_mapsplnkfreq(indxpixlrofi, numbsideoutp=4)
            plot_plnk(gdat)

    
    gdat.numbfreqexpr = gdat.freqexpr.size

    chan = tdpy.mcmc.init(retr_llik, datapara, numbproc=gdat.numbproc, initsamp=thissamp, numbswep=numbswep, numbburn=numbburn, factthin=factthin, rtag=rtag, \
                                                  gdatextr=gdat, factpropeffi=3., verbtype=gdat.verbtype, pathdata=gdat.pathdata, pathimag=gdat.pathimag, optiprop=optiprop)
    listsampvarb = chan[0]
    listsamp = chan[1]
    listsampcalc = chan[2]
    listllik = chan[3]
    listaccp = chan[4]
    
    gdat.numbsamp = listsamp.shape[0]
    
    listflux = zeros((gdat.numbsamp, gdat.numbfreqexpr))
    
    statpara = zeros((numbpara, 3))
    statpara[:, 0] = percentile(listsampvarb, 10., axis=0)
    statpara[:, 1] = percentile(listsampvarb, 50., axis=0)
    statpara[:, 2] = percentile(listsampvarb, 90., axis=0)

    if makeplot:

        indxparaself = where(datapara.scal == 'self')[0]
        indxparalogt = where(datapara.scal == 'logt')[0]

        listsampvarbtran = empty_like(listsampvarb)
        strgparatran = empty(numbpara, dtype=object)
        mocksampvarbtran = zeros_like(mocksampvarb)
        datapara.scaltran = empty(numbpara, dtype=object)
        datapara.scaltran[:] = 'self'

        listsampvarbtran[:, indxparaself] = listsampvarb[:, indxparaself] / mocksampvarb[None, indxparaself] - 1.
        listsampvarbtran[:, indxparalogt] = log10(listsampvarb[:, indxparalogt] / mocksampvarb[None, indxparalogt])
        listsampvarbtran[:, :numbpara-3] *= 1e6

        strgparatran[:] = ''
        strgparatran[:numbpara-3] += r'$10^6 \times$ '
        strgparatran[indxparaself] += '(' + datapara.labl[indxparaself] + ' / ' + datapara.labl[indxparaself] + r'$^{mock}$ - 1)'
        strgparatran[indxparalogt] += 'log(' + datapara.labl[indxparalogt] + ' / ' + datapara.labl[indxparalogt] + r'$^{mock}$)' 
        
        if gdat.savepost:
            listfluxtotl = empty((gdat.numbsamp, numbfreqmodl))
            listfluxcmbr = empty((gdat.numbsamp, numbfreqmodl))
            listfluxdustcold = empty((gdat.numbsamp, numbfreqmodl))
            listfluxdustwarm = empty((gdat.numbsamp, numbfreqmodl))
            listfluxsync = empty((gdat.numbsamp, numbfreqmodl))
            listfluxfree = empty((gdat.numbsamp, numbfreqmodl))
            listfluxydis = empty((gdat.numbsamp, numbfreqmodl))
            listfluxdeca = empty((gdat.numbsamp, numbfreqmodl))
            listresiflux = empty((gdat.numbsamp, gdat.numbfreqexpr))
            for k in range(gdat.numbsamp):
                listfluxtotl[k, :] = listsampcalc[0][k]
                listfluxcmbr[k, :] = listsampcalc[1][k]
                listfluxdustcold[k, :] = listsampcalc[2][k]
                listfluxdustwarm[k, :] = listsampcalc[3][k]
                listfluxsync[k, :] = listsampcalc[4][k]
                listfluxfree[k, :] = listsampcalc[5][k]
                listfluxydis[k, :] = listsampcalc[6][k]
                listfluxdeca[k, :] = listsampcalc[7][k]
                listresiflux[k, :] = listsampcalc[8][k]

            plot_flux(gdat, plottype='post')

    return statpara
    
    
def retr_rtag(gdat):
    
    if gdat.exprtype == 'plnk':
        rtag = 'plnk'
    else:
        rtag = '%d_%03.1f_%03.1f_%03.1f_%03.1f' % (gdat.numbfreqexpr, log10(gdat.minmfreqexpr), log10(gdat.maxmfreqexpr), log10(gdat.fluxstdvinst * 1e3), -log10(gdat.fluxstdvfrac))
 
    return rtag


def cnfg_arca():

    # ARCADE frequency axis
    cnfg['freqarca'] = array([3., 5., 8., 10., 30., 90.]) * 1e9
    cnfg['fluxarca'] = 24.1 * (freqarca / 0.31e9)**(-2.6)


def cnfg_plnk_mock():

    statpara = init( \
                    exprtype='plnk', \
                   )
    

def cnfg_plnk_expr():
    
    statpara = init( \
                    datatype='inpt', \
                    maxmgang=10., \
                    lgalcntr=0., \
                    bgalcntr=0., \
                    exprtype='plnk', \
                    inclcmbrmono=False, \
                   )
    

def cnfg_pixi():
    
    
    statpara = init( \
                   )
    

def cnfg_pixi_mock_stdv():
    
    numbswep = 10000
        
    datapara = retr_datapara()
    numbpara = len(datapara.labl)
    
    strgpertpara = [r'$\nu_{min}$', r'$\nu_{max}$', r'$N_\nu$', r'$\sigma$', r'$\sigma_f$']
    numbpertpara = 4
    numbpert = 1
    listminmfreqexpr = logspace(log10(3e9), log10(3e11), numbpert)
    listmaxmfreqexpr = logspace(log10(1.5e12), log10(6e12), numbpert)
    listnumbfreqexpr = linspace(100, 700, numbpert)
    listfluxstdvinst = logspace(log10(5e0), log10(5e4), numbpert)
    listfluxstdvfrac = logspace(-10., -6., numbpert)
    
    statparagrid = zeros((numbpertpara, numbpert, numbpara, 3))
    for k in range(numbpertpara):
        
        minmfreqexpr = 3e10 # [Hz]
        maxmfreqexpr = 3e12 # [Hz]
        numbfreqexpr = 400
        fluxstdvinst = 5e0 # [Jy/sr]
        fluxstdvfrac = 1e-10
        
        for l in range(numbpert):
            
            if k == 0:
                gdat.minmfreqexpr = arrygdat.minmfreqexpr[l]
            if k == 1:
                gdat.maxmfreqexpr = arrygdat.maxmfreqexpr[l]
            if k == 2:
                gdat.numbfreqexpr = arrygdat.numbfreqexpr[l]
            if k == 3:
                gdat.fluxstdvinst = arrygdat.fluxstdvinst[l]
            if k == 4:
                gdat.fluxstdvfrac = arrygdat.fluxstdvfrac[l]

            statparagrid[k, l, :, :] = init( \
                                            numbswep=numbswep, \
                                            inclcmbrmono=True, \
                                            numbfreqexpr=numbfreqexpr, \
                                            minmfreqexpr=minmfreqexpr, \
                                            maxmfreqexpr=maxmfreqexpr, \
                                            fluxstdvinst=fluxstdvinst, \
                                            fluxstdvfrac=fluxstdvfrac \
                                           )

    for k in range(numbpertpara):
        
        path = gdat.pathimag + 'stdv%d.pdf' % k
        figr, axgr = plt.subplots(numbpara / 2, 2, figsize=(14, numbpara * 3))

        if k == 0:
            xaxi = arrygdat.minmfreqexpr
        if k == 1:
            xaxi = arrygdat.maxmfreqexpr
        if k == 2:
            xaxi = arrygdat.numbfreqexpr
        if k == 3:
            xaxi = arrygdat.fluxstdvinst
        if k == 4:
            xaxi = arrygdat.fluxstdvfrac
            
        for a, axrw in enumerate(axgr):
            for b, axis in enumerate(axrw):

                m = 2 * a + b
                
                axis.plot(xaxi, statparagrid[k, :, m, 0])
                axis.plot(xaxi, statparagrid[k, :, m, 1])
                axis.plot(xaxi, statparagrid[k, :, m, 2])
                axis.set_xlabel(strgpertpara[k])
                axis.set_ylabel(datapara.labl[m] + ' ' + unitpara[m])
                if datapara.scal[m] == 'logt':
                    axis.set_yscale('log')
                if k != 2:
                    axis.set_xscale('log')
                    
                axis.set_xlim([amin(xaxi), amax(xaxi)])
                axis.set_ylim([amin(statparagrid[k, :, m, :]), amax(statparagrid[k, :, m, :])])

        plt.savefig(path)
        plt.close(figr)


globals().get(sys.argv[1])()
