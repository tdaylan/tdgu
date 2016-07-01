from commimpt import *

def retr_modlflux(gdat, sampvarb):

    norm = sampvarb.reshape((gdat.numbback, gdat.numbener))
    modlflux = norm[:, :, None, None] * gdat.fluxback
   
    if False:
        print 'retr_modlflux'
        print 'norm'
        print norm
        print 'modlflux'
        print amin(modlflux[0, :, :, :]), amax(modlflux[0, :, :, :])
        print amin(modlflux[1, :, :, :]), amax(modlflux[1, :, :, :])
        print 'modlfluxtotl'
        print amin(modlfluxtotl), amax(modlfluxtotl)
        print 

    return modlflux


def retr_llik(sampvarb, gdat):

    modlflux = retr_modlflux(gdat, sampvarb)
    modlfluxtotl = sum(modlflux, axis=0)
    modlcnts = modlfluxtotl * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    llik = sum(gdat.datacnts * log(modlcnts) - modlcnts)
    sampcalc = []


    return llik, sampcalc


def retr_datapara(gdat):
    
    gdat.numbpara = gdat.numbener * gdat.numbback
    gdat.indxpara = arange(gdat.numbpara)

    datapara = tdpy.util.gdatstrt()

    datapara.indx = dict()
    datapara.minm = zeros(gdat.numbpara)
    datapara.maxm = zeros(gdat.numbpara)
    datapara.name = empty(gdat.numbpara, dtype=object)
    datapara.scal = empty(gdat.numbpara, dtype=object)
    datapara.labl = empty(gdat.numbpara, dtype=object)
    datapara.unit = empty(gdat.numbpara, dtype=object)
    datapara.vari = zeros(gdat.numbpara)
    datapara.true = zeros(gdat.numbpara)

    for n in gdat.indxpara:
        datapara.indx['norm%04d' % n] = n
        datapara.name[n] = 'norm%04d' % n
        datapara.minm[n] = 1e-3
        datapara.maxm[n] = 1e1
        datapara.scal[n] = 'logt'
        if n // gdat.numbener == 0:
            strg = 'isot'
        if n // gdat.numbener == 1:
            strg = 'plnk'
        if n // gdat.numbener == 2:
            strg = 'wise'
        datapara.labl[n] = '$A_{%d}^{%s}$' % (n % gdat.numbener, strg)
        datapara.unit[n] = ''
        datapara.vari[n] = 1e-1
        datapara.true[n] = None
    
    datapara.strg = datapara.labl + ' ' + datapara.unit
    
    return datapara


def retr_rtag(gdat):
    
    gdat.rtag = '%d_%d_%s' % (gdat.numbproc, gdat.numbswep, gdat.datatype)


def retr_ener(gdat):

    gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
    gdat.meanenerfull = sqrt(gdat.binsenerfull[1:] * gdat.binsenerfull[:-1])
    gdat.numbenerfull = gdat.binsenerfull.size - 1
    gdat.indxenerfull = arange(gdat.numbenerfull)

    gdat.binsener = gdat.binsenerfull[gdat.indxenerincl[0]:gdat.indxenerincl[-1] + 2]
    gdat.meanener = sqrt(gdat.binsener[1:] * gdat.binsener[:-1])
    gdat.diffener = gdat.binsener[1:] - gdat.binsener[:-1]
    gdat.numbener = gdat.meanener.size
    gdat.indxener = arange(gdat.numbener)
    

def plot_backspec(gdat, indxpixlmean):
    
    listlinestyl = [':', '--', '-.', '-']
    listcolr = ['b', 'g', 'r', 'black']
    listlabl = ['Isotropic', 'Planck Dust Radiance', 'NFW', 'Data']

    figr, axis = plt.subplots()
    
    numbvarb = gdat.numbback + 1
    listydat = empty((numbvarb, gdat.numbener))
    listyerr = zeros((2, numbvarb, gdat.numbener))
    
    for n in gdat.indxback:
        listydat[n, :] = postnormback[0, n, :] * gdat.backfluxmean[n]
        listyerr[:, n, :] = retr_errrvarb(postnormback[:, n, :]) * gdat.backfluxmean[n]
    listydat[-1, :] = gdat.datafluxmean
    listyerr[:, -1, :] = retr_errrvarb(postpntsfluxmean)
    
    xdat = gdat.meanener
    for k in range(numbvarb):
        ydat = gdat.meanener**2 * listydat[k, :]
        yerr = gdat.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5, ls=listlinestyl[k], color=listcolr[k], label=listlabl[k])

    # Fermi-LAT results
    # temp
    if False and gdat.trueinfo:
        if gdat.datatype == 'mock':
            pass
        else:
            if gdat.exprtype == 'ferm':
                listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
                listmrkr = ['o', 's', 'p', '*', 'D', '^']
                listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
                listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
                for k, name in enumerate(listname):
                    path = os.environ["PCAT_DATA_PATH"] + '/fermspec' + name + '.csv'
                    data = loadtxt(path)
                    enertemp = data[:, 0] # [GeV]
                    fluxtemp = data[:, 1] * 1e-3 # [GeV/cm^2/s/sr]
                    fluxtemp = interp(gdat.meanener, enertemp, fluxtemp)
                    #fluxtemp = interpolate.interp1d(enertemp, fluxtemp)(gdat.meanener)
                    axis.plot(gdat.meanener, fluxtemp, marker=listmrkr[k], color=listcolr[k], label=listlabl[k])

    axis.set_xlim([amin(gdat.binsener), amax(gdat.binsener)])
    axis.set_yscale('log')
    axis.set_xlabel('E [GeV]')
    axis.set_xscale('log')
    axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [GeV/cm$^2$/s/sr]')
    axis.legend()

    path = gdat.pathplot + 'backspec.pdf'
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)


def init( \
         numbproc=1, \
         numbswep=10000, \
         datatype='inpt', \
         verbtype=1, \
         makeplot=False, \
         strgexpr='fermflux_cmp0_igal.fits', \
         strgexpo='fermexpo_cmp0_igal.fits', \
         strgback=['isottemp', 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 'wssa_sample_1024.fits', 'darktemp'], \
         indxenerincl=arange(1, 4), \
         indxevttincl=arange(3, 4), \
         maxmgang=20.
        ):

    # construct the global object
    gdat = tdpy.util.gdatstrt()
    gdat.numbproc = numbproc
    gdat.numbswep = numbswep
    gdat.datatype = datatype
    gdat.verbtype = verbtype
    gdat.makeplot = makeplot
    gdat.strgexpr = strgexpr
    gdat.strgexpo = strgexpo
    gdat.strgback = strgback
    gdat.indxenerincl = indxenerincl
    gdat.indxevttincl = indxevttincl
    gdat.maxmgang = maxmgang

    # setup
    retr_rtag(gdat)
    
    gdat.numbback = len(strgback)
    gdat.indxback = arange(gdat.numbback)
    retr_ener(gdat)

    ## event type
    gdat.evttfull = array([4, 8, 16, 32])
    gdat.numbevttfull = gdat.evttfull.size
    gdat.indxevttfull = arange(gdat.numbevttfull)
    
    gdat.evtt = gdat.evttfull[gdat.indxevttincl]
    gdat.numbevtt = gdat.evtt.size
    gdat.indxevtt = arange(gdat.numbevtt)

    ## pixelization
    gdat.numbside = 256
    gdat.numbpixlfull = gdat.numbside**2 * 12
    gdat.lgalheal, gdat.bgalheal, gdat.numbpixl, gdat.apix = tdpy.util.retr_healgrid(gdat.numbside)
    gdat.indxpixlrofi = where((abs(gdat.lgalheal) < gdat.maxmgang) & (abs(gdat.bgalheal) < gdat.maxmgang))[0]
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.lgalgrid = gdat.lgalheal[gdat.indxpixlrofi]
    gdat.bgalgrid = gdat.bgalheal[gdat.indxpixlrofi]
       
    minmlgal = -gdat.maxmgang
    maxmlgal = gdat.maxmgang
    minmbgal = -gdat.maxmgang
    maxmbgal = gdat.maxmgang

    indxdatacubefilt = meshgrid(gdat.indxenerincl, gdat.indxpixlrofi, gdat.indxevttincl, indexing='ij')
    
    ## data
    if gdat.datatype == 'inpt':
        
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgexpr
        gdat.exprflux = pf.getdata(path)
        
        ### filter
        gdat.dataflux = gdat.exprflux[indxdatacubefilt]

    ## exposure
    path = os.environ["PCAT_DATA_PATH"] + '/' + gdat.strgexpo
    gdat.expo = pf.getdata(path)
    gdat.expo = gdat.expo[indxdatacubefilt]
   
    gdat.datacnts = gdat.dataflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]

    ## templates
    gdat.fluxback = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:

        if c == 0:
            strg = 'isottemp.fits'
        if c == 1:
            strg = 'plnkdust.fits'
        if c == 2:
            strg = 'wisestar.fits'
        if c == 3:
            strg = 'darktemp.fits'

        pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + '/' + strg
        if os.path.isfile(pathtemp):
            fluxback = pf.getdata(pathtemp)
        else:
            if c == 0:
                fluxbacktemp = ones((gdat.numbpixlfull))
            if c == 1:
                path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgback[c]
                fluxbacktemp = pf.getdata(path, 1)['RADIANCE']
                fluxbacktemp = hp.ud_grade(fluxbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 2:
                path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgback[c]
                fluxbacktemp = pf.getdata(path, 0)
                fluxbacktemp = hp.ud_grade(fluxbacktemp, gdat.numbside, order_in='RING', order_out='RING')
            if c == 3:
                fluxbacktemp = tdpy.util.retr_nfwp(1., gdat.numbside)

            fluxback = empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            for i in gdat.indxenerfull:
                for m in gdat.indxevttfull:
                    fluxback[i, :, m] = fluxbacktemp
            maxmmpol = 3 * gdat.numbside - 1
            fluxback = tdpy.util.smth_ferm(fluxback, gdat.meanenerfull, gdat.indxevttfull, maxmmpol=maxmmpol)

            # temp
            fluxback[where(fluxback < 0.)] = 0.

            pf.writeto(pathtemp, fluxback, clobber=True)

        fluxback = fluxback[indxdatacubefilt]
        fluxback *= mean(mean(gdat.dataflux, 1), 1)[:, None, None] / mean(fluxback)
        indxdatacubetemp = meshgrid(array([c]), gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        gdat.fluxback[indxdatacubetemp] = fluxback

    # plot the input spatial templates
    for c in gdat.indxback:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/%s/fluxback_%d%d%d.pdf' % (gdat.rtag, c, i, m)
                tdpy.util.plot_heal(path, gdat.fluxback[c, i, :, m], indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                                                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
            
    # get data structure
    datapara = retr_datapara(gdat)
    
    gdat.pathbase = os.environ["FERM_IGAL_DATA_PATH"]
    gdat.pathplot = gdat.pathbase + '/imag/'

    optiprop = True

    initsamp = rand(gdat.numbproc * gdat.numbpara).reshape((gdat.numbproc, gdat.numbpara))

    numbplotside = gdat.numbpara
    chan = tdpy.mcmc.init(retr_llik, datapara, numbproc=gdat.numbproc, numbswep=gdat.numbswep, initsamp=initsamp, gdatextr=gdat, \
                optiprop=optiprop, verbtype=gdat.verbtype, pathbase=gdat.pathbase, rtag=gdat.rtag, numbplotside=numbplotside)
    
    listsampvarb, listsamp, listsampcalc, listllik, listaccp, listindxparamodi, propeffi, levi, info, gmrbstat = chan
    numbsamp = listsamp.shape[0]

    gdat.medisampvarb = percentile(listsampvarb, 50., axis=0)
    medimodlflux = retr_modlflux(gdat, gdat.medisampvarb)
    medimodlfluxtotl = sum(medimodlflux, axis=0)
    mediresiflux = gdat.dataflux - medimodlfluxtotl

    for i in gdat.indxener:
        for m in gdat.indxevtt:
            
            path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/%s/dataflux_%d%d.pdf' % (gdat.rtag, i, m)
            tdpy.util.plot_heal(path, gdat.dataflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            for c in gdat.indxback:
                path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/%s/medimodlflux_%d%d%d.pdf' % (gdat.rtag, c, i, m)
                tdpy.util.plot_heal(path, medimodlflux[c, i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/%s/medimodlfluxtotl_%d%d.pdf' % (gdat.rtag, i, m)
            tdpy.util.plot_heal(path, medimodlfluxtotl[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/%s/mediresiflux_%d%d.pdf' % (gdat.rtag, i, m)
            tdpy.util.plot_heal(path, mediresiflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True, resi=True)
    
    # make plots of the spectra of spatially averaged background components
    plot_backspec(gdat, gdat.indxpixlrofi)
    
    indxpixlmean = where(abs(gdat.bgalheal) < 2.)[0]
    plot_backspec(gdat, indxpixlmean)
    
    indxpixlmean = where((abs(gdat.lgalheal) < 5.) & (abs(gdat.bgalheal) < 5.))[0]
    plot_backspec(gdat, indxpixlmean)


def cnfg_nomi():
    
    init( \
         numbswep=100000, \
         verbtype=1, \
         makeplot=True, \
        )


if __name__ == '__main__':
    
    cnfg_nomi()

