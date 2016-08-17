from __init__ import *

def retr_modlflux(gdat, sampvarb):

    norm = sampvarb.reshape((gdat.numbback, gdat.numbener))
    modlflux = norm[:, :, None, None] * gdat.fluxback
   
    return modlflux


def make_maps_rec7_back():
    
    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec7']
    gdat.enertype = ['back']
    tdpy.util.make_maps_main(gdat, os.environ["FERM_IGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec7', 'back', 'igal', os.environ["FERM_IGAL_DATA_PATH"], 256, 'tim0')


def make_maps_rec8_back():
    
    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec8']
    gdat.enertype = ['back']
    tdpy.util.make_maps_main(gdat, os.environ["FERM_IGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec8', 'back', 'igal', os.environ["FERM_IGAL_DATA_PATH"], 256, 'tim0')


def merg_maps():

    numbside = 256
    scalwght = 'self'
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])[2:4]

    pathplot = os.environ["FERM_IGAL_DATA_PATH"]
    # get input maps
    ## Planck map
    path = pathplot + '/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
    maps = pf.getdata(path, 1)['RADIANCE']
    mapsplnk = hp.ud_grade(maps, numbside, order_in='NESTED', order_out='RING')
    mapsplnk -= mean(mapsplnk)
    mapsplnk /= std(mapsplnk)
    almcplnk = hp.map2alm(mapsplnk)

    ## Fermi Diffuse Model
    # temp
    mapsfdfm = tdpy.util.retr_fdfm(binsener).flatten()
    mapsfdfm -= mean(mapsfdfm)
    mapsfdfm /= std(mapsfdfm)
    almcfdfm = hp.map2alm(mapsfdfm)

    # get the merged map
    mpol = arange(2. * numbside)
    maxmmpol = amax(mpol)

    wght = empty_like(almcplnk)

    vari = 100.**2
    wghtsing = 1. / sqrt(2. / pi / vari) * exp(-0.5 * mpol * (mpol + 1.) / vari**2)
    mpolgrid, temp = hp.Alm.getlm(lmax=maxmmpol)
    if scalwght == 'self':
        for l in mpol:
            wght[where(mpolgrid == l)] = wghtsing[l]
    
    # plot the weight
    figr, axis = plt.subplots()
    axis.loglog(mpol, wght, label='FDM')
    axis.loglog(mpol, 1. - wght, label='Planck')
    axis.set_ylabel('$w_l$')
    axis.set_xlabel('$l$')
    axis.legend()
    
    path = pathplot + '/wght.pdf'
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)
    
    
    
    almcoutp = almcfdfm * wght + almcplnk * (1. - wght)
    mapsmerg = hp.alm2map(almcoutp, numbside, verbose=False)
                
    rtag = '%04d_%s' % (numbside, scalwght)

    for plotigal in [True, False]:

        if plotigal:
            minmlgal = -20.
            maxmlgal = 20.
            minmbgal = -20.
            maxmbgal = 20.
            strg = rtag + 'igal'
        else:
            minmlgal = -180.
            maxmlgal = 180.
            minmbgal = -90.
            maxmbgal = 90.
            strg = rtag + 'full'
           
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/mapsfdfm%s.pdf' % strg
        tdpy.util.plot_maps(path, mapsfdfm, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
        
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/mapsplnk%s.pdf' % strg
        tdpy.util.plot_maps(path, mapsplnk, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
        
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/mapsmerg%s.pdf' % strg
        tdpy.util.plot_maps(path, mapsmerg, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
        
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/mapsresifdfm%s.pdf' % strg
        tdpy.util.plot_maps(path, mapsmerg - mapsfdfm, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, resi=True)
        
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/mapsresiplnk%s.pdf' % strg
        tdpy.util.plot_maps(path, mapsmerg - mapsplnk, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, resi=True)
        

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
        datapara.scal[n] = 'logt'
        if n // gdat.numbener == 0 or n // gdat.numbener == 1:
            datapara.minm[n] = 1e-4
            datapara.maxm[n] = 1e1
        else:
            datapara.minm[n] = 1e-9
            datapara.maxm[n] = 1e-3
        if n // gdat.numbener == 0:
            strg = 'isot'
        if n // gdat.numbener == 1:
            strg = 'fdfm'
        if n // gdat.numbener == 2:
            strg = 'plnk'
        if n // gdat.numbener == 3:
            strg = 'wise'
        if n // gdat.numbener == 5:
            strg = 'fink'
        if n // gdat.numbener == 4:
            strg = 'dark'
        datapara.labl[n] = '$A_{%d}^{%s}$' % (n % gdat.numbener, strg)
        datapara.unit[n] = ''
        datapara.vari[n] = 1e-1
        datapara.true[n] = None
    
    datapara.strg = datapara.labl + ' ' + datapara.unit
    
    return datapara


def retr_rtag(gdat):
    
    gdat.rtag = '%d_%d_%s' % (gdat.numbproc, gdat.numbswep, gdat.datatype)


def defn_gtbn():
    
    numbener = 30
    minmener = 0.1
    maxmener = 100.
    binsener = logspace(log10(minmener), log10(maxmener), numbener + 1)
    lowrener = binsener[:-1]
    upprener = binsener[1:]
    limtener = stack((lowrener, upprener), axis=1)
    path = os.environ["TDPY_DATA_PATH"] + '/gtbndefn_back.dat'
    savetxt(path, limtener, fmt='%10.5g')


def retr_ener(gdat):

    # temp
    gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
    gdat.meanenerfull = sqrt(gdat.binsenerfull[1:] * gdat.binsenerfull[:-1])
    gdat.numbenerfull = gdat.binsenerfull.size - 1
    gdat.indxenerfull = arange(gdat.numbenerfull)
    gdat.binsener = gdat.binsenerfull[gdat.indxenerincl[0]:gdat.indxenerincl[-1] + 2]
    gdat.meanener = sqrt(gdat.binsener[1:] * gdat.binsener[:-1])
    gdat.diffener = gdat.binsener[1:] - gdat.binsener[:-1]
    gdat.numbener = gdat.meanener.size
    gdat.indxener = arange(gdat.numbener)
    gdat.strgbinsener = ['%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]) for i in gdat.indxener]
    

def plot_backspec(gdat, indxpixlmean):
    
    listcolr = ['black', 'b', 'g', 'r', 'm', 'y'][:gdat.numbback+1]
    listlabl = ['Data', 'Isotropic', 'FDM', 'Planck', r'WISE 12$\mu$m', 'NFW'][:gdat.numbback+1]

    figr, axis = plt.subplots()
    
    numbvarb = gdat.numbback + 1
    listydat = empty((numbvarb, gdat.numbener))
    listyerr = zeros((2, numbvarb, gdat.numbener))
    
    listydat[0, :] = mean(sum(gdat.datacnts, 2) / sum(gdat.expo, 2), 1) / gdat.apix / gdat.diffener
    listyerr[:, 0, :] = mean(sqrt(sum(gdat.datacnts, 2)) / sum(gdat.expo, 2), 1) / gdat.apix / gdat.diffener
    for n in gdat.indxback:
        if n == 0 or n == 1:
            listydat[n+1, :] = gdat.postnormback[0, n, :] * mean(mean(gdat.fluxback[n, :, :, :], 1), 1)
        else:
            listydat[n+1, :] = gdat.postnormback[0, n, :]
        listyerr[:, n+1, :] = tdpy.util.retr_errrvarb(gdat.postnormback[:, n, :])
    
    xdat = gdat.meanener
    for k in range(numbvarb):
        ydat = gdat.meanener**2 * listydat[k, :]
        yerr = gdat.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5, color=listcolr[k], label=listlabl[k])

    # Fermi-LAT results
    # temp
    if False:
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
    axis.legend(loc=4, ncol=2)

    path = gdat.pathplot + 'backspec.pdf'
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)


def plot_psec(gdat, mapsplot):

    listlabl = ['Data', 'Isotropic']
    listcolr = ['black', 'b']
    
    listlabl.append('FDM, %s' % gdat.strgbinsener[1])
    listcolr.append('g')
    
    listlabl.extend(['Planck', r'WISE 12$\mu$m', 'NFW'])
    listcolr.extend(['r', 'm', 'y'])

    figr, axis = plt.subplots()
    mpol = arange(3 * gdat.numbside, dtype=int)
    print 'hey'
    for n in range(gdat.numbback):
        print 'n: ', n
        print 'mapsplot[n, :]'
        print amin(mapsplot[n, :]), amax(mapsplot[n, :]), mean(mapsplot[n, :]), std(mapsplot[n, :])
        psec = hp.anafast(mapsplot[n, :])
        print 

        axis.loglog(mpol, mpol * (mpol + 1.) * psec, color=listcolr[n], label=listlabl[n])

    axis.set_ylabel('$l(l+1)C_l$')
    axis.set_xlabel('$l$')
    
    axis.legend(loc=4, ncol=2)

    path = gdat.pathplot + 'psec.pdf'
    plt.tight_layout()
    plt.savefig(path)


def init( \
         numbproc=1, \
         numbswep=None, \
         datatype='inpt', \
         verbtype=1, \
         makeplot=False, \
         strgexpr='fermflux_cmp0_igal.fits', \
         strgexpo='fermexpo_cmp0_igal.fits', \
         strgback=['isotflux', 'fdfmflux', 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 'wssa_sample_1024.fits', 'lambda_sfd_ebv.fits', 'darktemp'], \
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

    boolsmth = False

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

    boolplotpsec = True

    indxdatacubefilt = meshgrid(gdat.indxenerincl, gdat.indxpixlrofi, gdat.indxevttincl, indexing='ij')
    
    # get data structure
    datapara = retr_datapara(gdat)
    
    if gdat.numbswep == None:
        gdat.numbswep = 4 * gdat.numbpara * 10000

    # setup
    retr_rtag(gdat)
    
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

    # power spectrum calculation
    gdat.numbmapsplot = gdat.numbback + 1
    gdat.mapsplot = empty((gdat.numbmapsplot, gdat.numbpixlfull))
    
    gdat.mapsplot[0, gdat.indxpixlrofi] = sum(gdat.datacnts[1, :, :], 1)

    ## templates
    gdat.fluxback = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:

        if boolsmth:
            strg = '_smth'
        else:
            strg = ''
        if c == 0:
            strg = 'isotflux'
        if c == 1:
            strg = 'fdfmflux' + strg
        if c == 2:
            strg = 'plnkdust' + strg
        if c == 3:
            strg = 'wisestar' + strg
        if c == 4:
            strg = 'finkdust' + strg
        if c == 5:
            strg = 'darktemp' + strg

        # temp -- ROI should be fixed at 40 X 40 degree^2
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + strg + '.fits'
        if os.path.isfile(path) or not boolplotpsec:
            fluxback = pf.getdata(path)
        else:
            
            if c == 0:
                fluxbackorig = tdpy.util.retr_isot(gdat.binsenerfull)
            if c == 1:
                fluxbackorig = tdpy.util.retr_fdfm(gdat.binsenerfull) 
            if c == 2:
                pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgback[c]
                fluxbackorig = pf.getdata(pathtemp, 1)['RADIANCE']
                fluxbackorig = hp.ud_grade(fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 3:
                pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgback[c]
                fluxbackorig = pf.getdata(pathtemp, 0)
                fluxbackorig = hp.ud_grade(fluxbackorig, gdat.numbside, order_in='RING', order_out='RING')
            if c == 4:
                pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgback[c]
                fluxbackorig = pf.getdata(pathtemp)['TEMPERATURE']
                fluxbackorig = hp.ud_grade(fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 5:
                fluxbackorig = tdpy.util.retr_nfwp(1., gdat.numbside)

            # normalize the templates
            if c != 0 and c != 1:
                fluxbackorig /= mean(fluxbackorig[gdat.indxpixlrofi])
            
            # smooth the templates
            fluxback = empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            for m in gdat.indxevttfull:
                if c == 0 or c == 1:
                    fluxback[:, :, m] = fluxbackorig
                else:
                    for i in gdat.indxenerfull:
                        fluxback[i, :, m] = fluxbackorig
            
            if boolsmth:
                fluxback = tdpy.util.smth_ferm(fluxback, gdat.meanenerfull, gdat.indxevttfull)
            # temp
            fluxback[where(fluxback < 0.)] = 0.

            pf.writeto(path, fluxback, clobber=True)

            # load the map to the array whose power spectrum will be calculated
            if c == 0:
                gdat.mapsplot[c+1, :] = fluxbackorig[1, :]
            elif c == 1:
                gdat.mapsplot[c+1, :] = fluxbackorig[1, :]
            else:
                gdat.mapsplot[c+1, :] = fluxbackorig
    
        # take only the energy bins, spatial pixels and event types of interest
        fluxback = fluxback[indxdatacubefilt]
        indxdatacubetemp = meshgrid(array([c]), gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        gdat.fluxback[indxdatacubetemp] = fluxback

    gdat.pathbase = os.environ["FERM_IGAL_DATA_PATH"]
    gdat.pathplot = gdat.pathbase + '/imag/%s/' % gdat.rtag
    cmnd = 'mkdir -p ' + gdat.pathplot
    os.system(cmnd)

    # plot the power spectra
    plot_psec(gdat, gdat.mapsplot)
    
    # plot the input spatial templates
    for c in gdat.indxback:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                path = gdat.pathplot + 'fluxback_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_maps(path, gdat.fluxback[c, i, :, m], indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                                                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
            
    initsamp = rand(gdat.numbproc * gdat.numbpara).reshape((gdat.numbproc, gdat.numbpara))

    numbplotside = gdat.numbpara
    chan = tdpy.mcmc.init(retr_llik, datapara, numbproc=gdat.numbproc, numbswep=gdat.numbswep, initsamp=initsamp, gdatextr=gdat, \
                verbtype=gdat.verbtype, pathbase=gdat.pathbase, rtag=gdat.rtag, numbplotside=numbplotside)
    
    listsampvarb, listsamp, listsampcalc, listllik, listaccp, listindxparamodi, propeffi, levi, info, gmrbstat = chan
    numbsamp = listsamp.shape[0]

    gdat.medisampvarb = percentile(listsampvarb, 50., axis=0)
    medimodlflux = retr_modlflux(gdat, gdat.medisampvarb)
    medimodlfluxtotl = sum(medimodlflux, axis=0)
    mediresiflux = gdat.dataflux - medimodlfluxtotl

    gdat.postsampvarb = tdpy.util.retr_postvarb(listsampvarb)

    gdat.postnormback = gdat.postsampvarb.reshape((3, gdat.numbback, gdat.numbener))

    for i in gdat.indxener:
        for m in gdat.indxevtt:
            
            path = gdat.pathplot + 'dataflux_%d%d.pdf' % (i, m)
            tdpy.util.plot_maps(path, gdat.dataflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            for c in gdat.indxback:
                path = gdat.pathplot + 'medimodlflux_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_maps(path, medimodlflux[c, i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = gdat.pathplot + 'medimodlfluxtotl_%d%d.pdf' % (i, m)
            tdpy.util.plot_maps(path, medimodlfluxtotl[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = gdat.pathplot + 'mediresiflux_%d%d.pdf' % (i, m)
            tdpy.util.plot_maps(path, mediresiflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True, resi=True)
    
    # make plots of the spectra of spatially averaged background components
    plot_backspec(gdat, gdat.indxpixlrofi)
    
    # temp
    #indxpixlmean = where(abs(gdat.bgalheal) < 2.)[0]
    #plot_backspec(gdat, indxpixlmean)
    #indxpixlmean = where((abs(gdat.lgalheal) < 5.) & (abs(gdat.bgalheal) < 5.))[0]
    #plot_backspec(gdat, indxpixlmean)


def pcat_ferm_expr_igal(strgexpr='fermflux_cmp0_igal.fits', strgexpo='fermexpo_cmp0_igal.fits'):
    
    pcat.main.init( \
              psfntype='doubking', \
              numbswep=2000000, \
              randinit=False, \
              maxmgang=20., \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(2, 4), \
              minmflux=3e-11, \
              maxmflux=3e-6, \
              pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
              regitype='igal', \
              strgback=['isotflux.fits', 'fdfmflux.fits'], \
              strgexpo=strgexpo, \
              datatype='inpt', \
              strgexpr=strgexpr, \
             )
    
    
def pcat_ferm_mock_igal_brok_arry():
     
    indxenerincl = arange(1, 4)
    indxevttincl = arange(2, 4)
    numbener = indxenerincl.size

    listfdfnbrek = array([1e-10, 3e-10, 1e-9, 3e-9, 1e-8])
    listfdfnsloplowr = array([1.9, 2.2, 2.8, 3.1, 3.4])
    mockfdfnsloplowr = array([1.])
    mockfdfnslopuppr = array([2.5])


    for fdfnbrek in listfdfnbrek:
        pcat.main.init(psfntype='doubking', \
                       numbswep=400000, \
                       numbburn=0, \
                       randinit=False, \
                       boolproppsfn=False, \
                       maxmgang=20., \
                       fdfntype='brok', \
                       indxenerincl=indxenerincl, \
                       indxevttincl=indxevttincl, \
                       maxmnumbpnts=array([600]), \
                       minmflux=3e-11, \
                       maxmflux=1e-7, \
                       regitype='igal', \
                       strgexpo='fermexpo_cmp0_igal.fits', \
                       strgback=['isotflux.fits', 'fdfmflux.fits'], \
                       datatype='mock', \
                       mockfdfntype='brok', \
                       mocknumbpnts=array([300]), \
                       mocknumbsideheal=256, \
                       mockfdfnslop=mockfdfnslop, \
                       mockfdfnsloplowr=mockfdfnsloplowr, \
                       mockfdfnslopuppr=mockfdfnslopuppr, \
                       mockfdfnbrek=array([fdfnbrek]), \
                       mocknormback=ones((2, numbener)), \
                       pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                      )


def pcat_ferm_mock_igal_brok():
     
    indxevttincl = arange(2, 4)
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = 3e-11
    maxmflux = 3e-7
    mockfdfnbrek = array([1e-8])
    mockfdfnsloplowr = array([1.6])
    mockfdfnslopuppr = array([2.6])
    
    pcat.main.init( \
                   psfntype='doubking', \
                   numbswep=10000, \
                   randinit=False, \
                   maxmgang=20., \
                   indxevttincl=indxevttincl, \
                   indxenerincl=indxenerincl, \
                   fdfntype='brok', \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits', 'fdfmflux.fits'], \
                   pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                   regitype='igal', \
                   
                   maxmnumbpnts=array([400]), \
                   minmflux=minmflux, \
                   maxmflux=maxmflux, \

                   datatype='mock', \
                   numbsideheal=256, \
                   mocknormback=ones((2, numbener)), \
                   mocknumbpnts=array([400]), \
                   mockspatdisttype=['unif', 'disc', 'gang'], \
                   mockfdfntype='brok', \
                   mockfdfnbrek=mockfdfnbrek, \
                   mockfdfnsloplowr=mockfdfnsloplowr, \
                   mockfdfnslopuppr=mockfdfnslopuppr, \
                   mocksdfnstdv=array([.5]), \
                   mocksdfnmean=array([2.]), \
                  )


def pcat_ferm_mock_igal():
     
    indxevttincl = arange(2, 4)
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = 3e-11
    maxmflux = 3e-7
    mockfdfnslop = array([2.6, 2.6, 3.5])
      
    pcat.main.init( \
                   psfntype='doubking', \
                   numbswep=100, \
                   randinit=False, \
                   indxevttincl=indxevttincl, \
                   indxenerincl=indxenerincl, \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits', 'fdfmflux.fits'], \
                   pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                   regitype='igal', \
                   
                   maxmnumbpnts=array([200, 200, 400]), \
                   maxmgang=20., \
                   minmflux=minmflux, \
                   maxmflux=maxmflux, \
                   fdfntype='powr', \
                   sdfnstdv=array([.5, .5, .5]), \
                   sdfnmean=array([2., 2., 2.]), \
                   
                   datatype='mock', \
                   numbsideheal=256, \
                   mocknormback=ones((2, numbener)), \
                   mocknumbpnts=array([100, 100, 200]), \
                   mockspatdisttype=['unif', 'disc', 'gang'], \
                   mockfdfntype='powr', \
                   mockfdfnslop=mockfdfnslop, \
                   mocksdfntype='gaus', \
                   mocksdfnstdv=array([.5, .5, .5]), \
                   mocksdfnmean=array([2., 2., 2.]), \
                  )


def cnfg_nomi():
    
    init( \
         verbtype=1, \
         makeplot=True, \
         #strgback=['isotflux', 'fdfmflux'], \
        )


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass
    cnfg_nomi()

