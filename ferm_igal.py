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

    ## templates
    gdat.fluxback = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:

        if c == 0:
            strg = 'isotflux'
        if c == 1:
            strg = 'fdfmflux'
        if c == 2:
            strg = 'plnkdust'
        if c == 3:
            strg = 'wisestar'
        if c == 4:
            strg = 'finkdust'
        if c == 5:
            strg = 'darktemp'

        # temp -- ROI should be fixed at 40 X 40 degree^2
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + strg + '.fits'
        if os.path.isfile(path):
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
            fluxback = tdpy.util.smth_ferm(fluxback, gdat.meanenerfull, gdat.indxevttfull)
            # temp
            fluxback[where(fluxback < 0.)] = 0.

            pf.writeto(path, fluxback, clobber=True)

        # take only the energy bins, spatial pixels and event types of interest
        fluxback = fluxback[indxdatacubefilt]
        indxdatacubetemp = meshgrid(array([c]), gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        gdat.fluxback[indxdatacubetemp] = fluxback
        
    gdat.pathbase = os.environ["FERM_IGAL_DATA_PATH"]
    gdat.pathplot = gdat.pathbase + '/imag/%s/' % gdat.rtag
    cmnd = 'mkdir -p ' + gdat.pathplot
    os.system(cmnd)

    # plot the input spatial templates
    for c in gdat.indxback:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                path = gdat.pathplot + 'fluxback_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_heal(path, gdat.fluxback[c, i, :, m], indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
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
            tdpy.util.plot_heal(path, gdat.dataflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            for c in gdat.indxback:
                path = gdat.pathplot + 'medimodlflux_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_heal(path, medimodlflux[c, i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = gdat.pathplot + 'medimodlfluxtotl_%d%d.pdf' % (i, m)
            tdpy.util.plot_heal(path, medimodlfluxtotl[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = gdat.pathplot + 'mediresiflux_%d%d.pdf' % (i, m)
            tdpy.util.plot_heal(path, mediresiflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
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
              maxmnormback=array([5., 5.]), \
              minmnormback=array([.2, .2]), \
              strgback=['isotflux.fits', 'fdfmflux.fits'], \
              strgexpo=strgexpo, \
              datatype='inpt', \
              strgexpr=strgexpr, \
             )
    
    
def pcat_ferm_mock_igal_brok():
     
    indxenerincl = arange(1, 4)
    indxevttincl = arange(2, 4)
    numbener = indxenerincl.size

    mockfdfnslop = array([2.5])
    listfdfnbrek = array([1e-10, 3e-10, 1e-9, 3e-9, 1e-8])
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
                       maxmnormback=array([2., 2.]), \
                       minmnormback=array([0.5, 0.5]), \
                       strgexpo='fermexpo_cmp0_igal.fits', \
                       strgback=['isotflux.fits', 'fdfmflux.fits'], \
                       datatype='mock', \
                       mockfdfntype='brok', \
                       mocknumbpnts=array([300]), \
                       numbsideheal=256, \
                       mockfdfnslop=mockfdfnslop, \
                       mockfdfnsloplowr=mockfdfnsloplowr, \
                       mockfdfnslopuppr=mockfdfnslopuppr, \
                       mockfdfnbrek=array([fdfnbrek]), \
                       mocknormback=ones((2, numbener)), \
                       pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                      )


def pcat_ferm_mock_igal():
     
    indxevttincl = arange(2, 4)
    indxenerincl = arange(1, 4)
    numbener = indxenerincl.size

    minmflux = 5e-11
    maxmflux = 3e-7
    mockfdfnslop = array([1.9])
      
    pcat.main.init( \
                   psfntype='doubking', \
                   numbswep=5000, \
                   randinit=False, \
                   maxmgang=20., \
                   indxevttincl=indxevttincl, \
                   indxenerincl=indxenerincl, \
                   numbsideheal=256, \
                   mocknumbpnts=array([800]), \
                   maxmnumbpnts=array([1200]), \
                   minmflux=minmflux, \
                   maxmflux=maxmflux, \
                   mocknormback=ones((2, numbener)), \
                   maxmnormback=array([2., 2.]), \
                   mockfdfnslop=mockfdfnslop, \
                   minmnormback=array([0.5, 0.5]), \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits', 'fdfmflux.fits'], \
                   pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                   regitype='igal', \
                   datatype='mock', \
                  )


def cnfg_nomi():
    
    init( \
         verbtype=1, \
         makeplot=True, \
         strgback=['isotflux', 'fdfmflux'], \
        )


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass
    cnfg_nomi()

