from __init__ import *

def retr_axes():

    reco = 8
    evtc = 128
    
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    meanener = sqrt(binsener[1:] * binsener[0:-1])
    diffener = binsener[1:] - binsener[0:-1]
    numbener = meanener.size
    indxener = arange(numbener)
    minmener = amin(binsener)
    maxmener = amax(binsener)
    
    global indxevtt
    evtt = array([4, 16, 32, 64])
    numbevtt = evtt.size
    indxevtt = arange(numbevtt)
    
    nside = 256
    numbpixl = nside**2 * 12
    apix = 4. * pi / numbpixl
    
    return reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix


def make_maps_pss7pnts():
    
    gdat = tdpy.util.gdatstrt()
    
    numbtime = 4
    gdat.rtag = ['pss7pntscmp%d' % (numbtime - t - 1) for t in range(numbtime)]
    gdat.reco = [7 for t in range(numbtime)]
    gdat.evtc = [2 for t in range(numbtime)]
    gdat.strgtime = ['tmin=239155201 tmax=364953603' for t in range(numbtime)]
    gdat.weekinit = [9 for t in range(numbtime)]
    gdat.weekfinl = [218 for t in range(numbtime)]
    gdat.listtimefrac = [0.25, 0.50, 0.75, 1.]
    gdat.photpath = ['p7v6c' for t in range(numbtime)]
    gdat.strgregi = [' ra=INDEF dec=INDEF rad=INDEF ' for n in range(numbtime)]
    gdat.strgener = ['gtbndefn_pnts.fits']

    gdat.evtt = [4, 8, 16, 32]
    gdat.test = True

    make_maps_main(gdat)


def prep_maps(reco, enertype, regitype):
    
    if enertype == 'back':
        numbener = 30
        minmener = 0.1
        maxmener = 100.
        binsener = logspace(log10(minmener), log10(maxmener), numbener)
    else:
        binsener = array([0.1, 0.3, 1., 3., 10., 100.])

    evtt = array([4, 8, 16, 32])
    
    liststrgener = ['ENERGY1', 'ENERGY2', 'ENERGY3', 'ENERGY4', 'ENERGY5']
    liststrgchan = ['CHANNEL1', 'CHANNEL2', 'CHANNEL3', 'CHANNEL4', 'CHANNEL5']
    listdatatype = ['cmp0', 'cmp1', 'cmp2', 'cmp3', 'full']
    
    cnts = zeros((numbener, numbpixl, numbevtt))
    expo = zeros((numbener, numbpixl, numbevtt))
    flux = zeros((numbener, numbpixl, numbevtt))
    
    for k, datatype in enumerate(listdatatype):
    
        for m in indxevtt:

            if reco == 7:
                if m < 2:
                    continue
                elif m == 2:
                    thisevtt = 2
                elif m == 3:
                    thisevtt = 1
            else:
                thisevtt = evtt[m]

            path = os.environ["PCAT_DATA_PATH"] + '/expo_evtt%03d_pss%d_%s.fits' % (thisevtt, reco, enertype)
            expoarry = pf.getdata(path, 1)
            for i in indxener:
                expo[i, :, m] = expoarry[liststrgener[i]]

            path = os.environ["PCAT_DATA_PATH"] + '/cnts_evtt%03d_pss%d_%s.fits' % (thisevtt, reco, enertype)
            cntsarry = pf.getdata(path)
            for i in indxener:
                cnts[i, :, m] = cntsarry[liststrgchan[i]]

        indxexpo = where(expo > 0.) 
        flux[indxexpo] = cnts[indxexpo] / expo[indxexpo] / apix
        flux /= diffener[:, None, None]

        if regitype == 'ngal':
            for i in indxener:
                for m in indxevtt:
                    
                    if reco == 7:
                        if m < 2:
                            continue
                        elif m == 2:
                            thisevtt = 2
                        elif m == 3:
                            thisevtt = 1
                    else:
                        thisevtt = evtt[m]

                    almc = hp.map2alm(flux[i, :, m])
                    hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                    flux[i, :, m] = hp.alm2map(almc, nside)

                    almc = hp.map2alm(expo[i, :, m])
                    hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                    expo[i, :, m] = hp.alm2map(almc, nside)

        path = os.environ["PCAT_DATA_PATH"] + '/fermexpo_pss%d_%s.fits' % (reco, enertype)
        pf.writeto(path, expo, clobber=True)

        path = os.environ["PCAT_DATA_PATH"] + '/fermflux_pss%d_%s.fits' % (reco, enertype)
        pf.writeto(path, flux, clobber=True)


def writ_fdfm():
    
    nside = 256
    numbpixl = 12 * nside**2
    numbener = 3
    numbevtt = 4

    binsener = array([0.3, 1., 3., 10.])
    
    fermfdfmfluxigaltemp = tdpy.util.retr_fdfm(binsener, nside)

    fermfdfmfluxngaltemp = zeros((numbener, numbpixl))
    for i in range(numbener):
        almc = hp.map2alm(fermfdfmfluxigaltemp[i, :])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        fermfdfmfluxngaltemp[i, :] = hp.alm2map(almc, nside)

    fermfdfmfluxngal = zeros((numbener, numbpixl, numbevtt))
    fermfdfmfluxigal = zeros((numbener, numbpixl, numbevtt))
    for m in range(numbevtt):
        fermfdfmfluxigal[:, :, m] = fermfdfmfluxigaltemp
        fermfdfmfluxngal[:, :, m] = fermfdfmfluxngaltemp

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    pf.writeto(path, fermfdfmfluxigal, clobber=True)

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_ngal.fits'
    pf.writeto(path, fermfdfmfluxngal, clobber=True)


def prep_dust():

    minmlgal = -20.
    maxmlgal = 20.
    minmbgal = -20.
    maxmbgal = 20.

    path = os.environ["PCAT_DATA_PATH"] + '/fermfdfmflux_igal.fits'
    fdfmflux = pf.getdata(path)

    path = os.environ["PCAT_DATA_PATH"] + '/lambda_sfd_ebv.fits'
    dustigal = pf.getdata(path)['TEMPERATURE']
    numbside = int(sqrt(dustigal.size / 12))
    dustigal = hp.reorder(dustigal, n2r=True)

    dustigal = dustigal[None, :, None] * mean(fdfmflux, 1)[:, None, :] / mean(dustigal)

    path = os.environ["PCAT_DATA_PATH"] + '/fermdustflux_igal.fits'
    pf.writeto(path, dustigal, clobber=True)

    dustngal = empty_like(dustigal)
    numbener = 3
    for i in range(numbener):
        almc = hp.map2alm(dustigal[i, :, 0])
        hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
        dustngal[i, :, :] = hp.alm2map(almc, numbside)[:, None]

    path = os.environ["PCAT_DATA_PATH"] + '/fermdustflux_igal.fits'
    pf.writeto(path, dustigal, clobber=True)
    
    path = os.environ["PCAT_DATA_PATH"] + '/fermdustflux_ngal.fits'
    pf.writeto(path, dustngal, clobber=True)
   

def plot_spec():

    # Fermi-LAT best-fit components at the NGP
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

    
def cnfg_ferm_info():
    
    minmflux = array([3e-10, 1e-10, 3e-11, 1e-11])
    numbruns = minmflux.size
    maxmnumbpnts = zeros(numbruns, dtype=int) + 1000
    numbswep = zeros(numbruns, dtype=int) + 2000000
    numbburn = numbswep / 2
    
    numbiter = minmflux.size

    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)
    
    strgexpo='fermexpo_cmp0_ngal.fits'
    strgexpr='fermflux_cmp0_ngal.fits'

    indxenerincl = arange(2, 3)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size

    for k in range(numbiter):
        
        gridchan = pcat.main.init( \
                                  psfntype='doubking', \
                                  numbswep=numbswep[k], \
                                  numbburn=numbburn[k], \
                                  probprop=array([0.01, 0.01, 0., 0., 1., 1., 0, 0, 1., 1., 1., 1.], dtype=float), \
                                  trueinfo=True, \
                                  randinit=False, \
                                  makeplot=True, \
                                  maxmgang=10., \
                                  maxmnumbpnts=array([maxmnumbpnts[k]]), \
                                  indxenerincl=indxenerincl, \
                                  indxevttincl=indxevttincl, \
                                  minmflux=minmflux[k], \
                                  maxmflux=1e-7, \
                                  regitype='ngal', \
                                  pathdata=os.environ["FERM_NGAL_DATA_PATH"], \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  maxmnormback=array([5., 5.]), \
                                  minmnormback=array([0.2, 0.2]), \
                                  strgexpo=strgexpo, \
                                  datatype='inpt', \
                                  strgexpr=strgexpr, \
                                 )
        
        listlevi[k] = gridchan[-2]
        listinfo[k] = gridchan[-1]

    plot_minmfluxinfo(minmflux, listinfo, listlevi)


def intr_ferm_expr_ngal( \
                        strgexpr='fermflux_cmp0_ngal.fits', \
                        strgexpo='fermexpo_cmp0_ngal.fits', \
                       ): 

    karg = {}
    karg['psfntype'] = 'doubking'
    karg['numbswep'] = 2000000
    karg['randinit'] = False
    # temp
    karg['boolproppsfn'] = False
    karg['maxmgang'] = 20.
    karg['initfdfnslop'] = array([1.9])
    karg['initfdfnnorm'] = array([300])
    karg['maxmnumbpnts'] = array([500])
    karg['indxenerincl'] = arange(1, 4)
    karg['indxevttincl'] = arange(2, 4)
    karg['minmflux'] = 3e-11
    karg['maxmflux'] = 1e-7
    karg['regitype'] = 'ngal'
    karg['pathdata'] = os.environ["FERM_NGAL_DATA_PATH"]
    karg['maxmnormback'] = array([2., 2.])
    karg['minmnormback'] = array([0.5, 0.5])
    karg['strgback'] = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
    karg['strgexpo'] = strgexpo
    karg['datatype'] = 'inpt'
    karg['strgexpr'] = strgexpr

    return karg


def cnfg_ferm_expr_ngal():
    karg = intr_ferm_expr_ngal()
    pcat.main.init(**karg)


def cnfg_ferm_expr_ngal_cmp1():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_cmp1_ngal.fits', strgexpo='fermexpo_cmp1_ngal.fits')
    pcat.main.init(**karg)


def cnfg_ferm_expr_ngal_cmp2():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_cmp2_ngal.fits', strgexpo='fermexpo_cmp2_ngal.fits')
    pcat.main.init(**karg)


def cnfg_ferm_expr_ngal_cmp3():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_cmp3_ngal.fits', strgexpo='fermexpo_cmp3_ngal.fits')
    pcat.main.init(**karg)


def cnfg_ferm_mock_ngal():
     
    indxenerincl = arange(1, 4)
    indxevttincl = arange(2, 4)
    numbener = indxenerincl.size

    minmflux = 3e-11
    maxmflux = 1e-7
    mockfdfnslop = array([1.9])
    
    pcat.main.init(psfntype='doubking', \
                   numbswep=200000, \
                   randinit=False, \
                   trueinfo=True, \
                   maxmgang=20., \
                   indxenerincl=indxenerincl, \
                   indxevttincl=indxevttincl, \
                   mocknumbpnts=array([300]), \
                   maxmnumbpnts=array([600]), \
                   minmflux=minmflux, \
                   maxmflux=maxmflux, \
                   regitype='ngal', \
                   maxmnormback=array([2., 2.]), \
                   minmnormback=array([0.5, 0.5]), \
                   pathdata=os.environ["FERM_NGAL_DATA_PATH"], \
                   strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                   strgexpo='fermexpo_cmp0_ngal.fits', \
                   datatype='mock', \
                   numbsideheal=256, \
                   mockfdfnslop=mockfdfnslop, \
                   mocknormback=ones((2, numbener)), \
                  )

    
if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

