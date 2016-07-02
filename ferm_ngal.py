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


def prep_maps():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    liststrgener = ['ENERGY1', 'ENERGY2', 'ENERGY3', 'ENERGY4', 'ENERGY5']
    liststrgchan = ['CHANNEL1', 'CHANNEL2', 'CHANNEL3', 'CHANNEL4', 'CHANNEL5']
    listdatatype = ['cmp0', 'cmp1', 'cmp2', 'cmp3', 'full']
    listregitype = ['igal', 'ngal']
    
    cnts = zeros((numbener, numbpixl, numbevtt))
    expo = zeros((numbener, numbpixl, numbevtt))
    flux = zeros((numbener, numbpixl, numbevtt))
    
    for k, datatype in enumerate(listdatatype):
    
        for m in indxevtt:

            if datatype != 'full':
                if m < 2:
                    continue
                elif m == 2:
                    thisevtt = 2
                elif m == 3:
                    thisevtt = 1
            else:
                thisevtt = evtt[m]

            path = os.environ["PCAT_DATA_PATH"] + '/expo_evtt%03d_%s.fits' % (thisevtt, datatype)
            expoarry = pf.getdata(path, 1)
            for i in indxener:
                expo[i, :, m] = expoarry[liststrgener[i]]

            path = os.environ["PCAT_DATA_PATH"] + '/cnts_evtt%03d_%s.fits' % (thisevtt, datatype)
            cntsarry = pf.getdata(path)
            for i in indxener:
                cnts[i, :, m] = cntsarry[liststrgchan[i]]

        indxexpo = where(expo > 0.) 
        flux[indxexpo] = cnts[indxexpo] / expo[indxexpo] / apix
        flux /= diffener[:, None, None]

        for regitype in listregitype:
            if regitype == 'ngal':
                for i in indxener:
                    for m in indxevtt:
                        
                        if datatype != 'full':
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

            path = os.environ["PCAT_DATA_PATH"] + '/fermexpo_%s_%s.fits' % (datatype, regitype)
            pf.writeto(path, expo, clobber=True)

            path = os.environ["PCAT_DATA_PATH"] + '/fermflux_%s_%s.fits' % (datatype, regitype)
            pf.writeto(path, flux, clobber=True)


def writ_isot():
    
    global reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix
    reco, evtc, numbevtt, numbevtt, evtt, numbener, minmener, maxmener, binsener, meanener, diffener, indxener, nside, numbpixl, apix = retr_axes()

    # isotropic background
    path = os.environ["PCAT_DATA_PATH"] + '/iso_P8R2_ULTRACLEAN_V6_v06.txt'
    isotdata = loadtxt(path)
    enerisot = isotdata[:, 0] * 1e-3 # [GeV]
    isotflux = isotdata[:, 1] * 1e3 # [1/cm^2/s/sr/GeV]
    isotfluxheal = empty((numbener, numbpixl, numbevtt))
    nsampbins = 10
    enersamp = logspace(log10(amin(binsener)), log10(amax(binsener)), nsampbins * numbener)
    isotflux = interpolate.interp1d(enerisot, isotflux)(enersamp)
    for i in range(numbener):
        isotfluxheal[i, :, :] = trapz(isotflux[i*nsampbins:(i+1)*nsampbins], enersamp[i*nsampbins:(i+1)*nsampbins]) / diffener[i]
        
    path = os.environ["PCAT_DATA_PATH"] + '/fermisotflux.fits'
    pf.writeto(path, isotfluxheal, clobber=True)


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
    
    

