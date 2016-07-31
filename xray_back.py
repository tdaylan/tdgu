from __init__ import *


def prep_maps():

    minmindx = 0
    maxmindx = 1600

    binsener = array([0.5, 2., 7.])
    diffener = binsener[1:] - binsener[:-1]

    pixlsize = deg2rad(0.984 / 3600.)
    apix = pixlsize**2
    
    # flux
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_softw.fits'
    fluxtemp = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    flux = empty((2, fluxtemp.shape[0], fluxtemp.shape[1], 1))
    expo = empty((2, fluxtemp.shape[0], fluxtemp.shape[1], 1))
    flux[0, :, :, 0] = fluxtemp
    
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_hardw.fits'
    fluxtemp = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    flux[1, :, :, 0] = fluxtemp

    # exposure
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_soft.fits'
    expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_hard.fits'
    expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
   
    flux /= expo * diffener[:, None, None, None] * apix
    flux *= 1e3

    numbener = flux.shape[0]

    fluxtemp = copy(flux)
    expotemp = copy(expo)
    flux = empty((2, 100, 100, 1))
    expo = empty((2, 100, 100, 1))
    for i in arange(numbener):
        flux[i, :, :, 0] = tdpy.util.rebn(fluxtemp[i, :, :, 0], (100, 100))
        expo[i, :, :, 0] = tdpy.util.rebn(expotemp[i, :, :, 0], (100, 100), totl=True)

    path = os.environ["TDGU_DATA_PATH"] + 'chanflux.fits'
    pf.writeto(path, flux, clobber=True)

    path = os.environ["TDGU_DATA_PATH"] + 'chanexpo.fits'
    pf.writeto(path, expo, clobber=True)


def pcat_chan():
    
    maxmgang = deg2rad(0.984 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              verbtype=2, \
                              numbswep=1000, \
                              factthin=1, \
                              numbburn=0, \
                              randinit=True, \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              pathdata=os.environ["TDGU_DATA_PATH"], \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              datatype='inpt', \
                              strgexpr='chanflux.fits', \
                              initnumbpnts=array([3]), \
                              initfluxdistslop=array([2]), \
                              maxmnumbpnts=array([3]), \
                              maxmgang=maxmgang, \
                              binsenerfull=array([0.5e-3, 2e-3, 7e-3]), \
                              exprtype='chan', \
                              minmflux=1e-10, \
                              maxmflux=1e-7, \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

