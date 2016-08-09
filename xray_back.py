from __init__ import *


def prep_maps():

    minmindx = 0
    maxmindx = 1600

    binsener = array([0.5, 2., 7.]) * 1e-3
    diffener = binsener[1:] - binsener[:-1]
    numbener = diffener.size
    indxener = arange(numbener)

    pixlsize = deg2rad(0.984 / 3600.)
    apix = pixlsize**2
    
    # counts
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_softw.fits'
    cntstemp = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    cnts = empty((2, cntstemp.shape[0], cntstemp.shape[1], 1))
    expo = empty((2, cntstemp.shape[0], cntstemp.shape[1], 1))
    cnts[0, :, :, 0] = cntstemp
    
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_hardw.fits'
    cntstemp = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    cnts[1, :, :, 0] = cntstemp

    # exposure
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_soft.fits'
    expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_hard.fits'
    expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
   
    numbpixloutp = 100
    cntstemp = copy(cnts)
    expotemp = copy(expo)
    cnts = empty((2, numbpixloutp, numbpixloutp, 1))
    expo = empty((2, numbpixloutp, numbpixloutp, 1))
    for i in arange(numbener):
        cnts[i, :, :, 0] = tdpy.util.rebn(cntstemp[i, :, :, 0], (numbpixloutp, numbpixloutp), totl=True)
        expo[i, :, :, 0] = tdpy.util.rebn(expotemp[i, :, :, 0], (numbpixloutp, numbpixloutp), totl=True)

    flux = zeros_like(cnts)
    for i in indxener:
        indxtemp = where(expo[i, :, :, 0] > 0.)
        flux[i, indxtemp[0], indxtemp[1], 0] = cnts[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
   
    print 'diffener'
    print diffener
    print 'apix'
    print apix
    print

    path = os.environ["TDGU_DATA_PATH"] + 'chanflux.fits'
    pf.writeto(path, flux, clobber=True)

    path = os.environ["TDGU_DATA_PATH"] + 'chanexpo.fits'
    pf.writeto(path, expo, clobber=True)


def pcat_chan():
    
    maxmgang = deg2rad(0.984 / 3600.) * 50.
    gridchan = pcat.main.init( \
                              #verbtype=2, \
                              numbswep=1000, \
                              numbswepplot=400, \
                              #numbburn=0, \
                              #factthin=100, \
                              randinit=True, \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              pathdata=os.environ["TDGU_DATA_PATH"], \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              datatype='inpt', \
                              strgexpr='chanflux.fits', \
                              initnumbpnts=array([2]), \
                              initfluxdistslop=array([2]), \
                              maxmnumbpnts=array([2]), \
                              maxmgang=maxmgang, \
                              binsenerfull=array([0.5e-3, 2e-3, 7e-3]), \
                              exprtype='chan', \
                              minmflux=1e-7, \
                              maxmflux=1e-3, \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

