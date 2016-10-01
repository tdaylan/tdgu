from __init__ import *


def prep_maps():

    minmindx = 600
    maxmindx = 1000
    numbsideoutp = 200

    numbside = maxmindx - minmindx

    binsener = array([0.5, 2., 7.]) * 1e-3
    diffener = binsener[1:] - binsener[:-1]
    numbener = diffener.size
    indxener = arange(numbener)

    numbener = 2
    numbevtt = 1
    
    cnts = empty((numbener, numbside, numbside, numbevtt))
    expo = empty((numbener, numbside, numbside, numbevtt))

    # counts
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_softw.fits'
    cnts[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_hardw.fits'
    cnts[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]

    # exposure
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_soft.fits'
    expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_hard.fits'
    expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
  
    
    pixlsize = deg2rad(0.984 / 3600.) * numbside / numbsideoutp
    apix = pixlsize**2
    
    if numbsideoutp != numbside:
        cntstemp = copy(cnts)
        expotemp = copy(expo)
        cnts = empty((numbener, numbsideoutp, numbsideoutp, numbevtt))
        expo = empty((numbener, numbsideoutp, numbsideoutp, numbevtt))
        for i in arange(numbener):
            cnts[i, :, :, 0] = tdpy.util.rebn(cntstemp[i, :, :, 0], (numbsideoutp, numbsideoutp), totl=True)
            expo[i, :, :, 0] = tdpy.util.rebn(expotemp[i, :, :, 0], (numbsideoutp, numbsideoutp), totl=True)

    print 'cnts'
    print cnts[0, :30, :30, 0].astype(int)

    flux = zeros_like(cnts)
    for i in indxener:
        indxtemp = where(expo[i, :, :, 0] > 0.)
        flux[i, indxtemp[0], indxtemp[1], 0] = cnts[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
   
    path = os.environ["TDGU_DATA_PATH"] + 'chanflux.fits'
    pf.writeto(path, flux, clobber=True)

    path = os.environ["TDGU_DATA_PATH"] + 'chanexpo.fits'
    pf.writeto(path, expo, clobber=True)


def chan_inpt():
    
    maxmgang = deg2rad(0.984 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              pathdata=os.environ["TDGU_DATA_PATH"], \
                              #verbtype=2, \
                              numbswep=100000, \
                              numbswepplot=30000, \
                              numbburn=0, \
                              #factthin=100, \
                              exprinfo=False, \
                              randinit=True, \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              datatype='inpt', \
                              strgexpr='chanflux.fits', \
                              #probprop=array([0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0., 1., 1., 1., 1.], dtype=float), \
                              #initfluxdistslop=array([1.]), \
                              #initnumbpnts=array([50]), \
                              maxmnumbpnts=array([200]), \
                              maxmgang=maxmgang, \
                              binsenerfull=array([0.5e-3, 2e-3, 7e-3]), \
                              minmnormback=array([1e0]), \
                              maxmnormback=array([1e12]), \
                              exprtype='chan', \
                              minmflux=1e-15, \
                              maxmflux=1e-13, \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

