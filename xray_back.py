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


def pcat_chan_mock():
    
    maxmgang = deg2rad(0.492 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              #numbswep=100, \
                              #numbswepplot=10, \
                              #verbtype=2, \
                              exprinfo=False, \
                              randinit=False, \
                              modlpsfntype='singgaus', \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              exprtype='chan', \
                              #probprop=array([0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0., 1., 1., 1., 1.], dtype=float), \
                              maxmnumbpnts=array([400]), \
                              #maxmnumbpnts=array([4]), \
                              maxmgang=maxmgang, \
                              minmnormback=array([5e2]), \
                              maxmnormback=array([1e4]), \
                              minmflux=1e-8, \
                              maxmflux=1e-5, \
                              datatype='mock', \
                              #mockpsfntype='doubgaus', \
                              mocknormback=ones((1, 2)), \
                              numbsidecart=200, \
                              mocknumbpnts=array([200]), \
                              #mocknumbpnts=array([2]), \
                             )

def pcat_chan_inpt():
    
    maxmgang = deg2rad(0.492 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              numbswep=200000, \
                              #numbswep=100, \
                              #numbswepplot=10, \
                              #verbtype=2, \
                              #numbburn=0, \
                              #factthin=20, \
                              exprinfo=False, \
                              randinit=True, \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              datatype='inpt', \
                              strgexpr='chanflux.fits', \
                              exprtype='chan', \
                              boolproppsfn=True, \
                              #probprop=array([0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.], dtype=float), \
                              #initfluxdistslop=array([1.]), \
                              maxmnumbpnts=array([400]), \
                              #maxmnumbpnts=array([4]), \
                              maxmgang=maxmgang, \
                              minmnormback=array([5e2]), \
                              maxmnormback=array([1e4]), \
                              minmflux=1e-8, \
                              maxmflux=1e-5, \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

