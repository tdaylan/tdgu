from __init__ import *

def prep_maps():

    minmindx = 700
    maxmindx = 900

    numbside = maxmindx - minmindx

    binsener = array([0.5, 2., 7.]) * 1e-3
    diffener = binsener[1:] - binsener[:-1]
    numbener = diffener.size
    indxener = arange(numbener)
    numbevtt = 1
    
    cnts = empty((numbener, numbside, numbside, numbevtt))
    expo = empty((numbener, numbside, numbside, numbevtt))

    # paths
    pathdata = tdpy.util.retr_path('tdgu', 'xray_back/', onlydata=True)
    path = pathdata + 'img_940k_softw.fits'
    cnts[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = pathdata + 'img_940k_hardw.fits'
    cnts[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]

    # exposure
    path = pathdata + 'expmap_940k_soft.fits'
    expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = pathdata + 'expmap_940k_hard.fits'
    expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
  
    
    pixlsize = deg2rad(0.492 / 3600.)
    apix = pixlsize**2
    
    print 'cnts'
    print cnts[0, :30, :30, 0]
    print cnts[0, :30, :30, 0].astype(int)

    flux = zeros_like(cnts)
    for i in indxener:
        indxtemp = where(expo[i, :, :, 0] > 0.)
        flux[i, indxtemp[0], indxtemp[1], 0] = cnts[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
    
    print 'apix'
    print apix

    pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    path = pathdatapcat + 'chanflux.fits'
    print 'Writing to %s' % path
    pf.writeto(path, flux, clobber=True)
    path = pathdatapcat + 'chanexpo.fits'
    pf.writeto(path, expo, clobber=True)


def pcat_chan_mock():
    
    maxmgang = deg2rad(0.492 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              exprinfo=False, \
                              randinit=False, \
                              modlpsfntype='singgaus', \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              modlvarioaxi=False, \
                              exprtype='chan', \
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
                             )

def pcat_chan_inpt():
    
    maxmgang = deg2rad(0.492 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              exprinfo=False, \
                              randinit=True, \
                              indxenerincl=arange(2), \
                              modlvarioaxi=False, \
                              scalmaps='asnh', \
                              indxevttincl=arange(1), \
                              strgback=['unit'], \
                              strgexpo='chanexpo.fits', \
                              datatype='inpt', \
                              initfluxdistslop=array([2.]), \
                              strgexpr='chanflux.fits', \
                              exprtype='chan', \
                              boolproppsfn=True, \
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

