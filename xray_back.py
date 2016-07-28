from __init__ import *


def prep_maps():

    minmindx = 1000 # 10 # 1600
    maxmindx = 1020 # 10 # 1600

    # flux
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_softw.fits'
    fluxtemp = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    flux = empty((2, fluxtemp.shape[0], fluxtemp.shape[1], 1))
    expo = empty((2, fluxtemp.shape[0], fluxtemp.shape[1], 1))
    flux[0, :, :, 0] = fluxtemp
    
    path = os.environ["TDGU_DATA_PATH"] + 'img_940k_hardw.fits'
    fluxtemp = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    flux[1, :, :, 0] = fluxtemp

    path = os.environ["TDGU_DATA_PATH"] + 'chanflux.fits'
    pf.writeto(path, flux, clobber=True)

    # exposure
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_soft.fits'
    expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
    path = os.environ["TDGU_DATA_PATH"] + 'expmap_940k_hard.fits'
    expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx:maxmindx, minmindx:maxmindx]
   
    print 'hey'
    print amin(flux), amax(flux), mean(flux), std(flux)
    print 

    print 'hey'
    print amin(expo), amax(expo), mean(expo), std(expo)
    print 

    path = os.environ["TDGU_DATA_PATH"] + 'chanexpo.fits'
    pf.writeto(path, expo, clobber=True)


def pcat_chan():
    
    gridchan = pcat.main.init( \
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
                              maxmnumbpnts=array([100]), \
                              binsenerfull=array([0.5e-3, 2e-3, 7e-3]), \
                              maxmgang=10., \
                              exprtype='chan', \
                              minmflux=1e-10, \
                              maxmflux=1e-7, \
                             )
        

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

