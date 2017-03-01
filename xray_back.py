from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_chandata

def prep_maps():

    binsener = array([0.5, 2., 8.])
    diffener = binsener[1:] - binsener[:-1]
    numbener = diffener.size
    indxener = arange(numbener)
    numbevtt = 1
    
    expomaps = [2, 4]
    
    lgalcntr = 0.492
    cntrindx = array([3260, 3280]) / 2
    numbpixlshft = 30
    cntrindx[0] += numbpixlshft
    numbmaps = len(expomaps)

    listnumbside = [200, 300, 1500]
    for numbside in listnumbside:
        minmindx = cntrindx - numbside / 2
        maxmindx = cntrindx + numbside / 2

        print 'numbside'
        print numbside
        print 'minmindx[0]'
        print minmindx[0]
        print 'maxmindx[0]'
        print maxmindx[0]
        print 'minmindx[1]'
        print minmindx[1]
        print 'maxmindx[1]'
        print maxmindx[1]
        print 
        for k in range(numbmaps):
            strgmaps = '%04d_%dmsc' % (numbside, expomaps[k])

            cnts = empty((numbener, numbside, numbside, numbevtt))
            expo = empty((numbener, numbside, numbside, numbevtt))
            cntsback = empty((numbener, numbside, numbside, numbevtt))

            # paths
            pathdata = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
            path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
            if expomaps[k] == 2:
                path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
            elif expomaps[k] == 4:
                path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
            temp = pf.getdata(path, 0)
            #print 'Shape of the map'
            #print temp.shape

            cnts[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            if expomaps[k] == 2:
                path = pathdata + 'CDFS-2Ms-2to8-asca-im-bin1-astwk.fits'
            elif expomaps[k] == 4:
                path = pathdata + 'CDFS-4Ms-2to8-asca-im-bin1.fits'
            cnts[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

            # exposure
            if expomaps[k] == 2:
                path = pathdata + 'CDFS-2Ms-0p5to2-bin1-astwk.emap'
            elif expomaps[k] == 4:
                path = pathdata + 'CDFS-4Ms-0p5to2-bin1.emap'
            expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            if expomaps[k] == 2:
                path = pathdata + 'CDFS-2Ms-2to8-bin1-astwk.emap'
            elif expomaps[k] == 4:
                path = pathdata + 'CDFS-4Ms-2to8-bin1.emap'
            expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

            # background
            if expomaps[k] == 2:
                path = pathdata + 'CDFS-4Ms-0p5to2-bin1.back'
            elif expomaps[k] == 4:
                path = pathdata + 'CDFS-2Ms-0p5to2-asca-bkg-bin1.fits'
            cntsback[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            if expomaps[k] == 2:
                path = pathdata + 'CDFS-4Ms-2to8-bin1.back'
            elif expomaps[k] == 4:
                path = pathdata + 'CDFS-2Ms-2to8-asca-bkg-bin1.fits'
            cntsback[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            
            pixlsize = deg2rad(0.492 / 3600.)
            apix = pixlsize**2
            
            flux = zeros_like(cnts)
            for i in indxener:
                indxtemp = where(expo[i, :, :, 0] > 0.)
                flux[i, indxtemp[0], indxtemp[1], 0] = cnts[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
            
            fluxback = zeros_like(cnts)
            for i in indxener:
                indxtemp = where(expo[i, :, :, 0] > 0.)
                fluxback[i, indxtemp[0], indxtemp[1], 0] = cntsback[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
            
            pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
            
            path = pathdatapcat + 'chanexpo_%s.fits' % strgmaps
            print 'Writing to %s...' % path
            pf.writeto(path, expo, clobber=True)

            path = pathdatapcat + 'chanflux_%s.fits' % strgmaps
            print 'Writing to %s...' % path
            pf.writeto(path, flux, clobber=True)
            
            path = pathdatapcat + 'chanfluxback_%s.fits' % strgmaps
            print 'Writing to %s...' % path
            pf.writeto(path, fluxback, clobber=True)
            
            path = pathdatapcat + 'chanfluxisot_%s.fits' % strgmaps
            print 'Writing to %s...' % path
            fluxisot = zeros_like(flux)
            fluxisot[0, :, :, 0] = mean(flux[0, :, :, 0])
            fluxisot[1, :, :, 0] = mean(flux[1, :, :, 0])
            pf.writeto(path, fluxisot, clobber=True)


def pcat_chan_mock_zero():
    
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              back=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              numbsidecart=300, \
                              diagmode=False, \
                              truenumbpnts=array([0]), \
                             )

def pcat_chan_mock():
    
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              back=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              numbsidecart=300, \
                              truenumbpnts=array([50]), \
                             )

def pcat_chan_mock_popl():
    
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              back=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              numbsidecart=300, \
                              truenumbpnts=array([50, 40]), \
                             )


def pcat_chan_inpt():
   
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              inittype='refr', \
                              back=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              numbsidecart=300, \
                              strgexprflux='chanflux_%04d_4msc.fits' % numbsidecart, \
                             )


globals().get(sys.argv[1])()
