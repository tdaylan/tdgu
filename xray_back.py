from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_chandata

def writ_maps(datatype='extr'):

    if datatype == 'home':
        binsener = array([0.5, 0.91, 1.66, 3.02, 5.49, 10.])
        expomaps = [4]
    else:
        expomaps = [2, 4]
        binsener = array([0.5, 2., 8.])
    
    print 'Producing CDF-S images for PCAT...'
    print 'datatype'
    print datatype

    pixlsize = deg2rad(0.492 / 3600.)
    apix = pixlsize**2
            
    pathdata = os.environ["TDGU_DATA_PATH"] + '/xray_back/data/'

    diffener = binsener[1:] - binsener[:-1]
    numbener = diffener.size
    indxener = arange(numbener)
    numbevtt = 1
    
    cntrindx = array([3260, 3280]) / 2
    numbpixlshft = 30
    cntrindx[0] += numbpixlshft
    numbmaps = len(expomaps)

    listnumbside = [30, 200, 300]
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

        for k in range(numbmaps):
            
            print 'expomaps[k]'
            print expomaps[k]
            
            strgmaps = '%s%dmsc%04d' % (datatype, expomaps[k], numbside)

            cnts = empty((numbener, numbside, numbside, numbevtt))
            expo = empty((numbener, numbside, numbside, numbevtt))
            cntsback = empty((numbener, numbside, numbside, numbevtt))
            
            if datatype == 'extr':
                # count map
                ## soft band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
                cnts[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                ## hard band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-2Ms-2to8-asca-im-bin1-astwk.fits'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-4Ms-2to8-asca-im-bin1.fits'
                cnts[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

                # exposure
                ## soft band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-2Ms-0p5to2-bin1-astwk.emap'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-4Ms-0p5to2-bin1.emap'
                expo[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                ## hard band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-2Ms-2to8-bin1-astwk.emap'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-4Ms-2to8-bin1.emap'
                expo[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            
                # background
                ## soft band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-4Ms-0p5to2-bin1.back'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-2Ms-0p5to2-asca-bkg-bin1.fits'
                cntsback[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                ## hard band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-4Ms-2to8-bin1.back'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-2Ms-2to8-asca-bkg-bin1.fits'
                cntsback[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            
            if datatype == 'home':
                for i in indxener:
                    # count map
                    path = pathdata + '%.2f-%.2f_thresh.img' % (binsener[i], binsener[i+1])
                    cnts[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

                    # exposure
                    path = pathdata + '%.2f-%.2f_thresh.expmap' % (binsener[i], binsener[i+1])
                    expo[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            
            flux = zeros_like(cnts)
            for i in indxener:
                indxtemp = where(expo[i, :, :, 0] > 0.)
                flux[i, indxtemp[0], indxtemp[1], 0] = cnts[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix

            fluxback = zeros_like(cnts)
            for i in indxener:
                indxtemp = where(expo[i, :, :, 0] > 0.)
                fluxback[i, indxtemp[0], indxtemp[1], 0] = cntsback[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
            
            for i in indxener:
                print 'i'
                print i
                print 'expo[i, :, :, 0]'
                summgene(expo[i, :, :, 0])
                print 'cnts[i, :, :, 0]'
                summgene(cnts[i, :, :, 0])
                print 'flux[i, :, :, 0]'
                summgene(flux[i, :, :, 0])

            pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
            
            path = pathdatapcat + 'chanexpo%s.fits' % strgmaps
            print 'Writing to %s...' % path
            pf.writeto(path, expo, clobber=True)

            path = pathdatapcat + 'chanflux%s.fits' % strgmaps
            print 'Writing to %s...' % path
            pf.writeto(path, flux, clobber=True)
            
            if datatype == 'extr':
                path = pathdatapcat + 'chanfluxback%s.fits' % strgmaps
                print 'Writing to %s...' % path
                pf.writeto(path, fluxback, clobber=True)
            
        print
        print
        print
        print


def pcat_chan_mock_zero():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=300, \
                              diagmode=False, \
                              truenumbpnts=array([0]), \
                             )

def pcat_chan_mock():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              #verbtype=2, \
                              diagmode=True, \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )

# test suites

def pcat_chan_mock_test():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=1000, \
                              factthin=200, \
                              #verbtype=2, \
                              shrtfram=True, \
                              inittype='refr', \
                              trueminmflux=1e-7, \
                              truenumbpnts=array([1]), \
                              truemaxmnumbpnts=array([6]), \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=300, \
                             )


def pcat_chan_mock_spmr():
   
    anglfact = 3600. * 180. / pi
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 30
    
    maxmgangdata = 0.492 / anglfact * numbsidecart / 2.
    
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              numbburn=0, \
                              probbrde=0., \
                              probtran=1., \
                              priofactdoff=0., \
                              checprio=True, \
                              indxenerincl=array([0]), \
                              verbtype=2, \
                              #makeplot=False, \
                              numbswepplot=1000, \
                              #makeplotfram=False, \
                              shrtfram=True, \
                              inittype='refr', \
                              truelgalimps=array([0.]), \
                              truebgalimps=array([0.]), \
                              truefluximps=array([5e-7]), \
                              truenumbpnts=array([1]), \
                              truemaxmnumbpnts=array([6]), \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              maxmgangdata=maxmgangdata, \
                              numbsidecart=numbsidecart, \
                             )


def pcat_chan_mock_popl():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=300, \
                              truenumbpnts=array([50, 40]), \
                             )


def pcat_chan_inpt():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprflux='chanflux_%04d_4msc.fits' % numbsidecart, \
                             )


globals().get(sys.argv[1])(*sys.argv[2:])
