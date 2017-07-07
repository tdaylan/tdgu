from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_chandata

def writ_data(datatype='extr'):

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
    meanener = sqrt(binsener[1:] * binsener[:-1])
    indxener = arange(numbener)
    numbevtt = 1
    

    numbmaps = len(expomaps)

    #listnumbside = [30, 200, 300]
    listnumbside = [300, 1000]
    for numbside in listnumbside:

        print 'numbside'
        print numbside

        for k in range(numbmaps):
            
            print 'expomaps[k]'
            print expomaps[k]
            
            strgmaps = '%s%dmsc%04d' % (datatype, expomaps[k], numbside)

            # determine map shape
            if k == 0:
                if datatype == 'extr':
                    # count map
                    ## soft band
                    if expomaps[k] == 2:
                        path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
                    elif expomaps[k] == 4:
                        path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
                    temp = pf.getdata(path, 0)
                if datatype == 'home':
                    for i in indxener:
                        # count map
                        path = pathdata + '%.2f-%.2f_thresh.img' % (binsener[i], binsener[i+1])
                        #path = pathdata + '%dmsc/flux%04d.fits' % (expomaps, i)
                        temp = pf.getdata(path, 0)

                numbsideyaxi = temp.shape[0]
                numbsidexaxi = temp.shape[1]
                cntrindx = array([numbsideyaxi, numbsidexaxi]) / 2
                numbsideshft = 0
                #cntrindx[0] += numbsideshft
                #cntrindx[1] -= numbsideshft
                minmindx = cntrindx - numbside / 2
                maxmindx = cntrindx + numbside / 2
                
                print 'numbsideyaxi'
                print numbsideyaxi
                print 'numbsidexaxi'
                print numbsidexaxi
                print 'minmindx[0]'
                print minmindx[0]
                print 'minmindx[1]'
                print minmindx[1]
                #print 'maxmindx[0]'
                #print maxmindx[0]
                #print 'maxmindx[1]'
                #print maxmindx[1]
                print
            
            cntp = empty((numbener, numbside, numbside, numbevtt))
            expo = empty((numbener, numbside, numbside, numbevtt))
            cntpback = empty((numbener, numbside, numbside, numbevtt))
            
            if datatype == 'extr':
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
                cntp[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

                ## hard band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-2Ms-2to8-asca-im-bin1-astwk.fits'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-4Ms-2to8-asca-im-bin1.fits'
                cntp[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

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
                cntpback[0, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                ## hard band
                if expomaps[k] == 2:
                    path = pathdata + 'CDFS-4Ms-2to8-bin1.back'
                elif expomaps[k] == 4:
                    path = pathdata + 'CDFS-2Ms-2to8-asca-bkg-bin1.fits'
                cntpback[1, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
            
            if datatype == 'home':
                for i in indxener:
                    # count map
                    path = pathdata + '%.2f-%.2f_thresh.img' % (binsener[i], binsener[i+1])
                    #path = pathdata + '%dmsc/flux%04d.fits' % (expomaps, i)
                    cntp[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

                    # exposure
                    path = pathdata + '%.2f-%.2f_thresh.expmap' % (binsener[i], binsener[i+1])
                    #path = pathdata + '%dmsc/expo%04d.fits' % (expomaps, i)
                    expo[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
             
            numbsideyaxi = pf.getdata(path, 0).shape[0]
            numbsidexaxi = pf.getdata(path, 0).shape[1]

            flux = zeros_like(cntp)
            for i in indxener:
                indxtemp = where(expo[i, :, :, 0] > 0.)
                flux[i, indxtemp[0], indxtemp[1], 0] = cntp[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
            
            fluxback = zeros_like(cntp)
            for i in indxener:
                indxtemp = where(expo[i, :, :, 0] > 0.)
                fluxback[i, indxtemp[0], indxtemp[1], 0] = cntpback[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
            
            if False:
                for i in indxener:
                    print 'i'
                    print i
                    print 'expo[i, :, :, 0]'
                    summgene(expo[i, :, :, 0])
                    print 'cntp[i, :, :, 0]'
                    summgene(cntp[i, :, :, 0])
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
    print 'hey'
    print 2.5e-10 / apix


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


def pcat_chan_mock():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              numbswepplot=3000, \
                              diagmode=True, \
                              strgexpo='chanexpo%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              trueminmflux=1e-8, \
                              truemaxmflux=1e-6, \
                              condcatl=False, \
                              truenumbpnts=array([20]), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )

def pcat_chan_inpt():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=50000, \
                              numbswepplot=5000, \
                              recostat=True, \
                              optihess=True, \
                              savestat=True, \
                              fittmaxmnumbpnts=array([40]), \
                              strgexpo='chanexpo%s.fits' % rtagdata, \
                              exprtype='chan', \
                              condcatl=False, \
                              numbsidecart=numbsidecart, \
                              strgexprflux='chanflux%s.fits' % rtagdata, \
                             )


globals().get(sys.argv[1])(*sys.argv[2:])
