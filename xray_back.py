from __init__ import *

def prep_maps():

    numbside = 300
    cntrindx = array([3260, 3280]) / 2
    minmindx = cntrindx - numbside / 2
    maxmindx = cntrindx + numbside / 2

    binsener = array([0.5, 2., 8.]) * 1e-3
    diffener = binsener[1:] - binsener[:-1]
    numbener = diffener.size
    indxener = arange(numbener)
    numbevtt = 1
    
    expomaps = [2, 4]
    numbmaps = len(expomaps)
    for k in range(numbmaps):
        strgmaps = '%04d_%dmsc' % (numbside, expomaps[k])

        cnts = empty((numbener, numbside, numbside, numbevtt))
        expo = empty((numbener, numbside, numbside, numbevtt))
        cntsback = empty((numbener, numbside, numbside, numbevtt))

        # paths
        pathdata = tdpy.util.retr_path('tdgu', 'xray_back/', onlydata=True)
        path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
        if expomaps[k] == 2:
            path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
        elif expomaps[k] == 4:
            path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
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
        

def pcat_chan_mock():
    
    maxmgang = deg2rad(0.492 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              numbswep=100, \
                              #numbswep=200000, \
                              #numbswepplot=30000, \
                              #numbswepplot=1, \
                              #factthin=1000, \
                              exprinfo=False, \
                              randinit=False, \
                              modlpsfntype='singgaus', \
                              indxenerincl=arange(2), \
                              indxevttincl=arange(1), \
                              #verbtype=2, \
                              strgback=['unit'], \
                              #strgback=['unit', 'chanfluxback_0200_4msc.fits'], \
                              strgexpo='chanexpo_0200_4msc.fits', \
                              exprtype='chan', \
                              #mockvarioaxi=False, \
                              #modlvarioaxi=False, \
                              maxmgang=maxmgang, \
                              scalmaps='asnh', \
                              #minmnormback=array([5e2]), \
                              #maxmnormback=array([1e4]), \
                              minmflux=3e-7, \
                              maxmflux=1e-4, \
                              datatype='mock', \
                              mocknormback=ones((1, 2)), \
                              numbsidecart=200, \
                              #mocknumbpnts=array([200]), \
                              maxmnumbpnts=array([2]), \
                              mocknumbpnts=array([2]), \
                             )

def pcat_chan_inpt():
    
    maxmgang = deg2rad(0.492 / 3600.) * 100.
    gridchan = pcat.main.init( \
                              numbswep=100, \
                              numbswepplot=1, \
                              exprinfo=False, \
                              randinit=True, \
                              indxenerincl=arange(2), \
                              mockvarioaxi=False, \
                              modlvarioaxi=False, \
                              scalmaps='asnh', \
                              #verbtype=2, \
                              indxevttincl=arange(1), \
                              #strgback=['chanfluxback_0300_4msc.fits'], \
                              strgback=['unit'], \
                              #strgback=['unit', 'chanfluxback_0200_4msc.fits'], \
                              #strgexpo='chanexpo_0300_4msc.fits', \
                              strgexpo='chanexpo_0200_4msc.fits', \
                              datatype='inpt', \
                              initfluxdistslop=array([2.]), \
                              #strgexpr='chanflux_0300_4msc.fits', \
                              strgexpr='chanflux_0200_4msc.fits', \
                              exprtype='chan', \
                              maxmnumbpnts=array([1]), \
                              maxmgang=maxmgang, \
                              minmnormback=array([1e4]), \
                              maxmnormback=array([1e5]), \
                              minmflux=1e-7, \
                              maxmflux=1e-4, \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

