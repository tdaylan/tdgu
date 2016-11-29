from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_chandata


def prep_maps():

    binsener = array([0.5, 2., 8.]) * 1e-3
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

    rasc = (3. + 32. / 60. + 28.06 / 3600.) * 15.
    decl = -27. - 48. / 60. - 26.4 / 3600. + 14.76 / 3600.
    coor = SkyCoord(ra=rasc, dec=decl, unit='deg', equinox='J2000')
    lgal = coor.galactic.l.degree
    bgal = coor.galactic.b.degree
    print 'lgal'
    print lgal
    print 'bgal'
    print bgal

    listnumbside = [200, 300, 1500]
    for numbside in listnumbside:
        minmindx = cntrindx - numbside / 2
        maxmindx = cntrindx + numbside / 2

        declshft = 0.492 * (minmindx[0] + numbside / 2. - 3260. / 2) / 3600.
        rascshft = 0.492 * (minmindx[1] + numbside / 2. - 3280. / 2) / 3600.
        declcntr = decl + declshft 
        rasccntr = rasc + rascshft 
        coor = SkyCoord(ra=rasccntr, dec=declcntr, unit='deg', equinox='J2000')
        lgalcntr = coor.galactic.l.degree
        bgalcntr = coor.galactic.b.degree
        
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
        print 'bgalcntr'
        print bgalcntr
        print 'lgalcntr'
        print lgalcntr
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

        print
        print
            #print 'expo'
            #print mean(expo)
            #print 'diffener'
            #print diffener
            #print 'apix'
            #print apix
            #print 


def pcat_chan_mock():
    
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              factthin=100, \
                              strgback=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              datatype='mock', \
                              lgalcntr=deg2rad(223.562517912), \
                              bgalcntr=deg2rad(-54.4384411082), \
                              numbsidecart=numbsidecart, \
                              maxmnumbpnts=array([100]), \
                              verbtype=2, \
                              #proppsfp=False, \
                              #mocknumbpnts=array([5]), \
                              mockbacp=zeros((1, 2)), \
                             )

def pcat_chan_mock_popl():
    
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=210000, \
                              numbburn=10000, \
                              factthin=200, \
                              numbswepplot=20000, \
                              strgback=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              datatype='mock', \
                              numbsidecart=numbsidecart, \
                              maxmnumbpnts=array([5, 4]), \
                              mocknumbpnts=array([5, 4]), \
                             )


def pcattemp(args, gdat):
    
    anglcatlrttr = args[0]
    lgalcntr = args[1]
    bgalcntr = args[2]
    
    print 'pcattemp'
    print 'anglcatlrttr'
    print anglcatlrttr
    print 'lgalcntr'
    print lgalcntr
    print 'bgalcntr'
    print bgalcntr

    rttr = hp.rotator.Rotator(rot=[rad2deg(lgalcntr), rad2deg(bgalcntr), anglcatlrttr], deg=True, eulertype='ZYX')
    exprbgalrttr, exprlgalrttr = rttr(pi / 2. - gdat.exprbgal, gdat.exprlgal)
    exprbgalrttr = pi / 2. - exprbgalrttr
     
    indx = tdpy.util.corr_catl(gdat.lgalmaxm, gdat.bgalmaxm, exprlgalrttr, exprbgalrttr, anglassc=pi)
    dist = mean(sqrt((gdat.lgalmaxm[indx] - exprlgalrttr)**2 + (gdat.bgalmaxm[indx] - exprbgalrttr)**2)) * gdat.anglfact
    
    print 'gdat.lgalmaxm'
    print gdat.lgalmaxm
    print 'gdat.bgalmaxm'
    print gdat.bgalmaxm
    print 'exprlgalrttr'
    print exprlgalrttr
    print 'exprbgalrttr'
    print exprbgalrttr
    print 

#    dist = pcat.main.init( \
#                  verbtype=0, \
#                  makeplot=False, \
#                  numbswep=1300000, \
#                  numbburn=300000, \
#                  factthin=1000, \
#                  
#                  lgalcntr=deg2rad(lgalcntr), \
#                  bgalcntr=deg2rad(bgalcntr), \
#                  anglcatlrttr=anglcatlrttr, \
#                  
#                  # 4 ms includes shift
#                  #lgalcntr=deg2rad(223.557580277), \
#                  #bgalcntr=deg2rad(-54.4358432488), \
#                  # 4 ms v2
#                  #lgalcntr=deg2rad(223.564551147), \
#                  #bgalcntr=deg2rad(-54.4364535238), \
#                  # 4 ms
#                  #lgalcntr=deg2rad(223.57152222), \
#                  #bgalcntr=deg2rad(-54.4370634), \
#                  # 1 ms
#                  #lgalcntr=deg2rad(223.57318365), \
#                  #bgalcntr=deg2rad(-54.43741081), \
#                  # mean
#                  #lgalcntr=deg2rad(223.562517912), \
#                  #bgalcntr=deg2rad(-54.4384411082), \
#                  strgback=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
#                  strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
#                  exprtype='chan', \
#                  datatype='inpt', \
#                  numbsidecart=numbsidecart, \
#                  strgexpr='chanflux_%04d_4msc.fits' % numbsidecart, \
#                  #maxmnumbpnts=array([10]), \
#                 )
#
    return dist


def pcat_chan_catl():
    
    gdat = tdpy.util.gdatstrt()
    
    numbsidecart = 1500
    gdat.anglfact = 3600. * pi / 180.
    maxmgang = 0.492 * numbsidecart / 2. / gdat.anglfact
    minmlgal = -maxmgang
    maxmlgal = maxmgang
    minmbgal = -maxmgang
    maxmbgal = maxmgang
    binslgalcart = linspace(minmlgal, maxmlgal, numbsidecart + 1)
    binsbgalcart = linspace(minmbgal, maxmbgal, numbsidecart + 1)
    gdat.lgalcart = (binslgalcart[0:-1] + binslgalcart[1:]) / 2.
    gdat.bgalcart = (binsbgalcart[0:-1] + binsbgalcart[1:]) / 2.
    gdat.numbener = 2    
    gdat.numbevtt = 1
    gdat.pathdata = os.environ["PCAT_DATA_PATH"] + '/data/'
    retr_chandata(gdat) 
    
    path = gdat.pathdata + 'inpt/chanflux_%04d_4msc.fits' % numbsidecart
    data = pf.getdata(path)
    gdat.indxxaximaxm, gdat.indxyaximaxm = tdpy.util.retr_indximagmaxm(data)
   
    gdat.lgalmaxm = gdat.lgalcart[gdat.indxxaximaxm]
    gdat.bgalmaxm = gdat.bgalcart[gdat.indxyaximaxm]

    thissampvarb = array([170., deg2rad(223.5703), deg2rad(-54.4265)])
    stdvpara = 1e-3 * ones(thissampvarb.size)
    tdpy.util.minm(thissampvarb, pcattemp, stdvpara=stdvpara, gdat=gdat, verbtype=2)


def pcat_chan_inpt():
   
    numbsidecart = 1500
    gridchan = pcat.main.init( \
                              numbswep=1300000, \
                              numbburn=300000, \
                              factthin=1000, \
                              strgback=['chanfluxisot_%04d_4msc.fits' % numbsidecart], \
                              strgexpo='chanexpo_%04d_4msc.fits' % numbsidecart, \
                              exprtype='chan', \
                              datatype='inpt', \
                              numbsidecart=numbsidecart, \
                              strgexpr='chanflux_%04d_4msc.fits' % numbsidecart, \
                             )


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

