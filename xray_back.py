from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_refrchaninit

def writ_chan():

    print 'Writing CDF-S dataset for PCAT...'
    
    pixlsize = deg2rad(0.492 / 3600.)
    apix = pixlsize**2
            
    numbevtt = 1
  
    listnumbside = [300, 600]
    listdatatype = ['home', 'extr']
    for datatype in listdatatype:
        print 'datatype'
        print datatype

        if datatype == 'home':
            binsener = array([0.5, 0.91, 1.66, 3.02, 5.49, 10.])
            expomaps = [2, 4, 7]
            strgener = ['0.5', '0.91028', '1.65723', '3.01709', '5.4928', '10.0']
        else:
            expomaps = [2, 4]
            binsener = array([0.5, 2., 8.])
        
        pathdata = os.environ["TDGU_DATA_PATH"] + '/xray_back/data/'

        diffener = binsener[1:] - binsener[:-1]
        numbener = diffener.size
        meanener = sqrt(binsener[1:] * binsener[:-1])
        indxener = arange(numbener)

        numbmaps = len(expomaps)

        for numbside in listnumbside:

            print 'numbside'
            print numbside

            for k in range(numbmaps):
                
                strgmaps = '%s%dmsc%04d' % (datatype, expomaps[k], numbside)

                print 'expomaps[k]'
                print expomaps[k]
                
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
                            path = '/n/fink1/rfeder/obsids/full/merged_%dMs/rest_fov/%d/%s-%s_flux.img' % (expomaps[k], i, strgener[i], strgener[i+1])
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
                
                cntp = zeros((numbener, numbside, numbside, numbevtt))
                expo = zeros((numbener, numbside, numbside, numbevtt))
                cntpback = empty((numbener, numbside, numbside, numbevtt))
                sbrt = zeros_like(cntp)
                sbrtback = zeros_like(cntp)
                
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
                
                    for i in indxener:
                        indxtemp = where(expo[i, :, :, 0] > 0.)
                        sbrtback[i, indxtemp[0], indxtemp[1], 0] = cntpback[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
                
                if datatype == 'home':
                    strgvarb = ['thresh.expmap', 'thresh.img']
                    strgvarbmine = ['expo', 'sbrt']
                    for i in indxener:
                        for a in range(2):
                            path = '/n/fink1/rfeder/obsids/full/merged_%dMs/rest_fov/%d/%s-%s_%s' % (expomaps[k], i, strgener[i], strgener[i+1], strgvarb[a])
                            if a == 0:
                                expo[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                            if a == 1:
                                cntp[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                numbsideyaxi = pf.getdata(path, 0).shape[0]
                numbsidexaxi = pf.getdata(path, 0).shape[1]

                for i in indxener:
                    indxtemp = where(expo[i, :, :, 0] > 0.)
                    sbrt[i, indxtemp[0], indxtemp[1], 0] = cntp[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
                
                if True:
                    for i in indxener:
                        print 'i'
                        print i
                        print 'expo[i, :, :, 0]'
                        summgene(expo[i, :, :, 0])
                        print 'cntp[i, :, :, 0]'
                        summgene(cntp[i, :, :, 0])
                        print 'sbrt[i, :, :, 0]'
                        summgene(sbrt[i, :, :, 0])

                pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
                
                path = pathdatapcat + 'expochan%s.fits' % strgmaps
                print 'Writing to %s...' % path
                pf.writeto(path, expo, clobber=True)

                path = pathdatapcat + 'sbrtchan%s.fits' % strgmaps
                print 'Writing to %s...' % path
                pf.writeto(path, sbrt, clobber=True)
                
                if datatype == 'extr':
                    path = pathdatapcat + 'sbrtchanback%s.fits' % strgmaps
                    print 'Writing to %s...' % path
                    pf.writeto(path, sbrtback, clobber=True)
                
                print
                print
                print
                print
                print


# test suites

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
                              truelgal=array([0.]), \
                              truebgal=array([0.]), \
                              truesbrt=array([5e-7]), \
                              numbelemreg0pop0=1, \
                              minmnumbelemreg0pop0=1, \
                              maxmnumbelemreg0pop0=6, \
                              strgexpo=1e9, \
                              #'expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
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
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=300, \
                              numbpopl=2, \
                              numbelemreg0pop0=50, \
                              numbelemreg0pop1=50, \
                             )


def pcat_chan_mock_spec():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    gridchan = pcat.main.init( \
                              numbswep=1000, \
                              numbswepplot=500, \
                              #makeplot=False, \
                              #strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              #verbtype=2, \
                              maxmgangdata=1200./3600./180.*pi, \
                              spectype=['gaus'], \
                              strgexpo=1e10, \
                              elemtype='line', \
                              asscrefr=False, \
                              #verbtype=2, \
                              exprtype='chan', \
                              anlytype='spec', \
                              numbsidecart=1, \
                              maxmnumbelemreg0pop0=100, \
                              #maxmnumbelemreg0pop0=0, \
                              numbelemreg0pop0=30, \
                             )


def pcat_chan_mock_zero():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbswepplot=10000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              numbelemreg0pop0=0, \
                             )


def pcat_chan_inpt_assc():
    
    datatype = 'home'
    strgexpomaps = '7msc'
    numbsidecart = 300
    namestat = 'pcat_chan_inpt_' + datatype + '%04d' % numbsidecart
    anlytype = datatype + strgexpomaps
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=1000, \
                              factthin=100, \
                              numbswepplot=20000, \
                              anlytype=anlytype, \
                              namerecostat=namestat, \
                              recostat='reco', \
                              #savestat=True, \
                              #namesavestat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


# science suites

def pcat_chan_mock():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbburn=0, \
                              factthin=100, \
                              #verbtype=2, \
                              #makeplot=False, \
                              #makeplotinit=False, \
                              makeplotfram=False, \
                              inittype='refr', \
                              #killexpo=True, \
                              numbelemreg0pop0=100, \
                              maxmnumbelemreg0pop0=200, \
                              numbswepplot=20000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )


def pcat_chan_inpt_extr2msc():
    
    datatype = 'extr'
    strgexpomaps = '2msc'
    numbsidecart = 300
    namestat = 'pcat_chan_inpt_' + datatype + '%04d' % numbsidecart
    anlytype = datatype + strgexpomaps
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbswepplot=10000, \
                              anlytype=anlytype, \
                              inittype='reco', \
                              namerecostat=namestat, \
                              namesavestat=namestat, \
                              savestat=True, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


def pcat_chan_inpt_extr4msc():
    
    datatype = 'extr'
    strgexpomaps = '4msc'
    numbsidecart = 300
    namestat = 'pcat_chan_inpt_' + datatype + '%04d' % numbsidecart
    anlytype = datatype + strgexpomaps
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbswepplot=10000, \
                              anlytype=anlytype, \
                              inittype='reco', \
                              namerecostat=namestat, \
                              namesavestat=namestat, \
                              savestat=True, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


def pcat_chan_inpt_home2msc():
    
    datatype = 'home'
    strgexpomaps = '2msc'
    numbsidecart = 300
    namestat = 'pcat_chan_inpt_' + datatype + '%04d' % numbsidecart
    anlytype = datatype + strgexpomaps
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=2000000, \
                              numbswepplot=10000, \
                              anlytype=datatype, \
                              #inittype='reco', \
                              #namerecostat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


def pcat_chan_inpt_home4msc():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    namestat = 'pcat_chan_inpt_' + datatype + '%04d' % numbsidecart
    anlytype = datatype + strgexpomaps
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=2000000, \
                              numbswepplot=10000, \
                              anlytype=anlytype, \
                              #inittype='reco', \
                              #namerecostat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


def pcat_chan_inpt_home7msc():
    
    datatype = 'home'
    strgexpomaps = '7msc'
    numbsidecart = 300
    namestat = 'pcat_chan_inpt_' + datatype + '%04d' % numbsidecart
    anlytype = datatype + strgexpomaps
    rtagdata = '%s%s%04d' % (datatype, strgexpomaps, numbsidecart)
    gridchan = pcat.main.init( \
                              numbswep=2000000, \
                              numbswepplot=10000, \
                              diagmode=True, \
                              anlytype=anlytype, \
                              #inittype='reco', \
                              #namerecostat=namestat, \
                              namesavestat=namestat, \
                              savestat=True, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


globals().get(sys.argv[1])(*sys.argv[2:])
