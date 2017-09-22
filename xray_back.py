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

def pcat_chan_mock_spmr(nameconfexec=None):
   
    anglfact = 3600. * 180. / pi
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 15
    
    maxmgangdata = 0.492 / anglfact * numbsidecart / 2.

    dictargs = {}
    dictargs['numbswep'] = 100000
    dictargs['numbburn'] = 0
    dictargs['factthin'] = 100
    dictargs['exprtype'] = 'chan'
    dictargs['inittype'] = 'refr'
    dictargs['probbrde'] = 0.
    dictargs['probtran'] = 1.
    dictargs['indxenerincl'] = array([0])
    dictargs['numbswepplot'] = 1000
    dictargs['truelgalreg0pop00000'] = 0.
    dictargs['truebgalreg0pop00000'] = 0.
    dictargs['truesbrt'] = array([5e-7])
    dictargs['numbelemreg0pop0'] = 1
    dictargs['minmnumbelemreg0pop0'] = 1
    dictargs['maxmnumbelemreg0pop0'] = 2
    dictargs['strgexpo'] = 1e9
    dictargs['maxmgangdata'] = maxmgangdata
    dictargs['numbsidecart'] = numbsidecart
    dictargs['makeplotinit'] = False
    dictargs['shrtfram'] = True
    
    listnameconf = ['pars', 'genebrgt', 'genefain']
    dictargsvari = {}
    for nameconf in listnameconf:
        dictargsvari[nameconf] = {}
    dictargsvari['pars']['priofactdoff'] = 1.
    dictargsvari['genebrgt']['truefluxreg0pop00000'] = 1e-7
    dictargsvari['genefain']['truefluxreg0pop00000'] = 3e-8
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  nameconfexec=nameconfexec, \
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
                              numbswep=10000, \
                              numbburn=0, \
                              factthin=10, \
                              numbswepplot=3000, \
                              #makeplot=False, \
                              #strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              maxmgangdata=1200./3600./180.*pi, \
                              #diagmode=True, \
                              spectype=['gaus'], \
                              strgexpo=1e2, \
                              elemtype='line', \
                              asscrefr=False, \
                              inittype='refr', \
                              #verbtype=2, \
                              exprtype='chan', \
                              anlytype='spec', \
                              numbsidecart=1, \
                              #makeplotinit=False, \
                              #makeplotfram=False, \
                              #propcomp=False, \
                              #prophypr=False, \
                              #propbacp=False, \
                              #probtran=1., \
                              #propcomp=False, \
                              #probtran=0., \
                              #propbacp=False, \
                              maxmnumbelemreg0pop0=10, \
                              #maxmnumbelemreg0pop0=0, \
                              numbelemreg0pop0=5, \
                              #numbelemreg0pop0=0, \
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
def pcat_chan_mock_syst(nameconfexec=None):
   
    dictargs = {}
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 10
    dictargs['numbelemreg0pop0'] = 300
    dictargs['trueminmflux'] = 3e-10
    dictargs['makeplotinit'] = False
    dictargs['shrtfram'] = True
    
    listnameconf = ['nomi', 'unrsbadd', 'unrsfine']
    dictargsvari = {}
    for nameconf in listnameconf:
        dictargsvari[nameconf] = {}
    dictargsvari['nomi']['fittminmflux'] = 3e-10
    dictargsvari['unrsbadd']['fittminmflux'] = 3e-9
    dictargsvari['unrsfine']['fittminmflux'] = 3e-8
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  nameconfexec=nameconfexec, \
                                 )
    

def pcat_chan_mock():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbburn=20000, \
                              factthin=800, \
                              inittype='refr', \
                              numbelemreg0pop0=100, \
                              numbswepplot=10000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )


def pcat_chan_mock():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbburn=10000, \
                              factthin=900, \
                              inittype='refr', \
                              numbelemreg0pop0=100, \
                              numbswepplot=10000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )


def pcat_chan_mock_maxm():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100, \
                              numbburn=0, \
                              evoltype='maxmllik', \
                              factthin=1, \
                              #verbtype=2, \
                              #makeplot=False, \
                              #makeplotinit=False, \
                              #shrtfram=True, \
                              #makeplotfram=False, \
                              inittype='refr', \
                              #killexpo=True, \
                              verbtype=2, \
                              makeplot=False, \
                              numbelemreg0pop0=2, \
                              maxmnumbelemreg0pop0=3, \
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
                              makeplot=False, \
                              #inittype='reco', \
                              #namerecostat=namestat, \
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
                              numbswep=100000, \
                              numbburn=0, \
                              numbswepplot=10000, \
                              anlytype=anlytype, \
                              shrtfram=True, \
                              makeplotinit=False, \
                              inittype='reco', \
                              namerecostat=namestat, \
                              maxmnumbelemreg0pop0=10, \
                              namesavestat=namestat, \
                              savestat=True, \
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
                              #namesavestat=namestat, \
                              #savestat=True, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


globals().get(sys.argv[1])(*sys.argv[2:])
