from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_chandata

def writ_data():

    print 'Producing CDF-S images for PCAT...'
    
    pixlsize = deg2rad(0.492 / 3600.)
    apix = pixlsize**2
            
    numbevtt = 1
  
	# read raw files 
    strgproc = os.uname()[1]
    if strgproc == 'fink1.rc.fas.harvard.edu' or strgproc == 'fink2.rc.fas.harvard.edu':
        path = os.environ["TDGU_DATA_PATH"] + '/xray_back/data'
        strgvarb = ['thresh.expmap', 'flux.img']
        strgvarbmine = ['expo', 'sbrt']
        for a in range(2):
	        for expo in [2, 4, 7]:
	            for i in range(5):
	                cmnd = 'mkdir -p %s/%dmsc' % (path, expo)
	                print cmnd
	                os.system(cmnd)
	                cmnd = 'cp /n/fink1/rfeder/obsids/full/merged_%dMs/merged_%dMs_%d/%dMs_%d_%s %s/%dmsc/%s%04d.fits' \
	                                                                            % (expo, expo, i, expo, i, strgvarb[a], path, expo, strgvarbmine[a], i)
	                print cmnd
	                #os.system(cmnd)

    #listnumbside = [300, 480]
    listnumbside = [300]
    listdatatype = ['home', 'extr']
    for datatype in listdatatype:
        print 'datatype'
        print datatype

        if datatype == 'home':
            binsener = array([0.5, 0.91, 1.66, 3.02, 5.49, 10.])
            expomaps = [2, 4, 7]
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
                            #path = pathdata + '%dmsc/sbrt%04d.fits' % (expomaps, i)
                            temp = pf.getdata(path, 0)

                    numbsideyaxi = temp.shape[0]
                    numbsidexaxi = temp.shape[1]
                    cntrindx = array([numbsideyaxi, numbsidexaxi]) / 2
                    numbsideshft = 0
                    #cntrindx[0] += numbsideshft
                    #cntrindx[1] -= numbsideshft
                    minmindx = cntrindx - numbside / 2
                    maxmindx = cntrindx + numbside / 2
                    
                    print 'minmindx[0]'
                    print minmindx[0]
                    print 'minmindx[1]'
                    print minmindx[1]
                
                print 'expomaps[k]'
                print expomaps[k]
                
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
                        #path = pathdata + '%.2f-%.2f_thresh.img' % (binsener[i], binsener[i+1])
                        path = pathdata + '%dmsc/sbrt%04d.fits' % (expomaps[k], i)
                        cntp[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]

                        # exposure
                        #path = pathdata + '%.2f-%.2f_thresh.expmap' % (binsener[i], binsener[i+1])
                        path = pathdata + '%dmsc/expo%04d.fits' % (expomaps[k], i)
                        expo[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                
                print 'cntp'
                summgene(cntp)
                
                numbsideyaxi = pf.getdata(path, 0).shape[0]
                numbsidexaxi = pf.getdata(path, 0).shape[1]

                sbrt = zeros_like(cntp)
                for i in indxener:
                    indxtemp = where(expo[i, :, :, 0] > 0.)
                    sbrt[i, indxtemp[0], indxtemp[1], 0] = cntp[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
                
                sbrtback = zeros_like(cntp)
                for i in indxener:
                    indxtemp = where(expo[i, :, :, 0] > 0.)
                    sbrtback[i, indxtemp[0], indxtemp[1], 0] = cntpback[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
                
                if False:
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
                              trueminmsbrt=1e-7, \
                              truenumbpnts=array([1]), \
                              truemaxmnumbpnts=array([6]), \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
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
                              truesbrtimps=array([5e-7]), \
                              truenumbpnts=array([1]), \
                              truemaxmnumbpnts=array([6]), \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
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
                              truenumbpnts=array([50, 40]), \
                             )


def pcat_chan_mock_zero():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=100000, \
                              numbswepplot=10000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              condcatl=False, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                              truenumbpnts=array([0]), \
                             )


def pcat_chan_mock_assc():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=10000, \
                              inittype='refr', \
                              makeplotinit=False, \
                              numbswepplot=10000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              condcatl=False, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )

def pcat_chan_mock():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=2000000, \
                              numbswepplot=20000, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              condcatl=False, \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )

def pcat_chan_mock_delt():
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbswep=1, \
                              numbswepplot=3000, \
                              mockonly=True, \
                              #verbtype=2, \
                              truesbrtdistslop=1.1, \
                              diagmode=True, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              trueminmsbrt=5e-9, \
                              truemaxmsbrt=5e-7, \
                              condcatl=False, \
                              truenumbpnts=array([1000]), \
                              truemaxmnumbpnts=array([1000]), \
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
                              numbswep=2000000, \
                              numbswepplot=20000, \
                              optihess=True, \
                              anlytype=anlytype, \
                              recostat=namestat, \
                              savestat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              condcatl=False, \
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
                              numbswep=200000, \
                              numbswepplot=20000, \
                              #makeplot=False, \
                              optihess=True, \
                              anlytype=anlytype, \
                              recostat=namestat, \
                              savestat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              condcatl=False, \
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
                              numbswep=5000000, \
                              factthin=1000, \
                              numbswepplot=20000, \
                              optihess=True, \
                              anlytype=datatype, \
                              recostat=namestat, \
                              #savestat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              condcatl=False, \
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
                              numbswep=5000000, \
                              factthin=1000, \
                              numbswepplot=20000, \
                              optihess=True, \
                              anlytype=anlytype, \
                              recostat=namestat, \
                              #savestat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              condcatl=False, \
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
                              numbswep=5000000, \
                              factthin=1000, \
                              numbswepplot=20000, \
                              optihess=True, \
                              anlytype=anlytype, \
                              recostat=namestat, \
                              #savestat=namestat, \
                              strgexpo='expochan%s.fits' % rtagdata, \
                              exprtype='chan', \
                              condcatl=False, \
                              numbsidecart=numbsidecart, \
                              strgexprsbrt='sbrtchan%s.fits' % rtagdata, \
                             )


globals().get(sys.argv[1])(*sys.argv[2:])
