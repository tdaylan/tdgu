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
                            #path = '/n/fink1/rfeder/xray_pcat/obsids/full/merged_%dMs/rest_fov/%d/%s-%s_flux.img' % (expomaps[k], i, strgener[i], strgener[i+1])
                            path = '/n/fink1/rfeder/xray_pcat/merged_flux_10_02/%dMs/%s-%s_flux.img' % (expomaps[k], strgener[i], strgener[i+1])
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
                            if a == 0:
                                path = '/n/fink1/rfeder/xray_pcat/merged_flux_10_02/%dMs/%s-%s_flux.img' % (expomaps[k], strgener[i], strgener[i+1])
                                expo[i, :, :, 0] = pf.getdata(path, 0)[minmindx[0]:maxmindx[0], minmindx[1]:maxmindx[1]]
                            if a == 1:
                                path = '/n/fink1/rfeder/xray_pcat/merged_flux_10_02/%dMs/%s-%s_exposure.img' % (expomaps[k], strgener[i], strgener[i+1])
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

def pcat_chan_mock_spmr(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 15
    
    maxmgangdata = 0.492 / anglfact * numbsidecart / 2.

    dictargs = {}
    dictargs['numbswep'] = 100000
    dictargs['exprtype'] = 'chan'
    dictargs['inittype'] = 'refr'
    dictargs['truelgalpop0reg00000'] = 0.
    dictargs['truebgalpop0reg00000'] = 0.
    dictargs['truefluxpop0reg00000'] = 1e-7
    dictargs['truesbrt'] = array([5e-7])
    dictargs['numbelempop0reg0'] = 1
    dictargs['minmnumbelempop0reg0'] = 1
    dictargs['priofactdoff'] = 0.
    dictargs['strgexpo'] = 1e9
    dictargs['maxmgangdata'] = maxmgangdata
    dictargs['numbsidecart'] = numbsidecart
    dictargs['probtran'] = 1.
    dictargs['probspmr'] = 1.
    dictargs['indxenerincl'] = array([0])
    
    listnamecnfgextn = ['free', 'nomi', 'pars', 'genebrgt', 'genefain', 'psfn']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    #dictargsvari['free']['minmnumbelempop0reg0'] = 0.
    #dictargsvari['free']['maxmnumbelempop0reg0'] = 0.
    dictargsvari['free']['inittype'] = 'rand'
    dictargsvari['free']['probtran'] = 0.4
    dictargsvari['free']['probspmr'] = 0.3
    
    dictargsvari['pars']['priofactdoff'] = 1.
    
    dictargsvari['brgt']['truefluxpop0reg00000'] = 3e-7
    
    dictargsvari['fain']['truefluxpop0reg00000'] = 3e-8
    
    dictargsvari['psfn']['truefluxpop0reg00000'] = 3e-8
    dictargsvari['psfn']['probtran'] = 0.
    dictargsvari['psfn']['propbacp'] = False
    dictargsvari['psfn']['propcomp'] = False
    dictargsvari['psfn']['elemspatevaltype'] = ['full']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_chan_mock_popl(strgcnfgextnexec=None):
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    
    gridchan = pcat.main.init( \
                              exprtype='chan', \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              numbsidecart=300, \
                              numbpopl=2, \
                              numbelempop0reg0=50, \
                              numbelempop1reg0=50, \
                             )


def pcat_chan_mock_spec(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['spatdisttype'] = ['line']
    dictargs['spectype'] = ['gaus']
    dictargs['strgexpo'] = 1e8
    dictargs['elemtype'] = ['lghtline']
    dictargs['inittype'] = 'refr'
    #dictargs['anlytype'] = 'spec'
    # assume a pixel with side 100 arcsec
    dictargs['maxmgangdata'] = 100. / anglfact
    dictargs['numbsidecart'] = 1
    
    # temp
    dictargs['numbelempop0reg0'] = 1
    dictargs['maxmnumbelempop0reg0'] = 2
    dictargs['numbsamp'] = 2000
    dictargs['sqzeprop'] = True
    
    # true < thrs < modl -- trad 
    # true < modl < thrs -- pcat
    
    # thrs < true < modl -- trad
    # modl < true < thrs -- pcat
    
    # thrs < modl < true -- trad
    # modl < thrs < true -- pcat
    listnamecnfgextn = ['nomi', 'truevlow']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['truevlow']['fittminmflux'] = 3e-7
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

# science suites
def pcat_chan_mock(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['numbelempop0reg0'] = 100
    # temp
    #dictargs['strgexpo'] = 'expochanhome4msc0300.fits'
    dictargs['strgexpo'] = 1e9
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 100
    dictargs['priofactdoff'] = 0.
    
    listnamecnfgextn = ['nomi', 'truevlow', 'trueloww', 'truehigh', 'truenone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['truevlow']['trueminmflux'] = 3e-10
    dictargsvari['trueloww']['trueminmflux'] = 1e-9
    dictargsvari['truehigh']['trueminmflux'] = 1e-8
    dictargsvari['truenone']['truenumbelempop0reg0'] = 0
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_chan_mock_maxmllik(strgcnfgextnexec=None):
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              numbburn=0, \
                              evoltype='maxmllik', \
                              factthin=1, \
                              #makeplot=False, \
                              #makeplotinit=False, \
                              #makeplotfram=False, \
                              inittype='refr', \
                              #killexpo=True, \
                              makeplot=False, \
                              numbelempop0reg0=2, \
                              maxmnumbelempop0reg0=3, \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                              numbsidecart=numbsidecart, \
                             )


def pcat_chan_inpt(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['namerecostat'] = 'extr7msc0600'
    
    # temp
    dictargs['diagmode'] = False
    #dictargs['inittype'] = 'reco'
    #dictargs['anlytypedata'] = maxmgangdata 
    #dictargs['numbsidecart'] = numbsidecart 
    #dictargs['initnumbelempop0reg0'] = 1
    #dictargs['maxmnumbelempop0reg0'] = 1
    #dictargs['shrtfram'] = False
    dictargs['numbswep'] = 1
    dictargs['numbsamp'] = 1
    #dictargs['verbtype'] = 2
    dictargs['optitype'] = 'none'
    dictargs['elemspatevaltype'] = ['full']
    # temp
    dictargs['priofactdoff'] = 0.
    
    #dictargs['numbsamp'] = 1
    
    #listnamecnfgextn = ['home2msc0600', 'home4msc0600', 'home7msc0600']
    #listnamecnfgextn = ['extr2msc0300', 'extr4msc0300', 'home2msc0300', 'home4msc0300', 'home7msc0300']
    listnamecnfgextn = ['home2msc0300', 'home4msc0300', 'home7msc0300', 'home2msc0600', 'home4msc0600', 'home7msc0600']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
   
    for namecnfgextn in listnamecnfgextn:
        numbsidecart, strgexpo, strgexprsbrt, namestat, anlytype = retr_argschan(namecnfgextn[:4], namecnfgextn[4:8], int(namecnfgextn[8:]))
        
        # temp
        if '2msc' in anlytype:
            strgexpo = 500. * 2e6
        if '4msc' in anlytype:
            strgexpo = 500. * 4e6
        if '7msc' in anlytype:
            strgexpo = 500. * 7e6
        maxmgangdata = 0.492 / anglfact * numbsidecart / 2.
        dictargsvari[anlytype]['anlytype'] = anlytype
        dictargsvari[anlytype]['maxmgangdata'] = maxmgangdata 
        dictargsvari[anlytype]['numbsidecart'] = numbsidecart
        dictargsvari[anlytype]['strgexpo'] = strgexpo
        dictargsvari[anlytype]['strgexprsbrt'] = strgexprsbrt
        if namecnfgextn[:8] == 'extr4msc' or namecnfgextn[:8] == 'home7msc':
            dictargsvari[anlytype]['savestat'] = True
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
   

def retr_argschan(datatype, strgexpomaps, numbsidecart):
    
    anlytype = datatype + strgexpomaps + '%04d' % numbsidecart
    namestat = 'pcat_chan_inpt_' + anlytype
    strgexpo = 'expochan%s.fits' % anlytype
    strgexprsbrt = 'sbrtchan%s.fits' % anlytype

    return numbsidecart, strgexpo, strgexprsbrt, namestat, anlytype


globals().get(sys.argv[1])(*sys.argv[2:])
