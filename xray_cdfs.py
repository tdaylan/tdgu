from __init__ import *
from astropy.coordinates import SkyCoord
from pcat.util import retr_refrchaninit

def writ_chan():

    print 'Writing CDF-S dataset for PCAT...'
    
    pixlsize = deg2rad(0.492 / 3600.)
    apix = pixlsize**2
    anglfact = 3600. * 180. / pi
            
    numbevtt = 1
  
    #listnumbsidetile = [100, 300]
    #listnumbsidecntr = [600, 300]
    listnumbsidetile = [600, 300, 100]
    listnumbsidecntr = [600, 300, 600]
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

        for k in range(numbmaps):
            
            print 'expomaps[k]'
            print expomaps[k]
                
            if k == 0:
                #mapstemp = rand(600**2).reshape((600, 600))
                if datatype == 'extr':
                    # count map
                    if expomaps[k] == 2:
                        path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
                    elif expomaps[k] == 4:
                        path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
                    mapstemp = pf.getdata(path, 0)
                if datatype == 'home':
                    # count map
                    path = '/n/fink1/rfeder/xray_pcat/merged_flux_10_02/%dMs/%s-%s_flux.img' % (expomaps[k], strgener[0], strgener[1])
                    mapstemp = pf.getdata(path, 0)
            
                numbsideyaxi = mapstemp.shape[0]
                numbsidexaxi = mapstemp.shape[1]
                indxcntrxaxi = numbsidexaxi / 2
                indxcntryaxi = numbsideyaxi / 2
                print 'numbsideyaxi'
                print numbsideyaxi
                print 'numbsidexaxi'
                print numbsidexaxi
            
            for p, (numbsidetile, numbsidecntr) in enumerate(zip(listnumbsidetile, listnumbsidecntr)):

                print 'numbsidetile'
                print numbsidetile
                print 'numbsidecntr'
                print numbsidecntr
                
                facttile = numbsidecntr / numbsidetile
                numbtile = facttile**2

                if numbtile == 1:
                    numbpixloffsyaxi = (numbsideyaxi - numbsidecntr) / 2
                    numbpixloffsxaxi = (numbsidexaxi - numbsidecntr) / 2
                else:
                    numbpixloffsyaxi = 0
                    numbpixloffsxaxi = 0
                print 'facttile'
                print 'numbtile'
                print numbtile
                print facttile
                print

                figr, axis = plt.subplots(figsize=(12, 12))
                cntpfull = zeros((numbener, numbsidecntr, numbsidecntr, numbevtt))
                for t in range(numbtile):
                    
                    if numbtile == 1:
                        strgtile = 'none'
                    else:
                        strgtile = '%04d' % t
                    strgmaps = '%s%dmsc%04d%s' % (datatype, expomaps[k], numbsidecntr, strgtile)

                    # determine map shape
                    print 'numbpixloffsxaxi'
                    print numbpixloffsxaxi
                    print 'numbpixloffsyaxi'
                    print numbpixloffsyaxi
                    minmindxxaxi = (t // facttile) * numbsidetile + numbpixloffsxaxi
                    maxmindxxaxi = (t // facttile + 1) * numbsidetile + numbpixloffsxaxi
                    minmindxyaxi = (t % facttile) * numbsidetile + numbpixloffsyaxi
                    maxmindxyaxi = (t % facttile + 1) * numbsidetile + numbpixloffsyaxi
                     
                    print 'minmindxxaxi'
                    print minmindxxaxi
                    print 'maxmindxxaxi'
                    print maxmindxxaxi
                    print 'minmindxyaxi'
                    print minmindxyaxi
                    print 'maxmindxyaxi'
                    print maxmindxyaxi
                    print

                    cntp = zeros((numbener, numbsidetile, numbsidetile, numbevtt))
                    expo = zeros((numbener, numbsidetile, numbsidetile, numbevtt))
                    cntpback = empty((numbener, numbsidetile, numbsidetile, numbevtt))
                    sbrt = zeros_like(cntp)
                    sbrtback = zeros_like(cntp)
                    
                    if datatype == 'extr':
                        if expomaps[k] == 2:
                            path = pathdata + 'CDFS-2Ms-0p5to2-asca-im-bin1-astwk.fits'
                        elif expomaps[k] == 4:
                            path = pathdata + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
                        cntp[0, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]

                        ## hard band
                        if expomaps[k] == 2:
                            path = pathdata + 'CDFS-2Ms-2to8-asca-im-bin1-astwk.fits'
                        elif expomaps[k] == 4:
                            path = pathdata + 'CDFS-4Ms-2to8-asca-im-bin1.fits'
                        cntp[1, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]

                        # exposure
                        ## soft band
                        if expomaps[k] == 2:
                            path = pathdata + 'CDFS-2Ms-0p5to2-bin1-astwk.emap'
                        elif expomaps[k] == 4:
                            path = pathdata + 'CDFS-4Ms-0p5to2-bin1.emap'
                        expo[0, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                        ## hard band
                        if expomaps[k] == 2:
                            path = pathdata + 'CDFS-2Ms-2to8-bin1-astwk.emap'
                        elif expomaps[k] == 4:
                            path = pathdata + 'CDFS-4Ms-2to8-bin1.emap'
                        expo[1, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                    
                        # background
                        ## soft band
                        if expomaps[k] == 2:
                            path = pathdata + 'CDFS-4Ms-0p5to2-bin1.back'
                        elif expomaps[k] == 4:
                            path = pathdata + 'CDFS-2Ms-0p5to2-asca-bkg-bin1.fits'
                        cntpback[0, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                        ## hard band
                        if expomaps[k] == 2:
                            path = pathdata + 'CDFS-4Ms-2to8-bin1.back'
                        elif expomaps[k] == 4:
                            path = pathdata + 'CDFS-2Ms-2to8-asca-bkg-bin1.fits'
                        cntpback[1, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                    
                        for i in indxener:
                            indxtemp = where(expo[i, :, :, 0] > 0.)
                            sbrtback[i, indxtemp[0], indxtemp[1], 0] = cntpback[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
                    
                    if datatype == 'home':
                        strgvarb = ['thresh.expmap', 'thresh.img']
                        strgvarbmine = ['expo', 'sbrt']
                        for i in indxener:
                            for a in range(2):
                                if a == 0:
                                    path = '/n/fink1/rfeder/xray_pcat/cdfs/merged_%sMs/rest_fov/%d/%s-%s_thresh.expmap' % (expomaps[k], i, strgener[i], strgener[i+1])
                                    
                                    #expo[i, :, :, 0] = zeros((600, 600))[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                                    expo[i, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                                if a == 1:
                                    path = '/n/fink1/rfeder/xray_pcat/merged_flux_10_02/%dMs/%s-%s_flux.img' % (expomaps[k], strgener[i], strgener[i+1])
                                    cntp[i, :, :, 0] = pf.getdata(path, 0)[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                                    #cntp[i, :, :, 0] = zeros((600, 600))[minmindxyaxi:maxmindxyaxi, minmindxxaxi:maxmindxxaxi]
                                    cntp[i, :, :, 0] *= expo[i, :, :, 0]

                    cntpfull[:, minmindxyaxi-numbpixloffsyaxi:maxmindxyaxi-numbpixloffsyaxi, minmindxxaxi-numbpixloffsxaxi:maxmindxxaxi-numbpixloffsxaxi, :] = cntp
                    
                    axis.axvline(pixlsize * (minmindxxaxi - numbsidecntr / 2) * anglfact, lw=1)
                    axis.axhline(pixlsize * (minmindxyaxi - numbsidecntr / 2) * anglfact, lw=1)
                    axis.axvline(pixlsize * (maxmindxxaxi - numbsidecntr / 2) * anglfact, lw=1)
                    axis.axhline(pixlsize * (maxmindxyaxi - numbsidecntr / 2) * anglfact, lw=1)
                    
                    for i in indxener:
                        indxtemp = where(expo[i, :, :, 0] > 0.)
                        sbrt[i, indxtemp[0], indxtemp[1], 0] = cntp[i, indxtemp[0], indxtemp[1], 0] / expo[i, indxtemp[0], indxtemp[1], 0] / diffener[i] / apix
                    
                    #if True:
                    if datatype == 'home':
                        for i in indxener:
                            print 'i'
                            print i
                            print 'expo[i, :, :, 0]'
                            summgene(expo[i, :, :, 0])
                            print 'cntp[i, :, :, 0]'
                            summgene(cntp[i, :, :, 0])
                            print 'sbrt[i, :, :, 0]'
                            summgene(sbrt[i, :, :, 0])
                            if amax(cntp[i, :, :, 0]) == 0 and datatype == 'home':
                                print 'datatype'
                                print datatype
                                print 'indxener'
                                print indxener
                                raise Exception('')
                            print

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
                
                extt = anglfact * numbsidecntr * array([-1, 1, -1, 1]) * pixlsize / 2.
                imag = axis.imshow(arcsinh(cntpfull[0, :, :, 0]), extent=extt, cmap='Greys', origin='lower', interpolation='none')
                axis.set_xlabel(r'$\theta_1$ [$^\circ$]')
                axis.set_ylabel(r'$\theta_2$ [$^\circ$]')
                plt.colorbar(imag, ax=axis, fraction=0.03)
                plt.tight_layout()
                path = os.environ["TDGU_DATA_PATH"] + '/xray_back/imag/'
                strgmaps = '%s%dmsc%04d' % (datatype, expomaps[k], p)
                figr.savefig(path + 'cntpfull' + strgmaps + '.pdf')

# test suites

def pcat_chan_mock_spmr(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 15
    
    maxmgangdata = 0.492 / anglfact * numbsidecart / 2.

    dictargs = {}
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
    
    listnamecnfgextn = ['free', 'nomi', 'parsnomi', 'brgt', 'fain', 'psfn']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    #dictargsvari['free']['minmnumbelempop0reg0'] = 0.
    #dictargsvari['free']['maxmnumbelempop0reg0'] = 0.
    dictargsvari['free']['inittype'] = 'rand'
    dictargsvari['free']['probtran'] = 0.4
    dictargsvari['free']['probspmr'] = 0.3
    
    dictargsvari['parsnomi']['priofactdoff'] = 1.
    
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
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_chan_mock_popl(strgcnfgextnexec=None):
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              exprtype='chan', \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              numbpopl=2, \
                              numbelempop0reg0=50, \
                              numbelempop1reg0=50, \
                             )


# science suites
def pcat_chan_mock(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['trueelemtype'] = ['lghtpntsagnntrue']
    dictargs['fittelemtype'] = ['lghtpntsagnnassc']
    dictargs['priofactdoff'] = 0.
    
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 100
    
    #dictargs['numbswep'] = 1000
    #dictargs['numbsamp'] = 10

    dictargs['truenumbelempop0reg0'] = 100
    dictargs['maxmnumbelempop0reg0'] = 400
    
    listnamecnfgextn = ['nomi']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['nomi']['checprio'] = True
                
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                    
                                  forcprev=True, \
                                  #execpara=True, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_chan_mock_trueminmdlos(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['trueelemtype'] = ['lghtpntsagnntrue']
    dictargs['fittelemtype'] = ['lghtpntsagnnassc']
    dictargs['priofactdoff'] = 0.
    
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['maxmnumbelempop0reg0'] = 400
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['loww']['trueminmdlos'] = 3e6
    dictargsvari['nomi']['trueminmdlos'] = 1e7
    dictargsvari['high']['trueminmdlos'] = 3e7
                
    lablxaxi = r'$d_{los,min}$ [Mpc]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(1e-6 * dictargsvari[namecnfgextn]['trueminmdlos']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  #execpara=True, \
                                  
                                  namexaxi='trueminmdlos', \
                                  lablxaxi=lablxaxi, \
                                  scalxaxi=scalxaxi, \
                                  listtickxaxi=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_chan_mock_truemaxmdlos(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['trueelemtype'] = ['lghtpntsagnntrue']
    dictargs['fittelemtype'] = ['lghtpntsagnnassc']
    dictargs['priofactdoff'] = 0.
    
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['maxmnumbelempop0reg0'] = 400
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['loww']['truemaxmdlos'] = 3e8
    dictargsvari['nomi']['truemaxmdlos'] = 1e9
    dictargsvari['high']['truemaxmdlos'] = 3e9
                
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  #execpara=True, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_chan_mock_trueminmlum0(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['trueelemtype'] = ['lghtpntsagnntrue']
    dictargs['fittelemtype'] = ['lghtpntsagnnassc']
    dictargs['priofactdoff'] = 0.
    
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['maxmnumbelempop0reg0'] = 400
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['loww']['trueminmlum0'] = 3e42
    dictargsvari['nomi']['trueminmlum0'] = 1e43
    dictargsvari['high']['trueminmlum0'] = 3e43
                
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  #execpara=True, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_chan_mock_truemaxmlum0(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['strgexpo'] = 'expochanhome7msc06000000.fits'
    dictargs['trueelemtype'] = ['lghtpntsagnntrue']
    dictargs['fittelemtype'] = ['lghtpntsagnnassc']
    dictargs['priofactdoff'] = 0.
    
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['maxmnumbelempop0reg0'] = 400
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['loww']['truemaxmlum0'] = 3e45
    dictargsvari['nomi']['truemaxmlum0'] = 1e46
    dictargsvari['high']['truemaxmlum0'] = 3e46
                
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  #execpara=True, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_chan_mock_mapo(strgcnfgextnexec=None):
    
    datatype = 'home'
    strgexpomaps = '4msc'
    numbsidecart = 300
    gridchan = pcat.main.init( \
                              evoltype='maxmllik', \
                              strgexpo='expochan%s%s%04d.fits' % (datatype, strgexpomaps, numbsidecart), \
                              exprtype='chan', \
                             )


def pcat_chan_inpt(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    
    # temp
    dictargs['inittype'] = 'reco'
    #dictargs['anlytypedata'] = maxmgangdata 
    #dictargs['numbsidecart'] = numbsidecart 
    #dictargs['initnumbelempop0reg0'] = 1
    #dictargs['maxmnumbelempop0reg0'] = 1
    #dictargs['propcomp'] = False
    #dictargs['probtran'] = 0.
    #dictargs['spectype'] = ['colr']
    
    if os.uname()[1] == 'fink1.rc.fas.harvard.edu' or os.uname()[1] == 'fink2.rc.fas.harvard.edu':
        dictargs['rtagmock'] = '20180301_131911_pcat_chan_mock_nomi_1000000'
    else:
        dictargs['rtagmock'] = '20180205_184023_pcat_chan_mock_nomi_100000'
    dictargs['savestat'] = True
    dictargs['priofactdoff'] = 0.
    
    listnamecnfgextn = []
    listnamecnfgextn += ['home2msc0300none', 'home4msc0300none', 'home7msc0300none']
    #listnamecnfgextn += ['home7msc0600none']
    #for k in range(36):
    #    listnamecnfgextn.append('home7msc0600%04d' % k)

    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
   
    for namecnfgextn in listnamecnfgextn:
        numbsidecart, strgexpo, strgexprsbrt, namestat, anlytype = retr_argschan(namecnfgextn[:4], namecnfgextn[4:8], int(namecnfgextn[8:12]), namecnfgextn[12:16])
   
        dictargsvari[anlytype]['namerecostat'] = 'pcat_chan_inpt_' + namecnfgextn
    
        if namecnfgextn == 'home7msc0300none':
            dictargsvari[anlytype]['checprio'] = True
        
        # temp
        if namecnfgextn[8:] == '0600none':
            dictargsvari[anlytype]['numbsamp'] = 1
            dictargsvari[anlytype]['optitype'] = 'none'
            dictargsvari[anlytype]['elemspatevaltype'] = ['full']
    
        if namecnfgextn[8:12] == '0600' and namecnfgextn[12:16].isdigit():
            dictargsvari[anlytype]['maxmnumbelempop0reg0'] = 30
            
        if namecnfgextn[12:16] == 'none':
            maxmgangdata = 0.492 / anglfact * numbsidecart / 2.
        else:
            maxmgangdata = 0.492 / anglfact * 100. / 2.
        
        dictargsvari[anlytype]['anlytype'] = anlytype
        dictargsvari[anlytype]['maxmgangdata'] = maxmgangdata 
        dictargsvari[anlytype]['strgexpo'] = strgexpo
        dictargsvari[anlytype]['strgexprsbrt'] = strgexprsbrt
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  execpara=True, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
   

def retr_argschan(datatype, strgexpomaps, numbsidecart, strgtile):
    
    anlytype = datatype + strgexpomaps + '%04d' % numbsidecart + strgtile
    namestat = 'pcat_chan_inpt_' + anlytype
    strgexpo = 'expochan%s.fits' % anlytype
    strgexprsbrt = 'sbrtchan%s.fits' % anlytype

    return numbsidecart, strgexpo, strgexprsbrt, namestat, anlytype


globals().get(sys.argv[1])(*sys.argv[2:])
