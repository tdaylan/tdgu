from __init__ import *

def writ_ferm_raww():
    
    gdat = tdpy.util.gdatstrt()
    
    gdat.test = False
        
    #defn_gtbn()

    gdat.pathdata = os.environ["TDGU_DATA_PATH"] + '/ferm_igal/data/'
    #gdat.recotype = ['rec7', 'rec7', 'rec8', 'rec8', 'rec8']
    #gdat.enertype = ['pnts', 'pnts', 'pnts', 'back', 'back']
    #gdat.strgtime = ['tmin=239155201 tmax=364953603', 'tmin=INDEF tmax=INDEF', 'tmin=INDEF tmax=INDEF', 'tmin=INDEF tmax=INDEF', 'tmin=INDEF tmax=INDEF']
    #gdat.numbside = [256, 256, 256, 256, 512]
    
    gdat.recotype = ['rec8', 'rec8']
    gdat.enertype = ['pnts', 'back']
    gdat.strgtime = ['tmin=INDEF tmax=INDEF', 'tmin=INDEF tmax=INDEF']
    gdat.numbside = [256, 256]
    
    #gdat.recotype = ['rec8']
    #gdat.enertype = ['pnts']
    #gdat.strgtime = ['tmin=INDEF tmax=INDEF']
    #gdat.numbside = [256]
    
    numbproc = len(gdat.recotype)
    indxproc = arange(numbproc)

    gdat.evtc = []
    gdat.photpath = []
    gdat.weekinit = []
    gdat.weekfinl = []
    for n in range(numbproc):
        if gdat.recotype[n] == 'rec7':
            gdat.evtc.append(2)
            gdat.photpath.append('p7v6c')
            gdat.weekinit.append(9)
            gdat.weekfinl.append(11)
            gdat.weekfinl.append(5000)
        if gdat.recotype[n] == 'rec8':
            gdat.evtc.append(128)
            gdat.photpath.append('photon')
            gdat.weekinit.append(11)
            #gdat.weekfinl.append(12)
            gdat.weekfinl.append(700)
   
    gdat.strgener = [[] for k in indxproc]
    for k in indxproc:
        if gdat.enertype[k] == '3fgl':
            gdat.strgener[k] = 'ebinalg=FILE ebinfile=$TDGU_DATA_PATH/ferm_igal/data/gtbndefn_%s.fits' % gdat.enertype[k]
        if gdat.enertype[k] == 'pnts':
            gdat.strgener[k] = 'ebinalg=LOG EMIN=300 EMAX=10000 ENUMBINS=3'
        if gdat.enertype[k] == 'back':
            gdat.strgener[k] = 'ebinalg=LOG EMIN=300 EMAX=300000 ENUMBINS=30'
    
    indxproc = arange(numbproc)
    
    if numbproc == 1 or gdat.test:
        writ_ferm_raww_work(gdat, 0)
    else:
        # process pool
        pool = mp.Pool(numbproc)

        # spawn the processes
        writ_ferm_part = functools.partial(writ_ferm_raww_work, gdat)
        pool.map(writ_ferm_part, indxproc)
        pool.close()
        pool.join()


def writ_ferm_raww_work(gdat, indxprocwork):

    rtag = '%s%s%04d' % (gdat.recotype[indxprocwork], gdat.enertype[indxprocwork], gdat.numbside[indxprocwork])
    
    # make file lists
    infl = gdat.pathdata + 'phot_%s.txt' % rtag
    spac = gdat.pathdata + 'spac_%s.txt' % rtag
    
    gdat.evtt, gdat.numbevtt, gdat.indxevtt = tdpy.util.retr_evttferm(gdat.recotype[indxprocwork])
        
    numbweek = (gdat.weekfinl[indxprocwork] - gdat.weekinit[indxprocwork])
    listweek = floor(linspace(gdat.weekinit[indxprocwork], gdat.weekfinl[indxprocwork] - 1, numbweek)).astype(int)
    cmnd = 'rm -f ' + infl
    os.system(cmnd)
    cmnd = 'rm -f ' + spac
    os.system(cmnd)
    for week in listweek:
        cmnd = 'ls -d -1 $FERMI_DATA/weekly/spacecraft/*_w%03d_* >> ' % week + spac
        os.system(cmnd)
        cmnd = 'ls -d -1 $FERMI_DATA/weekly/%s/*_w%03d_* >> ' % (gdat.photpath[indxprocwork], week) + infl
        os.system(cmnd)
    for m in gdat.indxevtt:

        if gdat.recotype[indxprocwork] == 'rec7':
            strgirfn = 'P7REP_SOURCE_V15'
        if gdat.recotype[indxprocwork] == 'rec8':
            strgirfn = 'P8R2_SOURCE_V6'

        thisevtt = gdat.evtt[m]
        if gdat.recotype[indxprocwork] == 'rec7':
            strgpsfn = 'convtype=%d' % thisevtt
        if gdat.recotype[indxprocwork] == 'rec8':
            strgpsfn = 'evtype=%d' % thisevtt
         
        sele = gdat.pathdata + 'fermsele%04d%s.fits' % (thisevtt, rtag)
        filt = gdat.pathdata + 'fermfilt%04d%s.fits' % (thisevtt, rtag)
        live = gdat.pathdata + 'fermlive%04d%s.fits' % (thisevtt, rtag)
        cntp = gdat.pathdata + 'cntpferm%04d%s.fits' % (thisevtt, rtag)
        expo = gdat.pathdata + 'expoferm%04d%s.fits' % (thisevtt, rtag)
        psfn = gdat.pathdata + 'psfnferm%04d%s.fits' % (thisevtt, rtag)

        cmnd = 'gtselect infile=' + infl + ' outfile=' + sele + ' ra=INDEF dec=INDEF rad=INDEF ' + \
            gdat.strgtime[indxprocwork] + ' emin=300 emax=300000 zmax=90 evclass=%d %s' % (gdat.evtc[indxprocwork], strgpsfn)
        
        if os.path.isfile(cntp) and os.path.isfile(expo) and not gdat.test:
            continue
     
        print cmnd
        print ''
        if not gdat.test and not os.path.isfile(sele):
            os.system(cmnd)

        cmnd = 'gtmktime evfile=' + sele + ' scfile=' + spac + ' filter="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52"' + ' outfile=' + filt + ' roicut=no'
        print cmnd
        print ''
        if not gdat.test and not os.path.isfile(filt):
            os.system(cmnd)

        cmnd = 'gtbin evfile=' + filt + ' scfile=NONE outfile=' + cntp + \
            ' %s ' % gdat.strgener[indxprocwork] + \
            'algorithm=HEALPIX hpx_ordering_scheme=RING coordsys=GAL hpx_order=%d hpx_ebin=yes' % log2(gdat.numbside[indxprocwork])
        print cmnd
        print ''
        if not gdat.test and not os.path.isfile(cntp):
            os.system(cmnd)

        cmnd = 'gtltcube evfile=' + filt + ' scfile=' + spac + ' outfile=' + live + ' dcostheta=0.025 binsz=1'
        print cmnd
        print ''
        if not gdat.test and not os.path.isfile(live):
            os.system(cmnd)

        cmnd = 'gtexpcube2 infile=' + live + ' cmap=' + cntp + ' outfile=' + expo + ' irfs=CALDB evtype=%03d bincalc=CENTER' % thisevtt
        print cmnd
        print ''
        if not gdat.test and not os.path.isfile(expo):
            os.system(cmnd)

    #cmnd = 'rm %s %s %s %s %s' % (infl, spac, sele, filt, live)
    #os.system(cmnd)


def writ_ferm():

    recotype = 'rec8'
    regitype = 'ngal'
    numbside = 256
    enertype = 'pnts'
    
    pathinpt = os.environ["TDGU_DATA_PATH"] + '/ferm_igal/data/'
    pathoutp = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'

    if enertype == 'back':
        numbener = 30
        minmener = 0.3
        maxmener = 300.
        binsener = logspace(log10(minmener), log10(maxmener), numbener + 1)
    elif enertype == 'pnts':
        numbener = 3
        minmener = 0.3
        maxmener = 10.
        binsener = logspace(log10(minmener), log10(maxmener), numbener + 1)
    else:
        binsener = array([0.1, 0.3, 1., 3., 10., 100.])
        numbener = binsener.size - 1
    indxener = arange(numbener)
    diffener = binsener[1:] - binsener[:-1]

    numbside = 256
    evtt, numbevtt, indxevtt = tdpy.util.retr_evttferm(recotype)
    
    evtt = [16, 32]
    numbevtt = 2
    indxevtt = [0, 1]

    numbpixl = 12 * numbside**2
    apix = 4. * pi / numbpixl

    cntp = zeros((numbener, numbpixl, numbevtt))
    expo = zeros((numbener, numbpixl, numbevtt))
    sbrt = zeros((numbener, numbpixl, numbevtt))
    
    for m in indxevtt:
        
        thisevtt = evtt[m]

        # m == 1 exposure files are corrupted!
        if thisevtt == 8:
            continue

        path = pathinpt + '/expoferm%04d%s%s%04d.fits' % (thisevtt, recotype, enertype, numbside)
        listhdun = pf.open(path)

        expoarry = pf.getdata(path, 1)
        for i in indxener:
            expo[i, :, m] = expoarry['ENERGY%d' % (i + 1)]

        path = pathinpt + '/cntpferm%04d%s%s%04d.fits' % (thisevtt, recotype, enertype, numbside)
        cntparry = pf.getdata(path)
        for i in indxener:
            cntp[i, :, m] = cntparry['CHANNEL%d' % (i + 1)]

    indxexpo = where(expo > 0.) 
    sbrt[indxexpo] = cntp[indxexpo] / expo[indxexpo] / apix
    sbrt /= diffener[:, None, None]

    if regitype == 'ngal':
        for i in indxener:
            for m in indxevtt:
                
                if recotype == 'rec7':
                    if m < 2:
                        continue
                    elif m == 2:
                        thisevtt = 2
                    elif m == 3:
                        thisevtt = 1
                else:
                    thisevtt = evtt[m]

                almc = hp.map2alm(sbrt[i, :, m])
                hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                sbrt[i, :, m] = hp.alm2map(almc, numbside)

                almc = hp.map2alm(expo[i, :, m])
                hp.rotate_alm(almc, 0., 0.5 * pi, 0.)
                expo[i, :, m] = hp.alm2map(almc, numbside)
    
    for i in indxener:
        for m in indxevtt:
            if (sbrt[i, :, m] == 0.).all():
                print 'im'
                print i, m
                raise Exception('')
            if (expo[i, :, m] == 0.).all():
                print 'im'
                print i, m
                raise Exception('')

    path = pathoutp + '/expoferm%s%s%s%04d.fits' % (recotype, enertype, regitype, numbside)
    print 'Writing to %s...' % path
    print 'expo'
    summgene(expo)
    print
    pf.writeto(path, expo, clobber=True)

    path = pathoutp + '/sbrtferm%s%s%s%04d.fits' % (recotype, enertype, regitype, numbside)
    print 'Writing to %s...' % path
    print 'sbrt'
    summgene(sbrt)
    print
    pf.writeto(path, sbrt, clobber=True)


def retr_plnkmapsorig(gdat, strgmapsplnk):

    if strgmapsplnk == 'radi':
        mapsplnkorig = pf.getdata(gdat.pathdatatdgu + 'plnk/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 1)['RADIANCE']
        mapsplnkorig = hp.reorder(mapsplnkorig, n2r=True)
    else:
        
        mapsplnk = tdpy.util.retr_mapsplnkfreq(strgmapsplnk)

        # plot the Planck map
        tdpy.util.plot_maps(gdat.pathimag + 'mapsplnk%s.pdf' % strgmapsplnk, mapsplnk, satu=True)

        if gdat.subspnts:
            print 'Subtracting point sources...'

            # subtract PSs from the Planck maps
            ## read PCCS
            if strgmapsplnk[1] == '0':
                strg = 'R2.04'
            else:
                strg = 'R2.01'
            # temp
            dataplnk = pf.getdata(gdat.pathdata + 'plnk/COM_PCCS_%s_%s.fits' % (strgmapsplnk[1:], strg), 1)
            fluxpntsplnk = dataplnk['GAUFLUX'] * 1e-3 # [Jy]
            lgalpntsplnk = dataplnk['GLON'] # [deg]
            lgalpntsplnk = (lgalpntsplnk - 180.) % 360. - 180.
            bgalpntsplnk = dataplnk['GLAT'] # [deg]
            fwhmpntsplnk = dataplnk['GAU_FWHM_EFF'] / 60. # [deg]
            stdvpntsplnk = fwhmpntsplnk / 2. / sqrt(2. * log(2.)) # [deg]
            
            # filter out bad 
            indxpntsgood = where((stdvpntsplnk >= 0.) & (fluxpntsplnk > 0.))[0]
            lgalpntsplnk = lgalpntsplnk[indxpntsgood]
            bgalpntsplnk = bgalpntsplnk[indxpntsgood]
            fluxpntsplnk = fluxpntsplnk[indxpntsgood]
            stdvpntsplnk = stdvpntsplnk[indxpntsgood]
            
            # sort PS with respect to flux
            indxpntsplnk = argsort(fluxpntsplnk)
            lgalpntsplnk = lgalpntsplnk[indxpntsplnk]
            bgalpntsplnk = bgalpntsplnk[indxpntsplnk]
            fluxpntsplnk = fluxpntsplnk[indxpntsplnk]
            stdvpntsplnk = stdvpntsplnk[indxpntsplnk]
            
            # temp
            numbpntskeep = 3
            lgalpntsplnk = lgalpntsplnk[0:numbpntskeep]
            bgalpntsplnk = bgalpntsplnk[0:numbpntskeep]
            fluxpntsplnk = fluxpntsplnk[0:numbpntskeep]
            stdvpntsplnk = stdvpntsplnk[0:numbpntskeep]
            
            numbpntsplnk = fluxpntsplnk.size
            
            print 'Using %d PS from the PCCS...' % numbpntsplnk

            ## calculate PS map using PCCS
            numbsidepnts = int(sqrt(mapsplnk.size / 12))
            pathmapspntsplnk = gdat.pathdata + 'mapspntsplnk%s%04d.fits' % (strgmapsplnk, numbsidepnts)
            if os.path.isfile(pathmapspntsplnk):
                print 'Reading %s...' % pathmapspntsplnk
                mapspntsplnk = pf.getdata(pathmapspntsplnk)
            else:
                mapspntsplnk = tdpy.util.retr_mapspnts(lgalpntsplnk, bgalpntsplnk, stdvpntsplnk, fluxpntsplnk, verbtype=2, numbside=numbsidepnts)
                ## plot the PCSC map
                tdpy.util.plot_maps(gdat.pathimag + 'mapspntsplnk%s.pdf' % strgmapsplnk, mapspntsplnk, satu=True)
                pf.writeto(pathmapspntsplnk, mapspntsplnk, clobber=True)
            mapsorigplnk = mapsplnk - mapspntsplnk
        else:
            mapsorigplnk = mapsplnk

    return mapsorigplnk


def defn_gtbn():
    
    listenertype = ['pnts', 'back']
    for enertype in listenertype:
        # this does not work with healpix
        if enertype == 'ferm':
            lowrener = array([0.1, 0.3, 1., 3., 10.])
            upprener = array([0.3, 1., 3., 10., 100.])
        if enertype == 'pnts' or enertype == 'back':
            if enertype == 'pnts':
                numbener = 3
                minmener = 0.3
                maxmener = 10.
            if enertype == 'back':
                numbener = 30
                minmener = 0.1
                maxmener = 100.
            binsener = logspace(log10(minmener), log10(maxmener), numbener + 1)
            lowrener = binsener[:-1]
            upprener = binsener[1:]
        limtener = stack((lowrener, upprener), axis=1)
        pathinpt = os.environ["TDGU_DATA_PATH"] + '/ferm_igal/data/gtbndefn_%s.dat' % enertype
        pathoutp = os.environ["TDGU_DATA_PATH"] + '/ferm_igal/data/gtbndefn_%s.fits' % enertype
        print 'Writing to %s...' % pathinpt
        savetxt(pathinpt, limtener, fmt='%10.5g')
        os.system('gtbindef E %s GeV' % (pathoutp))


def merg_maps(numbside=256, mpolmerg=180., mpolsmth=360., strgmaps='radi'):

    # construct the global object
    gdat = tdpy.util.gdatstrt()
    
    # get the time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # run tag
    rtag = strgtimestmp + '_%04d_%04d_%04d' % (numbside, mpolmerg, mpolsmth)
    
    # paths
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'ferm_igal/', 'ferm_igal', rtag)

    # time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()

    timeinit = time.time()

    calcfactconv = False
    gdat.subspnts = True
   
    # analysis setup
    ## plots
    alph = 0.5
    plotfull = False

    ## Healpix grid
    lgalheal, bgalheal, numbpixl, apix = tdpy.util.retr_healgrid(numbside)

    ## axes
    ### Fermi-LAT energy
    binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    meanener = sqrt(binsener[1:] * binsener[:-1])
    numbener = meanener.size
    indxener = arange(numbener)

    ## constants
    consplnk = 6.63e-34 # [J s]
    consbolt = 1.38e-23 # [J/K]
    tempcmbr = 2.725 # [K]
    velolght = 3e8 # [m/s]

    ## multipole
    maxmmpol = 3. * numbside - 1.
    numbalmc = int(maxmmpol * (maxmmpol + 1.) / 2. + maxmmpol + 1)
    numbmpol = int(maxmmpol) + 1
    mpol = arange(numbmpol)
    mpolgrid, temp = hp.Alm.getlm(lmax=maxmmpol)

    # read unit conversion data provided by Planck
    factconvplnk = loadtxt(gdat.pathdata + 'plnkunitconv.dat')
    
    ## Fermi-LAT flux map
    path = gdat.pathdata + '/sbrtfermrec8pntsigal0256.fits'
    mapsfermorig = sum(pf.getdata(path), 2)
    numbpixlferm = mapsfermorig.shape[1]
    numbsideferm = int(sqrt(numbpixlferm / 12))
    mapsfermorig -= mean(mapsfermorig, 1)[:, None]
    mapsfermorig /= std(mapsfermorig, 1)[:, None]
    mapsferm = empty_like(mapsfermorig)
    for i in indxener:
        tdpy.util.plot_maps(gdat.pathimag + 'mapsfermorig%04d.pdf' % i, mapsfermorig[i, :], satu=True)
        mapsferm[i, :] = tdpy.util.smth(mapsfermorig[i, :], mpolsmth, mpol=True)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsferm%04d.pdf' % i, mapsferm[i, :], satu=True)

    # 3FGL flux map
    path = gdat.pathdata + 'gll_psc_v16.fit'
    datafgl3 = pf.getdata(path)
    lgalfgl3 = datafgl3['glon']
    lgalfgl3 = ((lgalfgl3 - 180.) % 360.) - 180.
    bgalfgl3 = datafgl3['glat']
    stdvfgl3 = tdpy.util.retr_fwhmferm(meanener, indxevtt) / 2.
    specfgl3 = stack((datafgl3['Flux100_300'], datafgl3['Flux300_1000'], datafgl3['Flux1000_3000'], \
                                                            datafgl3['Flux3000_10000'], datafgl3['Flux10000_100000'])) / gdat.diffener[:, None]
    
    # temp
    numbpntsfgl3 = 2
    indxfgl3brgt = argsort(specfgl3[2, :])
    lgalfgl3 = lgalfgl3[indxfgl3brgt][:numbpntsfgl3]
    bgalfgl3 = bgalfgl3[indxfgl3brgt][:numbpntsfgl3]
    specfgl3 = specfgl3[:, indxfgl3brgt][:, :numbpntsfgl3]
    mapspntsferm = empty((numbener, numbpixlferm, numbevtt))
    for i in indxener:
        for m in indxevtt:
            mapspntsferm[i, :, m] = tdpy.util.retr_mapspnts(lgalfgl3, bgalfgl3, stdvfgl3[i, m], specfgl3[i, :], verbtype=2, numbside=numbsideferm)

    if plotfull:
        # plot tranmission spectra
        figr, axis = plt.subplots()
        axis.set_ylabel('$T$')
        axis.set_xlabel('$f$ [GHz]')
        axis.set_ylim([1e-4, 1.])
        axis.set_xlim([10., 2e3])
        for k in range(numbfreqplnk):
            labl = '%d' % int(strgmapsplnk[k])
            axis.loglog(1e-9 * freqband, tranband, label=labl)
        axis.legend(loc=3, ncol=2)
        plt.tight_layout()
        path = gdat.pathimag + 'tran.pdf'
        plt.savefig(path)
        plt.close(figr)
        strgmapsplnk = ['0030', '0044', '0070', '0100', '0143', '0217', '0353', '0545', '0857', 'radi']
        for k in range(numbmapsplnk):
            print 'Map number ', k
            print 'Maps string: ', strgmapsplnk[k]
            writ_plnk(strgmapsplnk[k])
    
    # get Planck PS mask
    if False:
        print 'Reading the Planck mask...'
        path = gdat.pathdatatdgu + 'plnk/HFI_Mask_PointSrc_2048_R2.00.fits'
        mapsmask = pf.open(path)[1].data['F353']
        mapsmask = hp.reorder(mapsmask, n2r=True)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsmask.pdf', mapsmask)
    
    strgmapsplnk = 'radi'
    #strgmapsplnk = '0857'
    
    # get input maps
    ## Planck map
    print 'Smoothing the Planck map...'
    path = gdat.pathdata + 'mapsplnk.fits'
    if os.path.isfile(path):
        print 'Reading %s...' % path
        mapsplnk = pf.getdata(path)
    else:
        mapsplnkorig = retr_plnkmapsorig(gdat, strgmapsplnk)
        mapsplnkorig -= mean(mapsplnkorig)
        mapsplnkorig /= std(mapsplnkorig)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsplnkorig.pdf', mapsplnkorig, satu=True)
        # temp
        indxpixlmask = []
        mapsplnk, mapsalmc, mpolsmthplnk, wghtsmthplnk = tdpy.util.smth(mapsplnkorig, mpolsmth, mpol=True, retrfull=True, indxpixlmask=indxpixlmask, numbsideoutp=numbside)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsplnk.pdf', mapsplnk)
        pf.writeto(path, mapsplnk, clobber=True)
    almcplnktemp = hp.map2alm(mapsplnk)

    numbmapsplnk = len(strgmapsplnk)

    ## PS
    # temp
    if False:
        mapspnts = zeros((2, numbpixl))
        numbpnts = array([100, 100000])
        for k in arange(numbpnts.size):
            stdv = zeros(numbpnts[k]) + 0.5
            flux = zeros(numbpnts[k]) + 1.
            lgal = (rand(numbpnts[k]) - 0.5) * 360.
            bgal = (rand(numbpnts[k]) - 0.5) * 360.
            mapspnts[k, :] = tdpy.util.retr_mapspnts(lgal, bgal, stdv, flux, numbside=numbside)

    ## Gaussian noise map
    mapsgaus = sqrt(0.25 / 30.) * randn(numbpixl)

    ## Fermi Diffuse Model
    # temp
    print 'Reading the Fermi diffuse model...'
    mapsfdfmorig = tdpy.util.retr_sbrtfdfm(binsener, numbside=numbside)
    mapsfdfmorig -= mean(mapsfdfmorig, 1)[:, None]
    mapsfdfmorig /= std(mapsfdfmorig, 1)[:, None]
    mapsfdfm = empty_like(mapsfdfmorig)
    almcfdfmorig = empty((numbener, numbalmc), dtype=complex) 
    almcfdfm = empty((numbener, numbalmc), dtype=complex) 
    for i in arange(numbener):
        almcfdfmorig[i, :] = hp.map2alm(mapsfdfmorig[i, :])
        tdpy.util.plot_maps(gdat.pathimag + 'mapsfdfmorig%04d.pdf' % i, mapsfdfmorig[i, :], satu=True)
        mapsfdfm[i, :], almcfdfm[i, :], mpolsmthfdfm, wghtsmthfdfm = tdpy.util.smth(mapsfdfmorig[i,:], mpolsmth, mpol=True, retrfull=True)

    # compute the weight
    wghtsing = exp(-0.5 * mpol * (mpol + 1.) / mpolmerg**2)
    wght = empty(numbalmc) 
    for l in mpol.astype(int):
        wght[where(mpolgrid == l)] = wghtsing[l]
    
    # plot the weight
    figr, axis = plt.subplots()
    axis.loglog(mpol, wghtsing, label='FDM')
    axis.loglog(mpol, 1. - wghtsing, label='Planck')
    axis.loglog(mpolsmthfdfm, wghtsmthfdfm, label='Smoothing Kernel')
    #axis.loglog(mpolsmthplnk + 0.1, wghtsmthplnk, label='Smoothing Kernel')

    axis.axvline(numbside, ls='-', color='black', alpha=alph, label='$N_{side}$')
    axis.axvline(mpolsmth, ls='--', color='black', alpha=alph, label='$l_{smth}$')
    axis.axvline(mpolmerg, ls='-.', color='black', alpha=alph, label='$l_{merg}$')
    axis.set_ylabel('$w_l$')
    axis.set_xlabel('$l$')
    axis.set_ylim([1e-4, 1.])
    axis.set_xlim([amin(mpol), amax(mpol)])
    axis.legend(loc=2)
    plt.tight_layout()
    path = gdat.pathimag + 'wght.pdf'
    plt.savefig(path)
    plt.close(figr)
   

    # compute power spectra
    print 'Computing power spectra...'
    ## power spectrum prefactor
    factmpol = mpol * (2. * mpol + 1.) / 4. / pi
    
    ## indices of multipoles over which the variance of Fermi diffuse model and Planck map is matched
    indxmpoltemp = where(mpol < 10.)[0]
    
    psecplnkuncr = factmpol * hp.anafast(mapsplnk)
    psecgaus = factmpol * hp.anafast(mapsgaus)
    
    psecfermorig = empty((numbener, numbmpol))
    psecferm = empty((numbener, numbmpol))
    psecfdfmorig = empty((numbener, numbmpol))
    psecfdfm = empty((numbener, numbmpol))
    psecplnk = empty((numbener, numbmpol))
    almcplnk = empty((numbener, numbalmc), dtype=complex)
    
    for i in arange(numbener):
        psecfermorig[i, :] = factmpol * hp.anafast(mapsfermorig[i, :])
        psecferm[i, :] = factmpol * hp.anafast(mapsferm[i, :])
        psecfdfmorig[i, :] = factmpol * hp.anafast(mapsfdfmorig[i, :])
        psecfdfm[i, :] = factmpol * hp.anafast(mapsfdfm[i, :])
        ## correct the Planck variance
        factcorr = sum(psecfdfm[i, indxmpoltemp]) / sum(psecplnkuncr[indxmpoltemp])
        psecplnk[i, :] = factcorr * psecplnkuncr
        almcplnk[i, :] = sqrt(factcorr) * almcplnktemp

    # merge the maps
    mapsmerg = empty((numbener, numbpixl))
    for i in arange(numbener):
        almcmerg = almcfdfm[i, :] * wght + almcplnk[i, :] * (1. - wght)
        mapsmerg[i, :] = hp.alm2map(almcmerg, numbside, verbose=False)
    
    # calculate the power spectrum of the merged map
    psecmerg = empty((numbener, numbmpol))
    for i in arange(numbener):
        psecmerg[i, :] = factmpol * hp.anafast(mapsmerg[i, :])
    
    # plot the power spectra
    for i in arange(numbener):
        figr, axis = plt.subplots()
        
        axis.loglog(mpol, psecfermorig[i, :], label='Fermi-LAT Native', color='b')
        axis.loglog(mpol, psecferm[i, :], label='Fermi-LAT Smooth', color='b', alpha=0.5)
        axis.loglog(mpol, psecfdfmorig[i, :], label='FDM Native', color='c')
        axis.loglog(mpol, psecfdfm[i, :], label='FDM Smooth', color='c', alpha=0.5)
        axis.loglog(mpol, psecplnkuncr, label='Planck Smooth', color='r', alpha=0.5)
        axis.loglog(mpol, psecplnk[i, :], label='Planck Smooth Norm.', color='r')
        axis.loglog(mpol, psecmerg[i, :], label='Merged', color='m')
        
        axis.loglog(mpol, psecgaus, label='Uncorrelated Noise', alpha=0.1, ls='--', color='black')
#        axis.loglog(mpol, 3e-3 * exp(mpol / 1800.), label='Planck beam deconvolved', alpha=0.1, ls='-.', color='black')
    
        axis.axvline(numbside, ls='-', color='black', alpha=0.6, label='$N_{side}$')
        axis.axvline(mpolsmth, ls='--', color='black', alpha=0.6, label='$l_{smth}$')
        axis.axvline(mpolmerg, ls='-.', color='black', alpha=0.6, label='$l_{merg}$')
        
        axis.set_ylabel('$l(2l+1)C_l/4\pi$')
        axis.set_xlabel('$l$')
        axis.set_ylim([1e-3, 1.])
        axis.set_xlim([amin(mpol), amax(mpol)])
        axis.legend(loc=3, ncol=2)
        plt.tight_layout()
        path = gdat.pathimag + 'psec%04d.pdf' % i
        plt.savefig(path)
        plt.close(figr)

        for plotigal in [False, True]:

            if plotigal:
                minmlgal = -20.
                maxmlgal = 20.
                minmbgal = -20.
                maxmbgal = 20.
            else:
                minmlgal = -180.
                maxmlgal = 180.
                minmbgal = -90.
                maxmbgal = 90.
               
            path = gdat.pathimag + 'mapsfdfm%04d.pdf' % i
            tdpy.util.plot_maps(path, mapsfdfm[i, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            
            if i == 0:
                path = gdat.pathimag + 'mapsplnk.pdf'
                tdpy.util.plot_maps(path, mapsplnk, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            
            path = gdat.pathimag + 'mapsmerg%04d.pdf' % i
            tdpy.util.plot_maps(path, mapsmerg[i, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            
            path = gdat.pathimag + 'mapsresifdfm%04d.pdf' % i
            tdpy.util.plot_maps(path, mapsmerg[i, :] - mapsfdfm[i, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, resi=True, satu=True)
            
            path = gdat.pathimag + 'mapsresiplnk%04d.pdf' % i
            tdpy.util.plot_maps(path, mapsmerg[i, :] - mapsplnk, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, resi=True, satu=True)
        

def writ_ferm_back():

    gdat = tdpy.util.gdatstrt()
    
    writ=True
    
    strgexpr = 'sbrtfermrec8pntsigal0256.fits'
    
    maxmgangdata = 20.
    
    #listnameback = ['fdfm', 'dustsfdd', 'cmon', 'hydr', 'hasl', 'haslwise', 'wise', 'dark', 'pnlkdust', 'wisestar']
    listnameback = ['fdfm', 'dustsfdd', 'cmon', 'hydr', 'wise', 'dark']
    #listnameback = ['fdfm', 'dark']
    gdat.numbback = len(listnameback)
    gdat.indxback = arange(gdat.numbback)

    # construct the global object
    minmlgal = -maxmgangdata
    maxmlgal = maxmgangdata
    minmbgal = -maxmgangdata
    maxmbgal = maxmgangdata
    
    enertype = 'back'

    # axes
    if enertype == 'back':
        gdat.numbener = 30
        gdat.minmener = 0.3
        gdat.maxmener = 300.
    if enertype == 'pnts':
        gdat.numbener = 3
        gdat.minmener = 0.3
        gdat.maxmener = 10.
    gdat.binsener, gdat.meanener, gdat.diffener, gdat.numbener, gdat.indxener = tdpy.util.retr_axis(minm=gdat.minmener, maxm=gdat.maxmener, numb=gdat.numbener, scal='logt')
    
    gdat.strgbinsener = ['%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]) for i in gdat.indxener]
    
    ## pixelization
    gdat.numbside = 256
    gdat.lgalheal, gdat.bgalheal, gdat.numbpixl, gdat.apix = tdpy.util.retr_healgrid(gdat.numbside)
    gdat.indxpixlrofi = where((abs(gdat.lgalheal) < maxmgangdata) & (abs(gdat.bgalheal) < maxmgangdata))[0]
    gdat.indxpixlnorm = where((abs(gdat.lgalheal) < 10.) & (abs(gdat.bgalheal) < 10.))[0]
    gdat.indxpixl = arange(gdat.numbpixl)
       
    # setup
    gdat.rtag = ''

    # paths
    gdat.pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'ferm_igal/', 'ferm_igal/', gdat.rtag)
     
    ## data
    #path = gdat.pathdata + strgexpr
    #gdat.sbrtdata = pf.getdata(path)

    # power spectrum calculation
    gdat.numbmapsplot = gdat.numbback + 1
    gdat.mapsplot = empty((gdat.numbmapsplot, gdat.numbpixl))
    
    gdat.sbrtisot = tdpy.util.retr_isot(gdat.binsener)
    
    smth = True
    
    listrecotype = ['rec8']
    #listrecotype = ['manu', 'rec7', 'rec8']
    for recotype in listrecotype:
        gdat.evtt, gdat.numbevtt, gdat.indxevtt = tdpy.util.retr_evttferm(recotype)
        gdat.numbevtt = 2
        gdat.evtt = [16, 32]
        gdat.indxevtt = [0, 1]
        ## templates
        gdat.sbrtback = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        gdat.sbrtbacksmth = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        gdat.sbrtbacknorm = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for c, strg in enumerate(listnameback):

            # temp -- ROI should be fixed at 40 X 40 degree^2
            path = gdat.pathdatapcat + strg + '.fits'
            if False and os.path.isfile(path) and not writ:
                print 'Reading %s...' % path
                gdat.sbrtback[c, :, :, :] = pf.getdata(path)
            else:
                
                print 'listnameback[c]'
                print listnameback[c]
                print
                
                if strg.startswith('fdfm'):
                    sbrtbacktemp = tdpy.util.retr_sbrtfdfm(gdat.binsener) 
                elif strg == 'dustsfdd':
                    pathtemp = gdat.pathdata + 'lambda_fds_dust_94GHz.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = pf.getdata(pathtemp, 1)['TEMPERATURE']
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
                elif strg == 'dustplnk':
                    pathtemp = gdat.pathdata + 'plnk/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = pf.getdata(pathtemp, 1)['RADIANCE']
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
                elif strg == 'cmon':
                    pathtemp = gdat.pathdata + 'lambda_wco_dht2001.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = pf.getdata(pathtemp, 1)['TEMPERATURE']
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
                elif strg == 'hydr':
                    pathtemp = gdat.pathdata + 'lambda_combined_nh.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = pf.getdata(pathtemp, 1)['TEMPERATURE']
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
                elif strg == 'hasl':
                    pathtemp = gdat.pathdata + 'haslam408_dsds_Remazeilles2014.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = pf.getdata(pathtemp)['TEMPERATURE']
                    print 'sbrtbacktemp'
                    summgene(sbrtbacktemp)
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
                elif strg == 'haslwise':
                    pathtemp = gdat.pathdata + 'lambda_sfd_ebv.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='NESTED', order_out='RING')
                elif strg == 'wise':
                    pathtemp = gdat.pathdata + 'wssa_sample_1024.fits'
                    tdpy.util.read_fits(pathtemp, verbtype=2) 
                    sbrtbacktemp = pf.getdata(pathtemp, 0)
                    sbrtbacktemp = hp.ud_grade(sbrtbacktemp, gdat.numbside, order_in='RING', order_out='RING')
                elif strg == 'dark':
                    sbrtbacktemp = tdpy.util.retr_nfwp(1., gdat.numbside)
                else:
                    raise Exception('')
                # temp
                #gdat.sbrtback[where(gdat.sbrtback < 0.)] = 0.
                
                # make copies
                for m in gdat.indxevtt:
                    if listnameback[c].startswith('fdfm'):
                        gdat.sbrtback[c, :, :, m] = sbrtbacktemp
                    else:
                        for i in gdat.indxener:
                            gdat.sbrtback[c, i, :, m] = sbrtbacktemp

                print 'gdat.sbrtback[c, :, :, :]'
                summgene(gdat.sbrtback[c, :, :, :])
                
                # normalize
                for i in gdat.indxener:
                    for m in gdat.indxevtt: 
                        gdat.sbrtbacknorm[c, i, :, m] = gdat.sbrtback[c, i, :, m] / mean(gdat.sbrtback[c, i, gdat.indxpixlnorm, m])
                print 'gdat.sbrtbacknorm[c, :, :, :]'
                summgene(gdat.sbrtbacknorm[c, :, :, :])
                path = gdat.pathdatapcat + 'sbrt' + strg + enertype + '.fits'
                print 'Writing to %s...' % path
                pf.writeto(path, gdat.sbrtbacknorm[c, :, :, :], clobber=True)
                
                if smth:
                    print 'gdat.sbrtback'
                    summgene(gdat.sbrtback)
                    print
                    gdat.sbrtbacksmth[c, :, :, :] = tdpy.util.smth_ferm(gdat.sbrtback[c, :, :, :], gdat.meanener, recotype, kerntype='gaus')
                    
                    indxbadd = where(gdat.sbrtbacksmth[c, :, :, :] < 0.)
                    if indxbadd[0].size > 0:
                        print 'Smoothed template went negative. Cutting at 0.'
                        gdat.sbrtbacksmth[c, :, :, :][indxbadd] = 0.

                    print 'gdat.sbrtbacksmth[c, :, :, :]'
                    summgene(gdat.sbrtbacksmth[c, :, :, :])
                    path = gdat.pathdatapcat + 'sbrt' + strg + enertype + 'smth%s.fits' % recotype
                    #print 'Writing to %s...' % path
                    #pf.writeto(path, gdat.sbrtbacksmth[c, :, :, :], clobber=True)
                    
                    # normalize
                    for i in gdat.indxener:
                        for m in gdat.indxevtt: 
                            gdat.sbrtbacksmth[c, i, :, m] = gdat.sbrtbacksmth[c, i, :, m] / mean(gdat.sbrtbacksmth[c, i, gdat.indxpixlnorm, m])
                    
                    print 'gdat.sbrtbacksmth[c, :, :, :]'
                    summgene(gdat.sbrtbacksmth[c, :, :, :])
                    path = gdat.pathdatapcat + 'sbrt' + strg + enertype + 'smth%s.fits' % recotype
                    print 'Writing to %s...' % path
                    pf.writeto(path, gdat.sbrtbacksmth[c, :, :, :], clobber=True)
            
    #merg_maps(numbside=512)
    #merg_maps(mpolmerg=360.)
    #merg_maps(mpolmerg=90.)
    
    # load the map to the array whose power spectrum will be calculated
    gdat.mapsplot[1:, :] = gdat.sbrtback[:, 0, :, 0]
    
    # plot the power spectra
    listlabl = ['Data', 'Isotropic']
    listcolr = ['black', 'b']
    listlabl.append('FDM, %s' % gdat.strgbinsener[1])
    listcolr.append('g')
    listlabl.extend(['Planck', r'WISE 12$\mu$m', 'NFW'])
    listcolr.extend(['r', 'm', 'y'])

    figr, axis = plt.subplots()
    mpol = arange(3 * gdat.numbside, dtype=int)
    for n in range(gdat.numbback):
        psec = hp.anafast(gdat.mapsplot[n, :])
        axis.loglog(mpol, mpol * (mpol + 1.) * psec, color=listcolr[n], label=listlabl[n])
    axis.set_ylabel('$l(l+1)C_l$')
    axis.set_xlabel('$l$')
    axis.legend(loc=4, ncol=2)
    path = gdat.pathimag + 'psec.pdf'
    plt.tight_layout()
    plt.savefig(path)

    # plot the input spatial templates
    for c in gdat.indxback:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                path = gdat.pathimag + 'sbrtback_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_maps(path, gdat.sbrtback[c, i, :, m], indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixl, \
                                                                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
           
    # plot the spectra of spatially averaged background components
    
    listlabl = ['Data', 'Isotropic', 'FDM', 'PlanckDust', r'WISE 12$\mu$m', 'SFD', 'NFW']

    figr, axis = plt.subplots()
    
    numbvarb = gdat.numbback + 1
    listydat = empty((numbvarb, gdat.numbener))
    listyerr = zeros((2, numbvarb, gdat.numbener))
    
    xdat = gdat.meanener
    for k in range(numbvarb):
        ydat = gdat.meanener**2 * listydat[k, :]
        yerr = gdat.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5, label=listlabl[k])

    # Fermi-LAT results
    
    listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
    listmrkr = ['o', 's', 'p', '*', 'D', '^']
    listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
    listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
    for k, name in enumerate(listname):
        path = os.environ["TDGU_DATA_PATH"] + '/ferm_igal/data/fermspec' + name + '.csv'
        data = loadtxt(path)
        enertemp = data[:, 0] # [GeV]
        fluxtemp = data[:, 1] * 1e-3 # [GeV/cm^2/s/sr]
        fluxtemp = interp(gdat.meanener, enertemp, fluxtemp)
        #fluxtemp = interpolate.interp1d(enertemp, fluxtemp)(gdat.meanener)
        axis.plot(gdat.meanener, fluxtemp, marker=listmrkr[k], color=listcolr[k], label=listlabl[k])

    axis.set_xlim([amin(gdat.binsener), amax(gdat.binsener)])
    axis.set_yscale('log')
    axis.set_xlabel('E [GeV]')
    axis.set_xscale('log')
    axis.set_ylabel('$E^2dN/dAdtd\Omega dE$ [GeV/cm$^2$/s/sr]')
    axis.legend(loc=4, ncol=2)

    path = gdat.pathimag + 'backspec.pdf'
    plt.tight_layout()
    plt.savefig(path)
    plt.close(figr)

    # temp
    #indxpixlmean = where(abs(gdat.bgalheal) < 2.)[0]
    #plot_backspec(gdat, indxpixlmean)
    #indxpixlmean = where((abs(gdat.lgalheal) < 5.) & (abs(gdat.bgalheal) < 5.))[0]
    #plot_backspec(gdat, indxpixlmean)


def test_ferm_bubb(strgcnfgextnexec=None):
    
    anglfact = 180. / pi

    dictargs = {}
    dictargs['truenumbelempop0reg0'] = 0
    dictargs['maxmnumbelempop0reg0'] = 0
    
    dictargs['anlytype'] = 'rec8back'
    dictargs['minmlgalrofi'] = -25. / anglfact
    dictargs['maxmlgalrofi'] = 25. / anglfact
    dictargs['minmbgalrofi'] = -25. / anglfact
    dictargs['strgexpo'] = 'expofermrec8backigal0256.fits'
    dictargs['maxmbgalrofi'] = 25. / anglfact
    
    dictargs['psfnevaltype'] = 'init'
    dictargs['trueelemregitype'] = [True, True, True]
    
    dictargs['backtype'] = [[1., 'sbrtfdfmbacksmthrec8.fits', 'sbrtdarkbacksmthrec8.fits']]
    
    listnamecnfgextn = ['nomi']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['nomi']['forccart'] = False
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_ferm_mock_brek(strgcnfgextnexec=None):
     
    dictargs = {}
    dictargs['truemaxmnumbelempop0reg0'] = 100
    
    dictargs['forccart'] = True
    dictargs['pixltype'] = 'cart'
    dictargs['fluxdisttype'] = ['dpowslopbrek']
    dictargs['elemtype'] = ['lghtpnts']
    dictargs['numbsidecart'] = 100
    dictargs['truefluxdistslopupprpop0'] = 1.5
    dictargs['truefluxdistsloplowrpop0'] = 2.
    dictargs['truefluxdistbrekpop0'] = 3e-8
    #dictargs['listnameback'] = ['isot']
    dictargs['backtype'] = [[1e1]]
    
    # temp
    #dictargs['mockonly'] = True
    
    listnamecnfgextn = ['brekloww', 'brekhigh', 'nomi', 'brekstep', 'brekflat', 'brekdown']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['brekloww']['truefluxdistbrekpop0'] = 1e-8
    
    dictargsvari['brekhigh']['truefluxdistbrekpop0'] = 1e-7
    
    dictargsvari['brekstep']['truefluxdistsloplowrpop0'] = 2.5
    
    dictargsvari['brekflat']['truefluxdistsloplowrpop0'] = 1.5
    
    dictargsvari['brekdown']['truefluxdistsloplowrpop0'] = 1.1
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_ferm_igal_mock_flat(strgcnfgextnexec=None):
    
    anglfact = 180. / pi
    
    dictargs = {}
    dictargs['truemaxmnumbelempop0reg0'] = 400
    dictargs['truenumbelempop0reg0'] = 400
    
    dictargs['listnameback'] = ['isot']
    dictargs['backtype'] = [[10.]]
    dictargs['truenumbpopl'] = 1
    dictargs['refrlegdpopl'] = ['PS']
    dictargs['trueelemtype'] = ['lghtpnts']
    dictargs['maxmgangdata'] = 10. / anglfact
    dictargs['truespatdisttype'] = ['self']
    dictargs['spectype'] = ['powr']
    dictargs['psfnevaltype'] = 'kern'
    dictargs['trueelemregitype'] = [True]
    dictargs['proppsfp'] = False
    
    dictargs['fittnumbpopl'] = 1
    dictargs['fittelemtype'] = ['lghtpnts']
    dictargs['fittspatdisttype'] = ['self']
    #dictargs['fittspectype'] = ['colr']
    dictargs['fittmaxmnumbelempop0reg0'] = 1000
    
    #dictargs['strgexprsbrt'] = 'sbrtfermrec8pntsigal0256.fits'
    #dictargs['spectype'] = ['colr']
    #dictargs['listnameback'] = ['isot', 'fdfm', 'dark']
    #dictargs['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'sbrtdarkpntssmthrec8.fits']]
    #dictargs['psfnevaltype'] = 'kern'
    
    dictargs['forccart'] = True
    dictargs['pixltype'] = 'cart'
    dictargs['numbsidecart'] = 100
    
    #dictargs['forccart'] = True
    #dictargs['pixltype'] = 'cart'
    #dictargs['numbsidecart'] = 100
    
    dictargs['numbswep'] = 1000
    dictargs['inittype'] = 'refr'
    #dictargs['numbsamp'] = 500000
    dictargs['numbsamp'] = 10
    
    listnamecnfgextn = ['nomi', 'parsloww', 'parsnone', \
                       ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    #dictargsvari['nomi']['checprio'] = True
    
    dictargsvari['parsloww']['priofactdoff'] = 0.5

    dictargsvari['parsnone']['priofactdoff'] = 0.
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_ferm_igal_mock(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['truemaxmnumbelempop0reg0'] = 100
    dictargs['truemaxmnumbelempop1reg0'] = 100
    dictargs['truemaxmnumbelempop2reg0'] = 100
    dictargs['truenumbelempop0reg0'] = 100
    dictargs['truenumbelempop1reg0'] = 100
    dictargs['truenumbelempop2reg0'] = 100
    
    dictargs['truenumbpopl'] = 3
    dictargs['refrlegdpopl'] = ['AGN', 'Disk MSP', 'GC MSP']
    dictargs['trueelemtype'] = ['lghtpnts', 'lghtpnts', 'lghtpntspuls']
    dictargs['truefluxdistsloppop0'] = 2.6
    dictargs['truefluxdistsloppop1'] = 2.6
    dictargs['truefluxdistsloppop2'] = 3.5
    dictargs['truesinddistmeanpop0'] = 2.
    dictargs['truesinddistmeanpop1'] = 2.
    dictargs['truesinddistmeanpop2'] = 2.
    dictargs['truesinddiststdvpop0'] = 0.5
    dictargs['truesinddiststdvpop1'] = 0.5
    dictargs['truesinddiststdvpop2'] = 0.5
    dictargs['truespatdisttype'] = ['self', 'disc', 'glc3']
    dictargs['truespectype'] = ['powr', 'expc', 'expc']
    dictargs['psfnevaltype'] = 'kern'
    dictargs['trueelemregitype'] = [True, True, True]
    
    dictargs['proppsfp'] = False
    
    dictargs['fittnumbpopl'] = 1
    dictargs['fittelemtype'] = ['lghtpnts']
    dictargs['fittspatdisttype'] = ['self']
    dictargs['fittspectype'] = ['colr']
    dictargs['fittmaxmnumbelempop0reg0'] = 1000
    
    dictargs['forccart'] = True
    dictargs['pixltype'] = 'cart'
    dictargs['numbsidecart'] = 100
    
    #dictargs['numbswep'] = 1000000
    #dictargs['inittype'] = 'refr'
    #dictargs['numbsamp'] = 100
    
    listnamecnfgextn = ['nomi', 'parsloww', 'parsnone', 'truedark', \
                        #'backwfou', 'backfour', 'popl', 'penalpridiff', \
                        ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['nomi']['checprio'] = True
    
    dictargsvari['parsloww']['priofactdoff'] = 0.5

    dictargsvari['parsnone']['priofactdoff'] = 0.
    
    dictargsvari['truedark']['listnameback'] = ['isot', 'fdfm', 'dark']
    dictargsvari['truedark']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'sbrtdarkpntssmthrec8.fits']]
    dictargsvari['truedark']['truemaxmnumbelempop0reg0'] = 0
    dictargsvari['truedark']['truemaxmnumbelempop1reg0'] = 0
    dictargsvari['truedark']['truemaxmnumbelempop2reg0'] = 0
    dictargsvari['truedark']['truenumbelempop0reg0'] = 0
    dictargsvari['truedark']['truenumbelempop1reg0'] = 0
    dictargsvari['truedark']['truenumbelempop2reg0'] = 0
   
    #dictargsvari['backwfou']['fittbacktype'] = [['bfunwfou0004']]
    
    #dictargsvari['backfour']['fittbacktype'] = [['bfunfour0004']]
    
    #dictargsvari['popl']['fittnumbpopl'] = 3
    #dictargsvari['popl']['fittspatdisttype'] = ['self', 'disc', 'los3']
    #dictargsvari['popl']['fittspectype'] = ['powr', 'expc', 'expc']
    #dictargsvari['popl']['fittelemregitype'] = [True, True, True]
    #dictargsvari['popl']['fittmaxmnumbelempop1reg0'] = 1000
    #dictargsvari['popl']['fittmaxmnumbelempop2reg0'] = 1000
    
    #dictargsvari['penalpridiff']['penalpridiff'] = True
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                  forcprev=True, \
                                 )


def pcat_ferm_igal_inpt(strgcnfgextnexec=None):
   
    initpsfp = array([ \
                0.16186914,  2.56138951,  0.66329441,  4.99256131,  0.36684878, \
                0.17057817,  2.57968511,  0.85638406,  3.13328478,  0.61303915, \
                0.13568948,  2.42704043,  0.60715504,  3.94330250,  0.64194487, \
                0.20564394,  2.60527561,  0.73475433,  3.14712137,  0.49128431, \
                0.24889309,  5.00905942,  0.94812157,  2.59796587,  0.59966941, \
                0.25425057,  5.39808726,  0.91544108,  3.86351548,  0.63255054, \
               ])
    anglfact = 180. / pi

    # overplot NPTF results

    dictargs = {}
    dictargs['strgexprsbrt'] = 'sbrtfermrec8pntsigal0256.fits'
    dictargs['spectype'] = ['colr']
    dictargs['savestat'] = True
    dictargs['listnameback'] = ['isot', 'fdfm', 'dark']
    dictargs['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'sbrtdarkpntssmthrec8.fits']]
    
    dictargs['listprefhistfluxlabl'] = ['NPTF']
    dictargs['listprefhistfluxflux'] = [[1e-11, 1e-10, 3e-10, 5e-10]]
    dictargs['listprefhistfluxhist'] = [array([1e6, 1e7, 1e4, 1e6])]
    dictargs['listprefhistfluxtype'] = ['line']
   
    dictargs['initpsfp'] = initpsfp
    dictargs['inittype'] = 'reco'
    
    dictargs['fittmaxmnumbelempop0reg0'] = 1000
    
    #dictargs['propbacp'] = False
    #dictargs['propcomp'] = False
    #dictargs['probtran'] = 0.
    dictargs['proppsfp'] = False
    
    #dictargs['initpsfprefr'] = True
    dictargs['maxmgangdata'] = 10. / anglfact
    dictargs['psfnevaltype'] = 'kern'
    
    dictargs['forccart'] = True
    dictargs['pixltype'] = 'cart'
    dictargs['numbsidecart'] = 100
    
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 10
    
    #dictargs['probtran'] = 0.
    #dictargs['propcomp'] = False
    
    #dictargs['forcsavestat'] = True
    
    listnamecnfgextn = ['chec', 'nomi', 'parsloww', 'parslowr', 'parsvlow', 'parsnone', 'darknone', 'darknoneparsnone', 'mask', 'maskparsnone', \
                        #'darknone', 'exce', 'excefixd', 'rofilarg', 'mask', 'isotfixd', 'heal', 'rec7', \
                        #'backtemp', 'psfnvari', 'backsmth', \
                        #'backwfou', 'backfour', 'popl', \
                        ]
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['chec']['checprio'] = True
                   
    dictargsvari['parsloww']['priofactdoff'] = 0.75
    
    dictargsvari['parslowr']['priofactdoff'] = 0.5
    
    dictargsvari['parsvlow']['priofactdoff'] = 0.25
    
    dictargsvari['parsnone']['priofactdoff'] = 0.
    
    dictargsvari['darknone']['listnameback'] = ['isot', 'fdfm']
    dictargsvari['darknone']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits']]
    
    dictargsvari['darknoneparsnone']['listnameback'] = ['isot', 'fdfm']
    dictargsvari['darknoneparsnone']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits']]
    dictargsvari['darknoneparsnone']['priofactdoff'] = 0.
    
    dictargsvari['mask']['mask'] = array([-30., 30., -2., 2.]) / 180. * pi
    dictargsvari['mask']['listmask'] = [['hstr', -2. / anglfact, 2. / anglfact]]
    
    dictargsvari['maskparsnone']['mask'] = array([-30., 30., -2., 2.]) / 180. * pi
    dictargsvari['maskparsnone']['listmask'] = [['hstr', -2. / anglfact, 2. / anglfact]]
    dictargsvari['maskparsnone']['priofactdoff'] = 0.
    
    #dictargsvari['exce']['numbswep'] = 10000
    #dictargsvari['exce']['numbsamp'] = 100
    #dictargsvari['exce']['inittype'] = 'refr'
    #dictargsvari['exce']['listnameback'] = ['isot', 'fdfm']
    #dictargsvari['exce']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits']]
    #
    #dictargsvari['excefixd']['maxmgangdata'] = 5. / anglfact
    #dictargsvari['excefixd']['numbswep'] = 10000
    #dictargsvari['excefixd']['numbsamp'] = 100
    #dictargsvari['excefixd']['probtran'] = 0.
    #dictargsvari['excefixd']['probcomp'] = False
    #dictargsvari['excefixd']['inittype'] = 'reco'
    #dictargsvari['excefixd']['forcsavestat'] = True
    #dictargsvari['excefixd']['listnameback'] = ['isot', 'fdfm']
    #dictargsvari['excefixd']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits']]
    
    #dictargsvari['rofilarg']['maxmgangdata'] = 15. / anglfact
    
    #dictargsvari['isotfixd']['initbacpbac0ene0'] = 4e-6
    #dictargsvari['isotfixd']['propbacpbac0ene0'] = False
    
    #dictargsvari['heal']['pixltype'] = 'heal'
    #dictargsvari['heal']['forccart'] = False
    
    #dictargsvari['rec7']['strgexpo'] = 'expofermrec7pntsigal0256.fits'
                 
    #dictargsvari['backtemp']['listnameback'] = ['isot', 'hydr', 'cmon', 'dustsfdd', 'dark', 'wise']
    #dictargsvari['backtemp']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'hydrpntssmthrec8.fits', 'cmonpntssmthrec8.fits', \
    #                                                                                'dustsfddpntssmthrec8.fits', 'sbrtdarkpntssmthrec8.fits', 'wisepntssmthrec8.fits']]
    #
    #dictargsvari['psfnvari']['proppsfp'] = True
    #
    #dictargsvari['backsmth']['proppsfp'] = True
    #dictargsvari['backsmth']['psfnevaltype'] = 'full'
    #dictargsvari['backsmth']['backtype'] = [[1., 'sbrtfdfmpnts.fits', 'sbrtdarkpnts.fits']]
    
    #dictargsvari['popl']['numbpopl'] = 3
    
    #dictargsvari['backwfou']['listnameback'] = ['isot', 'fdfm', 'bfunwfou0004']
    #dictargsvari['backwfou']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'bfunwfou0004']]
    
    #dictargsvari['backfour']['listnameback'] = ['isot', 'fdfm', 'bfunfour0004']
    #dictargsvari['backfour']['backtype'] = [[1., 'sbrtfdfmpntssmthrec8.fits', 'bfunfour0002']]
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                  
                                  forcprev=True, \
                                 )
    

def test_ferm_inpt_ptch(strgcnfgextnexec=None):

    anglfact = 180. / pi
    pathdata = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    
    lgalcntr = 0.
    bgalcntr = 5. / anglfact
    
    # rotate the input map
    if pixltype == 'heal' and (lgalcntr != 0. or bgalcntr != 0.):
        strgcntr = 'cntr%04d%04d' % (rad2deg(lgalcntr), rad2deg(bgalcntr))
        
        listcubertttr = [sbrtdata, expo]
        for strgbacktemp in strgback:
            if isinstance(strgbacktemp, str):
                listvarb.append(strgbacktemp)
        numbcuberttr = len(listcuberttr)

        for k in range(numbcuberttr):
            path = pathdata + strgback[k] + strgcntr + '.fits'
            if False and os.path.isfile(path):
                print 'Reading %s...' % path
                maps = pf.getdata(path)
            else:
                pathorig = pathdata + strgback[k] + '.fits'
                maps = pf.getdata(pathorig)
                print 'Rotating the data cube in %s...' % pathorig
                numbener = maps.shape[0]
                numbevtt = maps.shape[2]
                numbside = int(sqrt(maps.shape[1] / 12))
                for i in range(numbener):
                    for m in range(numbevtt):
                        almc = hp.map2alm(maps[i, :, m])
                        # temp
                        #hp.rotate_alm(almc, 0., bgalcntr, 0.)
                        print 'i m'
                        print i, m
                        maps[i, :, m] = hp.rotate_alm(almc, lgalcntr, bgalcntr, 0.)
                        #maps[i, :, m] = hp.alm2map(almc, numbside)
                print 'Writing to %s...' % path
                pf.writeto(path, maps, clobber=True)
    
    pcat.main.init( \
              lgalcntr=lgalcntr, \
              bgalcntr=bgalcntr, \
              minmflux=3e-11, \
              maxmflux=3e-6, \
              backtype=[[1., 'sbrtfdfmpntssmth%s%s.fits' % (recotype, strgcntr)]], \
              strgexpo='expofermrec8pnts%s%s.fits' % (recotype, strgcntr), \
              strgexprsbrt='expofermrec8pntsigal0256%s%s.fits' % (recotype, strgcntr), \
             )
    
globals().get(sys.argv[1])(*sys.argv[2:])
