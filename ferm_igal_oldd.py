from __init__ import *

# writing data
def writ_maps_rec7_back():
    
    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec7']
    gdat.enertype = ['back']
    tdpy.util.writ_fdfm()
    tdpy.util.writ_maps_main(gdat, os.environ["FERM_IGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec7', 'back', 'igal', os.environ["FERM_IGAL_DATA_PATH"], 256, 'tim0')


def writ_maps_rec8_back():
    
    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec8']
    gdat.enertype = ['back']
    tdpy.util.writ_maps_main(gdat, os.environ["FERM_IGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec8', 'back', 'igal', os.environ["FERM_IGAL_DATA_PATH"], 256, 'tim0')


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
    
    numbener = 30
    minmener = 0.1
    maxmener = 100.
    binsener = logspace(log10(minmener), log10(maxmener), numbener + 1)
    lowrener = binsener[:-1]
    upprener = binsener[1:]
    limtener = stack((lowrener, upprener), axis=1)
    path = os.environ["TDPY_DATA_PATH"] + '/gtbndefn_back.dat'
    savetxt(path, limtener, fmt='%10.5g')


def merg_maps_arry():
    
    merg_maps(numbside=512)
    merg_maps(mpolmerg=360.)
    merg_maps(mpolmerg=90.)


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
   
    indxevttrofi = arange(4)

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
    indxevttrofi = arange(3, 4)
    numbevtt = indxevttrofi.size
    indxevtt = arange(numbevtt)

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
    path = gdat.pathdata + '/fermflux_cmp0_igal.fits'
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
    stdvfgl3 = tdpy.util.retr_fwhmferm(meanener, indxevttrofi) / 2.
    specfgl3 = stack((datafgl3['Flux100_300'], datafgl3['Flux300_1000'], datafgl3['Flux1000_3000'], \
                                                            datafgl3['Flux3000_10000'], datafgl3['Flux10000_100000'])) / gdat.diffenerfull[:, None]
    
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
    mapsfdfmorig = tdpy.util.retr_fdfm(binsener, numbside=numbside)
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
        

def writ_data():

    verbtype=1
    strgexpr='fermflux_cmp0_igal.fits'

    # construct the global object
    gdat = tdpy.util.gdatstrt()
    gdat.verbtype = verbtype
    gdat.strgexpr = strgexpr
    
    indxevttrofi=arange(3, 4)
    gdat.indxevttrofi = indxevttrofi
    
    maxmgangdata=20.
    gdat.maxmgangdata = maxmgangdata
    minmlgal = -maxmgangdata
    maxmlgal = maxmgangdata
    minmbgal = -maxmgangdata
    maxmbgal = maxmgangdata

    # axes
    gdat.binsener = array([0.1, 0.3, 1., 3., 10., 100.])
    gdat.binsener, gdat.meanener, gdat.diffener, gdat.numbener, gdat.indxener = tdpy.util.retr_axis(bins=gdat.binsener, scal='logt')
    gdat.strgbinsener = ['%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]) for i in gdat.indxener]
    
    ## event type
    gdat.indxevttfull = arange(4)
    gdat.numbevttfull = gdat.indxevttfull.size
    gdat.indxevttrofi = gdat.indxevttfull[gdat.indxevttrofi]
    gdat.numbevtt = gdat.indxevttrofi.size
    gdat.indxevtt = arange(gdat.numbevtt)

    ## pixelization
    gdat.numbside = 256
    gdat.numbpixlfull = gdat.numbside**2 * 12
    gdat.lgalheal, gdat.bgalheal, gdat.numbpixl, gdat.apix = tdpy.util.retr_healgrid(gdat.numbside)
    gdat.indxpixlnorm = where((abs(gdat.lgalheal) < 10.) & (abs(gdat.bgalheal) < 10.))[0]
    gdat.indxpixlrofi = where((abs(gdat.lgalheal) < gdat.maxmgangdata) & (abs(gdat.bgalheal) < gdat.maxmgangdata))[0]
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
       
    indxdatacubefilt = meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevttrofi, indexing='ij')
    
    # paths
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'ferm_igal/', 'ferm_igal/', 'inpt')
    
    ## data
    path = gdat.pathdata + gdat.strgexpr
    gdat.exprflux = pf.getdata(path)
    
    ### filter
    gdat.dataflux = gdat.exprflux[indxdatacubefilt]

    ## templates
    #strgback=['', '', '', 'plnk/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 'wssa_sample_1024.fits', 'lambda_sfd_ebv.fits', '']
    #listnameback = ['isotflux', 'fdfmflux', 'fdfmfluxnorm', 'plnkdust', 'wisestar', 'finkdust', 'darktemp']
    strgback=['']
    listnameback = ['fdfmflux']
    #strgback=['', '', '']
    #listnameback = ['isotflux', 'fdfmflux', 'darktemp']
    for nameback in deepcopy(listnameback):
        listnameback += [nameback + 'smth']
    numbback = len(listnameback)
    gdat.indxback = arange(numbback)
    gdat.fluxbackfull = empty((numbback, gdat.numbener, gdat.numbpixlfull, gdat.numbevttfull))
    
    # power spectrum calculation
    gdat.numbmapsplot = numbback + 1
    gdat.mapsplot = empty((gdat.numbmapsplot, gdat.numbpixlfull))
    gdat.mapsplot[0, gdat.indxpixlrofi] = sum(gdat.exprflux[1, gdat.indxpixlrofi, :], 1)

    for c, strg in enumerate(listnameback):

        # temp -- ROI should be fixed at 40 X 40 degree^2
        path = gdat.pathdata + strg + '.fits'
        if False and os.path.isfile(path):
            print 'Reading %s...' % path
            gdat.fluxbackfull[c, :, :, :] = pf.getdata(path)
        else:
            
            print 'c'
            print c
            print 'strg'
            print strg

            if strg == 'isotflux':
                gdat.fluxbackorig = tdpy.util.retr_isot(gdat.binsener)
            if strg.startswith('fdfmflux'):
                gdat.fluxbackorig = tdpy.util.retr_fdfm(gdat.binsener) 
            if strg == 'plnkdust':
                pathtemp = gdat.pathdata + strgback[c]
                gdat.fluxbackorig = pf.getdata(pathtemp, 1)['RADIANCE']
                gdat.fluxbackorig = hp.ud_grade(gdat.fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if strg == 'wisestar':
                pathtemp = gdat.pathdata + strgback[c]
                gdat.fluxbackorig = pf.getdata(pathtemp, 0)
                gdat.fluxbackorig = hp.ud_grade(gdat.fluxbackorig, gdat.numbside, order_in='RING', order_out='RING')
            if strg == 'finkdust':
                pathtemp = gdat.pathdata + strgback[c]
                gdat.fluxbackorig = pf.getdata(pathtemp)['TEMPERATURE']
                gdat.fluxbackorig = hp.ud_grade(gdat.fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if strg == 'darktemp':
                gdat.fluxbackorig = tdpy.util.retr_nfwp(1., gdat.numbside)

            if strg.endswith('smth'):
                
                indxbackorig = listnameback.index(strg[:-4])
                print 'c'
                print 'indxbackorig'
                print indxbackorig
                print 'gdat.fluxbackfull[indxbackorig, :, :, :]'
                summgene(gdat.fluxbackfull[indxbackorig, :, :, :])
                # smooth
                gdat.fluxbackfull[c, :, :, :] = tdpy.util.smth_ferm(gdat.fluxbackfull[indxbackorig, :, :, :], gdat.meanener, gdat.indxevttfull)
                print 'gdat.fluxbackfull[c, :, :, :]'
                summgene(gdat.fluxbackfull[c, :, :, :])
                # normalize
                for i in gdat.indxener:
                    for m in gdat.indxevttfull:
                        gdat.fluxbackfull[c, i, :, m] = gdat.fluxbackfull[c, i, :, m] / mean(gdat.fluxbackfull[c, i, gdat.indxpixlnorm, m])

            else:
                # make copies
                for m in gdat.indxevttfull:
                    if strg == 'isotflux' or strg.startswith('fdfmflux'):
                        gdat.fluxbackfull[c, :, :, m] = gdat.fluxbackorig
                    else:
                        for i in gdat.indxener:
                            gdat.fluxbackfull[c, i, :, m] = gdat.fluxbackorig
            
            # temp
            #gdat.fluxback[where(gdat.fluxback < 0.)] = 0.
            print 'Writing to %s...' % path
            pf.writeto(path, gdat.fluxbackfull[c, :, :, :], clobber=True)
            print

    # take only the energy bins, spatial pixels and event types of interest
    gdat.fluxback = empty((numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.fluxback[c, i, :, m] = gdat.fluxbackfull[c, gdat.indxener[i], gdat.indxpixlrofi, gdat.indxevttrofi[m]]

    # load the map to the array whose power spectrum will be calculated
    gdat.mapsplot[1:, gdat.indxpixlrofi] = gdat.fluxback[:, 0, :, 0]
    
    # plot the power spectra
    listlabl = ['Data', 'Isotropic']
    listcolr = ['black', 'b']
    listlabl.append('FDM, %s' % gdat.strgbinsener[1])
    listcolr.append('g')
    listlabl.extend(['Planck', r'WISE 12$\mu$m', 'NFW'])
    listcolr.extend(['r', 'm', 'y'])

    figr, axis = plt.subplots()
    mpol = arange(3 * gdat.numbside, dtype=int)
    for n in range(numbback):
        psec = hp.anafast(gdat.mapsplot[n, :])
        axis.loglog(mpol, mpol * (mpol + 1.) * psec)#, label=listlabl[n])
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
                path = gdat.pathimag + 'fluxback_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_maps(path, gdat.fluxback[c, i, :, m], indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                                                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
           
    # plot the spectra of spatially averaged background components
    
    listlabl = ['Data', 'Isotropic', 'FDM', 'PlanckDust', r'WISE 12$\mu$m', 'SFD', 'NFW']

    figr, axis = plt.subplots()
    
    numbvarb = numbback + 1
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


    



globals().get(sys.argv[1])()

