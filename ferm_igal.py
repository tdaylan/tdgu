from __init__ import *

def make_maps_rec7_back():
    
    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec7']
    gdat.enertype = ['back']
    tdpy.util.make_maps_main(gdat, os.environ["FERM_IGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec7', 'back', 'igal', os.environ["FERM_IGAL_DATA_PATH"], 256, 'tim0')


def make_maps_rec8_back():
    
    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec8']
    gdat.enertype = ['back']
    tdpy.util.make_maps_main(gdat, os.environ["FERM_IGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec8', 'back', 'igal', os.environ["FERM_IGAL_DATA_PATH"], 256, 'tim0')


def mergmaps_arry():
    
    mergmaps(numbside=512)
    mergmaps(mpolmerg=360.)
    mergmaps(mpolmerg=90.)


def retr_plnkmapsorig(gdat, strgmapsplnk):

    if strgmapsplnk == 'radi':
        mapsplnkorig = pf.getdata(gdat.pathdata + 'plnk/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 1)['RADIANCE']
        mapsplnkorig = hp.reorder(mapsplnkorig, n2r=True)
    else:
        
        # read sky maps
        if strgmapsplnk[1] == '0':
            strgfrst = 'plnk/LFI_SkyMap_' 
            strgseco = '-BPassCorrected-field-IQU_0256_R2.01_full.fits'
            strgcols = 'TEMPERATURE'
        elif strgmapsplnk[1] == '1' or strgmapsplnk[1] == '2' or strgmapsplnk[1] == '3':
            strgfrst = 'plnk/HFI_SkyMap_'
            strgseco = '-field-IQU_2048_R2.02_full.fits'
            strgcols = 'I_STOKES'
        else:
            strgfrst = 'plnk/HFI_SkyMap_'
            strgseco = '-field-Int_2048_R2.02_full.fits'
            strgcols = 'I_STOKES'
        strg = strgfrst + '%s' % strgmapsplnk[1:] + strgseco
        mapsplnk = pf.getdata(gdat.pathdata + strg, 1)[strgcols]
        mapsplnk = hp.reorder(mapsplnk, n2r=True)

        print 'Changing units...'
        # change units of the sky maps to Jy/sr
        if strgmapsplnk != '0545' and strgmapsplnk != '0857':
            ## from Kcmb
            if calcfactconv:
                # read Planck band transmission data
                if strgmapsplnk[1] == '0':
                    strg = 'LFI_RIMO_R2.50'
                    strgextn = 'BANDPASS_%s' % strgmapsplnk[1:]
                    freqband = pf.open(gdat.pathdata + '/%s.fits' % strg)[strgextn].data['WAVENUMBER'][1:] * 1e9
                else:
                    strg = 'plnk/HFI_RIMO_R2.00'
                    strgextn = 'BANDPASS_F%s' % strgmapsplnk[1:]
                    freqband = 1e2 * velolght * pf.open(gdat.pathdata + '/%s.fits' % strg)[strgextn].data['WAVENUMBER'][1:]
                tranband = pf.open(gdat.pathdata + '/%s.fits' % strg)[strgextn].data['TRANSMISSION'][1:]
                indxfreqbandgood = where(tranband > 1e-6)[0]
                indxfreqbandgood = arange(amin(indxfreqbandgood), amax(indxfreqbandgood) + 1)

                # calculate the unit conversion factor
                freqscal = consplnk * freqband[indxfreqbandgood] / consbolt / tempcmbr
                freqcntr = float(strgmapsplnk) * 1e9
                specdipo = 1e26 * 2. * (consplnk * freqband[indxfreqbandgood]**2 / velolght / tempcmbr)**2 / consbolt / (exp(freqscal) - 1.)**2 * exp(freqscal) # [Jy/sr]
                factconv = trapz(specdipo * tranband[indxfreqbandgood], freqband[indxfreqbandgood]) / \
                                                trapz(freqcntr * tranband[indxfreqbandgood] / freqband[indxfreqbandgood], freqband[indxfreqbandgood]) # [Jy/sr/Kcmb]
            else:
                # read the unit conversion factors provided by Planck
                factconv = factconvplnk[k, 1] * 1e6
        else:
            ## from MJy/sr
            factconv = 1e6
        mapsplnk *= factconv

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
            print 'total'
            print fwhmpntsplnk.size
            print 'good'
            print indxpntsgood.size
            print 'amin(stdvpntsplnk)'
            print amin(stdvpntsplnk)
            print 'stdvpntsplnk'
            print stdvpntsplnk
            print 'indxpntsgood'
            print indxpntsgood
            lgalpntsplnk = lgalpntsplnk[indxpntsgood]
            bgalpntsplnk = bgalpntsplnk[indxpntsgood]
            fluxpntsplnk = fluxpntsplnk[indxpntsgood]
            stdvpntsplnk = stdvpntsplnk[indxpntsgood]
            print 'amin(stdvpntsplnk)'
            print amin(stdvpntsplnk)
            print 'stdvpntsplnk'
            print stdvpntsplnk
            
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
            
            print 'lgalpntsplnk'
            print lgalpntsplnk
            print 'bgalpntsplnk'
            print bgalpntsplnk
            print 'fluxpntsplnk'
            print fluxpntsplnk
            print 'stdvpntsplnk'
            print stdvpntsplnk
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


def mergmaps(numbside=256, mpolmerg=180., mpolsmth=360., strgmaps='radi'):

    # construct the global object
    gdat = tdpy.util.gdatstrt()
    
    # get the time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # run tag
    rtag = '%s_%04d_%04d_%04d' % (strgtimestmp, numbside, mpolmerg, mpolsmth)
    
    # paths
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'ferm_igal/', 'ferm_igal/mergmaps/', rtag)

    # time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()

    timeinit = time.time()

    calcfactconv = False
    gdat.subspnts = True
    
    retr_axisener(gdat)
    
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
    mapsfermorig -= mean(mapsfermorig, 1)[:, None]
    mapsfermorig /= std(mapsfermorig, 1)[:, None]
    mapsferm = empty_like(mapsfermorig)
    for i in arange(numbener):
        tdpy.util.plot_maps(gdat.pathimag + 'mapsfermorig%04d.pdf' % i, mapsfermorig[i, :], satu=True)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsferm%04d.pdf' % i, mapsferm[i, :], satu=True)
        mapsferm[i, :] = tdpy.util.smth(mapsfermorig[i, :], mpolsmth, mpol=True)
    numbsideferm = int(sqrt(mapsfermorig.shape[1] / 12))

    # 3FGL flux map
    path = gdat.pathdata + 'gll_psc_v16.fit'
    datafgl3 = pf.getdata(path)
    lgalfgl3 = datafgl3['glon']
    lgalfgl3 = ((lgalfgl3 - 180.) % 360.) - 180.
    bgalfgl3 = datafgl3['glat']
    stdvfgl3 = tdpy.util.retr_fwhmferm() / 2.
    specfgl3 = stack((datafgl3['Flux100_300'], datafgl3['Flux300_1000'], datafgl3['Flux1000_3000'], \
                                                            datafgl3['Flux3000_10000'], datafgl3['Flux10000_100000'])) / gdat.diffenerfull[:, None]
    mapspntsferm = tdpy.util.retr_mapspnts(lgalfgl3, bgalfgl3, stdvfgl3, specfgl3, verbtype=2, numbside=numbsideferm)

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
        path = gdat.pathimag + 'tran_%s.pdf' % rtag
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
        path = gdat.pathdata + 'plnk/HFI_Mask_PointSrc_2048_R2.00.fits'
        mapsmask = pf.open(path)[1].data['F353']
        mapsmask = hp.reorder(mapsmask, n2r=True)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsmask_%s.pdf' % rtag, mapsmask)
    
    #strgmapsplnk = 'radi'
    strgmapsplnk = '0857'
    
    # get input maps
    ## Planck map
    print 'Smoothing the Planck map...'
    path = gdat.pathdata + 'mapsplnk_%s.fits' % rtag
    if os.path.isfile(path):
        print 'Reading %s...' % path
        mapsplnk = pf.getdata(path)
    else:
        mapsplnkorig = retr_plnkmapsorig(gdat, strgmapsplnk)
        mapsplnkorig -= mean(mapsplnkorig)
        mapsplnkorig /= std(mapsplnkorig)
        tdpy.util.plot_maps(gdat.pathimag + 'mapsplnkorig_%s.pdf' % rtag, mapsplnkorig, satu=True)
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
        

def retr_modlflux(gdat, sampvarb):

    norm = sampvarb.reshape((gdat.numbback, gdat.numbener))
    modlflux = norm[:, :, None, None] * gdat.fluxback
   
    return modlflux


def retr_llik(sampvarb, gdat, gdatintr):

    modlflux = retr_modlflux(gdat, sampvarb)
    modlfluxtotl = sum(modlflux, axis=0)
    modlcnts = modlfluxtotl * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    llik = sum(gdat.datacnts * log(modlcnts) - modlcnts)
    sampcalc = []

    return llik, sampcalc


def retr_datapara(gdat):
    
    gdat.numbpara = gdat.numbener * gdat.numbback
    gdat.indxpara = arange(gdat.numbpara)

    datapara = tdpy.util.gdatstrt()

    datapara.indx = dict()
    datapara.minm = zeros(gdat.numbpara)
    datapara.maxm = zeros(gdat.numbpara)
    datapara.name = empty(gdat.numbpara, dtype=object)
    datapara.scal = empty(gdat.numbpara, dtype=object)
    datapara.labl = empty(gdat.numbpara, dtype=object)
    datapara.unit = empty(gdat.numbpara, dtype=object)
    datapara.vari = zeros(gdat.numbpara)
    datapara.true = zeros(gdat.numbpara)

    for n in gdat.indxpara:
        datapara.indx['norm%04d' % n] = n
        datapara.name[n] = 'norm%04d' % n
        datapara.scal[n] = 'logt'
        if n // gdat.numbener == 0 or n // gdat.numbener == 1:
            datapara.minm[n] = 1e-4
            datapara.maxm[n] = 1e1
        else:
            datapara.minm[n] = 1e-9
            datapara.maxm[n] = 1e-3
        if n // gdat.numbener == 0:
            strg = 'isot'
        if n // gdat.numbener == 1:
            strg = 'fdfm'
        if n // gdat.numbener == 2:
            strg = 'plnk'
        if n // gdat.numbener == 3:
            strg = 'wise'
        if n // gdat.numbener == 5:
            strg = 'fink'
        if n // gdat.numbener == 4:
            strg = 'dark'
        datapara.labl[n] = '$A_{%d}^{%s}$' % (n % gdat.numbener, strg)
        datapara.unit[n] = ''
        datapara.vari[n] = 1e-1
        datapara.true[n] = None
    
    datapara.strg = datapara.labl + ' ' + datapara.unit
    
    return datapara


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


def retr_axisener(gdat):
    
    gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
    gdat.binsenerfull, gdat.meanenerfull, gdat.diffenerfull, gdat.numbenerfull, gdat.indxenerfull = tdpy.util.retr_axis(bins=gdat.binsenerfull, scal='logt')


def regrback( \
              numbproc=1, \
              numbswep=None, \
              datatype='inpt', \
              verbtype=1, \
              makeplot=False, \
              strgexpr='fermflux_cmp0_igal.fits', \
              strgexpo='fermexpo_cmp0_igal.fits', \
              #strgback=['isotflux', 'fdfmflux', 'plnk/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 'wssa_sample_1024.fits', 'lambda_sfd_ebv.fits', 'darktemp'], \
              strgback=['isotflux', 'fdfmflux'], \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(3, 4), \
              maxmgang=20.
             ):

    smthmaps = False
    optiprop = True

    # construct the global object
    gdat = tdpy.util.gdatstrt()
    gdat.numbproc = numbproc
    gdat.numbswep = numbswep
    gdat.datatype = datatype
    gdat.verbtype = verbtype
    gdat.makeplot = makeplot
    gdat.strgexpr = strgexpr
    gdat.strgexpo = strgexpo
    gdat.strgback = strgback
    gdat.indxenerincl = indxenerincl
    gdat.indxevttincl = indxevttincl
    gdat.maxmgang = maxmgang
    
    #theory of everything in a shell
    #ecedmmkj = tanshuuydiomw if hkshdud => hsgdsiiiis 
    #return hsudghin <= TRUE
    #for d(f) if ede = FALSE
    #return tansu => dghsishjklnbjks,  a =  nonsense(tansu=cooncon)


    # axes
    retr_axisener(gdat)
    gdat.binsener = gdat.binsenerfull[gdat.indxenerincl[0]:gdat.indxenerincl[-1] + 2]

    gdat.binsener, gdat.meanener, gdat.diffener, gdat.numbener, gdat.indxener = tdpy.util.retr_axis(bins=gdat.binsener, scal='logt')
    gdat.strgbinsener = ['%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]) for i in gdat.indxener]

    ## event type
    gdat.indxevttfull = arange(4)
    gdat.numbevttfull = gdat.indxevttfull.size
    gdat.indxevttrofi = gdat.indxevttfull[gdat.indxevttincl]
    gdat.numbevtt = gdat.indxevttrofi.size
    gdat.indxevtt = arange(gdat.numbevtt)

    gdat.numbback = len(strgback)
    gdat.indxback = arange(gdat.numbback)

    ## pixelization
    gdat.numbside = 256
    gdat.numbpixlfull = gdat.numbside**2 * 12
    gdat.lgalheal, gdat.bgalheal, gdat.numbpixl, gdat.apix = tdpy.util.retr_healgrid(gdat.numbside)
    gdat.indxpixlrofi = where((abs(gdat.lgalheal) < gdat.maxmgang) & (abs(gdat.bgalheal) < gdat.maxmgang))[0]
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
       
    minmlgal = -gdat.maxmgang
    maxmlgal = gdat.maxmgang 
    minmbgal = -gdat.maxmgang
    maxmbgal = gdat.maxmgang 

    indxdatacubefilt = meshgrid(gdat.indxenerincl, gdat.indxpixlrofi, gdat.indxevttincl, indexing='ij')
    
    # get data structure
    datapara = retr_datapara(gdat)
    
    if gdat.numbswep == None:
        gdat.numbswep = 4 * gdat.numbpara * 1000
    
    # get the time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # setup
    gdat.rtag = '%s_%s' % (strgtimestmp, gdat.datatype)
    
    # paths
    gdat.pathimag, gdat.pathdata = tdpy.util.retr_path('tdgu', 'ferm_igal/', 'ferm_igal/regrback/', gdat.rtag)

    ## data
    if gdat.datatype == 'inpt':
        
        path = gdat.pathdata + gdat.strgexpr
        gdat.exprflux = pf.getdata(path)
        
        ### filter
        gdat.dataflux = gdat.exprflux[indxdatacubefilt]

    ## exposure
    path = gdat.pathdata + gdat.strgexpo
    gdat.expo = pf.getdata(path)
    gdat.expo = gdat.expo[indxdatacubefilt]
   
    gdat.datacnts = gdat.dataflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]

    # power spectrum calculation
    gdat.numbmapsplot = gdat.numbback + 1
    gdat.mapsplot = empty((gdat.numbmapsplot, gdat.numbpixlfull))
    
    gdat.mapsplot[0, gdat.indxpixlrofi] = sum(gdat.datacnts[1, :, :], 1)

    ## templates
    gdat.fluxbackfull = empty((gdat.numbback, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
    for c in gdat.indxback:

        if c == 0:
            strg = 'isotflux'
        if c == 1:
            strg = 'fdfmflux'
        if c == 2:
            strg = 'plnkdust'
        if c == 3:
            strg = 'wisestar'
        if c == 4:
            strg = 'finkdust'
        if c == 5:
            strg = 'darktemp'
        if smthmaps:
            strg += 'smth'

        # temp -- ROI should be fixed at 40 X 40 degree^2
        path = gdat.pathdata + strg + '.fits'
        if os.path.isfile(path):
            print 'Reading %s...' % path
            gdat.fluxbackfull[c, :, :, :] = pf.getdata(path)
        else:
            
            if c == 0:
                gdat.fluxbackorig = tdpy.util.retr_isot(gdat.binsenerfull)
            if c == 1:
                gdat.fluxbackorig = tdpy.util.retr_fdfm(gdat.binsenerfull) 
            if c == 2:
                pathtemp = gdat.pathdata + gdat.strgback[c]
                gdat.fluxbackorig = pf.getdata(pathtemp, 1)['RADIANCE']
                gdat.fluxbackorig = hp.ud_grade(gdat.fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 3:
                pathtemp = gdat.pathdata + gdat.strgback[c]
                gdat.fluxbackorig = pf.getdata(pathtemp, 0)
                gdat.fluxbackorig = hp.ud_grade(gdat.fluxbackorig, gdat.numbside, order_in='RING', order_out='RING')
            if c == 4:
                pathtemp = gdat.pathdata + gdat.strgback[c]
                gdat.fluxbackorig = pf.getdata(pathtemp)['TEMPERATURE']
                gdat.fluxbackorig = hp.ud_grade(gdat.fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 5:
                gdat.fluxbackorig = tdpy.util.retr_nfwp(1., gdat.numbside)

            # normalize the templates
            #if c != 0 and c != 1:
            #    gdat.fluxbackorig /= mean(gdat.fluxbackorig[gdat.indxpixlrofi])
            
            # make copies of the maps
            for m in gdat.indxevttfull:
                if c == 0 or c == 1:
                    gdat.fluxbackfull[c, :, :, m] = gdat.fluxbackorig
                else:
                    for i in gdat.indxenerfull:
                        gdat.fluxbackfull[c, i, :, m] = gdat.fluxbackorig
            
            # smooth the templates
            if smthmaps:
                for c in gdat.indxback:
                    gdat.fluxbackfull[c, :, :, :] = tdpy.util.smth_ferm(gdat.fluxbackfull[c, :, :, :], gdat.meanenerfull, gdat.indxevttfull)
            
            # temp
            #gdat.fluxback[where(gdat.fluxback < 0.)] = 0.

            pf.writeto(path, gdat.fluxbackfull[c, :, :, :], clobber=True)

    # take only the energy bins, spatial pixels and event types of interest
    gdat.fluxback = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.fluxback[c, i, :, m] = gdat.fluxbackfull[c, gdat.indxenerincl[i], gdat.indxpixlrofi, gdat.indxevttincl[m]]

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
                path = gdat.pathimag + 'fluxback_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_maps(path, gdat.fluxback[c, i, :, m], indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                                                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
           
    initsamp = rand(gdat.numbproc * gdat.numbpara).reshape((gdat.numbproc, gdat.numbpara))

    numbplotside = gdat.numbpara
    chan = tdpy.mcmc.init(retr_llik, datapara, numbproc=gdat.numbproc, numbswep=gdat.numbswep, initsamp=initsamp, gdatextr=gdat, optiprop=optiprop, \
                verbtype=gdat.verbtype, pathdata=gdat.pathdata, pathimag=gdat.pathimag, rtag=gdat.rtag, numbplotside=numbplotside)
    
    listsampvarb, listsamp, listsampcalc, listllik, listaccp, listchro, listindxparamodi, propeffi, levi, info, gmrbstat = chan
    numbsamp = listsamp.shape[0]

    gdat.medisampvarb = percentile(listsampvarb, 50., axis=0)
    medimodlflux = retr_modlflux(gdat, gdat.medisampvarb)
    medimodlfluxtotl = sum(medimodlflux, axis=0)
    mediresiflux = gdat.dataflux - medimodlfluxtotl

    gdat.postsampvarb = tdpy.util.retr_postvarb(listsampvarb)

    gdat.postnormback = gdat.postsampvarb.reshape((3, gdat.numbback, gdat.numbener))

    for i in gdat.indxener:
        for m in gdat.indxevtt:
            
            path = gdat.pathimag + 'dataflux_%d%d.pdf' % (i, m)
            tdpy.util.plot_maps(path, gdat.dataflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            for c in gdat.indxback:
                path = gdat.pathimag + 'medimodlflux_%d%d%d.pdf' % (c, i, m)
                tdpy.util.plot_maps(path, medimodlflux[c, i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                                minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = gdat.pathimag + 'medimodlfluxtotl_%d%d.pdf' % (i, m)
            tdpy.util.plot_maps(path, medimodlfluxtotl[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True)
            path = gdat.pathimag + 'mediresiflux_%d%d.pdf' % (i, m)
            tdpy.util.plot_maps(path, mediresiflux[i, :, m] * 1e6, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, \
                            minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, satu=True, resi=True)
    
    # plot the spectra of spatially averaged background components
    
    listlabl = ['Data', 'Isotropic', 'FDM', 'PlanckDust', r'WISE 12$\mu$m', 'SFD', 'NFW']

    figr, axis = plt.subplots()
    
    numbvarb = gdat.numbback + 1
    listydat = empty((numbvarb, gdat.numbener))
    listyerr = zeros((2, numbvarb, gdat.numbener))
    
    listydat[0, :] = mean(sum(gdat.datacnts, 2) / sum(gdat.expo, 2), 1) / gdat.apix / gdat.diffener
    listyerr[:, 0, :] = mean(sqrt(sum(gdat.datacnts, 2)) / sum(gdat.expo, 2), 1) / gdat.apix / gdat.diffener
    for n in gdat.indxback:
        if n == 0 or n == 1:
            listydat[n+1, :] = gdat.postnormback[0, n, :] * mean(mean(gdat.fluxback[n, :, :, :], 1), 1)
        else:
            listydat[n+1, :] = gdat.postnormback[0, n, :]
        listyerr[:, n+1, :] = tdpy.util.retr_errrvarb(gdat.postnormback[:, n, :])
    
    xdat = gdat.meanener
    for k in range(numbvarb):
        ydat = gdat.meanener**2 * listydat[k, :]
        yerr = gdat.meanener**2 * listyerr[:, k, :]
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', markersize=5, label=listlabl[k])

    # Fermi-LAT results
    # temp
    if False:
        if gdat.datatype == 'mock':
            pass
        else:
            if gdat.exprtype == 'ferm':
                listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
                listmrkr = ['o', 's', 'p', '*', 'D', '^']
                listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
                listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
                for k, name in enumerate(listname):
                    path = os.environ["PCAT_DATA_PATH"] + '/fermspec' + name + '.csv'
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


def pcat_ferm_inpt_ptch():

    pathdata = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    lgalcntr = deg2rad(0)
    bgalcntr = deg2rad(45.)
    liststrg = ['fermflux_cmp0_igal', 'fermexpo_cmp0_igal', 'fdfmflux']
    numbmaps = len(liststrg)
    strgcntr = '_cntr%04d%04d' % (rad2deg(lgalcntr), rad2deg(bgalcntr))
    for k in range(numbmaps):
        path = pathdata + liststrg[k] + strgcntr + '.fits'
        if os.path.isfile(path):
            print 'Reading %s...' % path
            maps = pf.getdata(path)
        else:
            pathorig = pathdata + liststrg[k] + '.fits'
            maps = pf.getdata(pathorig)
            numbener = maps.shape[0]
            numbevtt = maps.shape[2]
            numbside = int(sqrt(maps.shape[1] / 12))
            for i in range(numbener):
                for m in range(numbevtt):
                    almc = hp.map2alm(maps[i, :, m])
                    # temp
                    #hp.rotate_alm(almc, 0., bgalcntr, 0.)
                    hp.rotate_alm(almc, lgalcntr, bgalcntr, 0.)
                    maps[i, :, m] = hp.alm2map(almc, numbside)
            pf.writeto(path, maps, clobber=True)
    
    pcat.main.init( \
              numbswep=100000, \
              randinit=True, \
              maxmgang=deg2rad(20.), \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(2, 4), \
              lgalcntr=lgalcntr, \
              bgalcntr=bgalcntr, \
              minmflux=3e-11, \
              maxmflux=3e-6, \
              strgback=['isotflux.fits', 'fdfmflux%s.fits' % strgcntr], \
              strgexpo='fermexpo_cmp0_igal%s.fits' % strgcntr, \
              datatype='inpt', \
              strgexpr='fermflux_cmp0_igal%s.fits' % strgcntr, \
             )
    
    
def pcat_ferm_inpt_igal(strgexpr='fermflux_cmp0_igal.fits', strgexpo='fermexpo_cmp0_igal.fits'):
    
    pcat.main.init( \
              numbswep=10000, \
              randinit=False, \
              maxmgang=deg2rad(20.), \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(2, 4), \
              minmflux=1e-9, \
              maxmflux=3e-6, \
              strgback=['isotflux.fits', 'fdfmflux.fits'], \
              strgexpo=strgexpo, \
              datatype='inpt', \
              strgexpr=strgexpr, \
             )
    
    
def pcat_ferm_mock_igal_brok():
     
    #mockfluxdistbrek = array([1e-10, 3e-10, 1e-9, 3e-9, 1e-8])
    mockfluxdistbrek = array([1e-9])
    listmockfluxdistsloplowr = array([1.9, 2.2, 2.8, 3.1, 3.4])
    numbiter = listmockfluxdistsloplowr.size
    for k in range(numbiter):
    
        pcat.main.init( \
                       numbswep=100, \
                       verbtype=2, \
                       randinit=False, \
                       exprinfo=False, \
                       indxevttincl=arange(2, 4), \
                       indxenerincl=arange(1, 4), \
                       strgexpo='fermexpo_cmp0_igal.fits', \
                       strgback=['isotflux.fits'], \
                       
                       maxmgang=deg2rad(20.), \
                       fluxdisttype=['brok'], \
                       
                       boolpropfluxdistbrek=False, \

                       maxmnumbpnts=array([4]), \
                       minmflux=3e-11, \
                       maxmflux=3e-7, \

                       datatype='mock', \
                       mocknumbpnts=array([4]), \
                       
                       mockfluxdistbrek=mockfluxdistbrek, \
                       mockfluxdistsloplowr=listmockfluxdistsloplowr[k], \
                       mockfluxdistslopuppr=array([1.6]), \
                       
                       mocksinddiststdv=array([.5]), \
                       mocksinddistmean=array([2.]), \
                      )


def pcat_ferm_mock_igal_popl():
     
    pcat.main.init( \
                   numbswep=10000, \
                   randinit=False, \
                   indxevttincl=arange(2, 4), \
                   indxenerincl=arange(1, 4), \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits'], \
                   
                   maxmgang=deg2rad(20.), \
                   fluxdisttype=['brok'], \
                   
                   maxmnumbpnts=array([400]), \
                   minmflux=3e-11, \
                   maxmflux=3e-7, \

                   datatype='mock', \
                   mocknumbpnts=array([400]), \
                   
                   mockspatdisttype=['unif', 'disc', 'gang'], \
                   
                   mocksinddiststdv=array([.5]), \
                   mocksinddistmean=array([2.]), \
                  )


def pcat_ferm_mock_igal():
     
    pcat.main.init( \
                   numbswep=10, \
                   verbtype=2, \
                   randinit=False, \
                   indxevttincl=arange(2, 4), \
                   indxenerincl=arange(1, 4), \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits', 'fdfmflux.fits'], \
                   #maxmnumbpnts=array([2, 2, 4]), \
                   #maxmnumbpnts=array([200, 200, 400]), \
                   maxmnumbpnts=array([2]), \
                   maxmgang=deg2rad(20.), \
                   minmflux=3e-11, \
                   maxmflux=3e-7, \
                   datatype='mock', \
                   #mocknumbpnts=array([1, 1, 2]), \
                   mocknumbpnts=array([3, 3, 3]), \
                   mockspatdisttype=['unif', 'disc', 'gang'], \
                   mockfluxdistslop=array([2.6, 2.6, 3.5]), \
                   mocksinddiststdv=array([.5, .5, .5]), \
                   mocksinddistmean=array([2., 2., 2.]), \
                  )


def regrback_arry():
    
    regrback( \
             verbtype=1, \
             makeplot=True, \
             #strgback=['isotflux', 'fdfmflux'], \
            )


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

