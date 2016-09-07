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


def merg_maps_arry():
    
    merg_maps(numbside=512)
    merg_maps(mpolmerg=360.)
    merg_maps(mpolmerg=90.)


def smth(maps, scalsmth, mpol=False, retrfull=False, numbsideoutp=None, indxpixlmask=None):

    if mpol:
        mpolsmth = scalsmth
    else:
        mpolsmth = 180. / scalsmth

    numbpixl = maps.size
    numbside = int(sqrt(numbpixl / 12))
    numbmpol = 3 * numbside
    maxmmpol = 3. * numbside - 1.
    mpolgrid, temp = hp.Alm.getlm(lmax=maxmmpol)
    mpol = arange(maxmmpol + 1.)
    
    if numbsideoutp == None:
        numbsideoutp = numbside
    
    if indxpixlmask != None:
        mapsoutp = copy(maps)
        mapsoutp[indxpixlmask] = hp.UNSEEN
        mapsoutp = hp.ma(mapsoutp)
        
        almctemp = hp.map2alm(maps)
        print 'almcmasknone'
        for k in range(10):
            print almctemp[k]
        print
    else:
        mapsoutp = maps
    
    almc = hp.map2alm(mapsoutp)
    print 'almcmask'
    for k in range(10):
        print almc[k]
    print

    wght = exp(-0.5 * (mpolgrid / mpolsmth)**2)
    almc *= wght

    mapsoutp = hp.alm2map(almc[where(mpolgrid < 3 * numbsideoutp)], numbsideoutp, verbose=False)

    if retrfull:
        return maps, almc, mpol, exp(-0.5 * (mpol / mpolsmth)**2)
    else:
        return maps 


def merg_maps(numbside=256, mpolmerg=180., mpolsmth=360., strgmaps='radi'):

    timeinit = time.time()

    calcfactconv = False
    subspnts = False
    
    # runtag
    rtag = '%04d%04d%04d' % (numbside, mpolmerg, mpolsmth)
    
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

    ## base path
    pathbase = os.environ["FERM_IGAL_DATA_PATH"]
    pathdata = pathbase + '/mergmaps/'
    pathimag = pathbase + '/imag/mergmaps/'
    pathimagngal = pathbase + '/imag/mergmaps/ngal'
    pathimagigal = pathbase + '/imag/mergmaps/igal'
    os.system('mkdir -p ' + pathimagngal + ' ' + pathimagigal)
    
    # read unit conversion data provided by Planck
    factconvplnk = loadtxt(pathdata + 'plnkunitconv.dat')
    
    if plotfull:
        strgmapsplnk = ['0030', '0044', '0070', '0100', '0143', '0217', '0353', '0545', '0857', 'radi']
    strgmapsplnk = ['radi']
    #strgmapsplnk = ['0857']
    numbmapsplnk = len(strgmapsplnk)
    for k in range(numbmapsplnk):

        print 'Map number ', k
        print 'Maps string: ', strgmapsplnk[k]
        
        if strgmapsplnk[k] == 'radi':
            mapsplnkorig = pf.getdata(pathbase + '/HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 1)['RADIANCE']
            mapsplnkorig = hp.reorder(mapsplnkorig, n2r=True)
        else:
            
            print 'Reading the raw map at %d...' % int(strgmapsplnk[k][1:])
            # read sky maps
            if strgmapsplnk[k][1] == '0':
                strgfrst = '/LFI_SkyMap_' 
                strgseco = '-BPassCorrected-field-IQU_0256_R2.01_full.fits'
                strgcols = 'TEMPERATURE'
            elif strgmapsplnk[k][1] == '1' or strgmapsplnk[k][1] == '2' or strgmapsplnk[k][1] == '3':
                strgfrst = '/HFI_SkyMap_'
                strgseco = '-field-IQU_2048_R2.02_full.fits'
                strgcols = 'I_STOKES'
            else:
                strgfrst = '/HFI_SkyMap_'
                strgseco = '-field-Int_2048_R2.02_full.fits'
                strgcols = 'I_STOKES'
            strg = strgfrst + '%s' % strgmapsplnk[k][1:] + strgseco
            mapsplnk = pf.getdata(pathbase + strg, 1)[strgcols]
            mapsplnk = hp.reorder(mapsplnk, n2r=True)

            print 'Changing units...'
            # change units of the sky maps to Jy/sr
            if strgmapsplnk[k] != '0545' and strgmapsplnk[k] != '0857':
                ## from Kcmb
                if calcfactconv:
                    # read Planck band transmission data
                    if strgmapsplnk[k][1] == '0':
                        strg = 'LFI_RIMO_R2.50'
                        strgextn = 'BANDPASS_%s' % strgmapsplnk[k][1:]
                        freqband = pf.open(pathbase + '/%s.fits' % strg)[strgextn].data['WAVENUMBER'][1:] * 1e9
                    else:
                        strg = 'HFI_RIMO_R2.00'
                        strgextn = 'BANDPASS_F%s' % strgmapsplnk[k][1:]
                        freqband = 1e2 * velolght * pf.open(pathbase + '/%s.fits' % strg)[strgextn].data['WAVENUMBER'][1:]
                    tranband = pf.open(pathbase + '/%s.fits' % strg)[strgextn].data['TRANSMISSION'][1:]
                    indxfreqbandgood = where(tranband > 1e-6)[0]
                    indxfreqbandgood = arange(amin(indxfreqbandgood), amax(indxfreqbandgood) + 1)

                    # calculate the unit conversion factor
                    freqscal = consplnk * freqband[indxfreqbandgood] / consbolt / tempcmbr
                    freqcntr = float(strgmapsplnk[k]) * 1e9
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
            tdpy.util.plot_maps(pathimag + 'mapsplnk%s.pdf' % strgmapsplnk[k], mapsplnk, satu=True)

            if subspnts:
                print 'Subtracting point sources...'

                # subtract PSs from the Planck maps
                ## read PCCS
                if strgmapsplnk[k][1] == '0':
                    strg = 'R2.04'
                else:
                    strg = 'R2.01'
                dataplnk = pf.getdata(pathbase + '/COM_PCCS_%s_%s.fits' % (strgmapsplnk[k][1:], strg), 1)
                fluxpntsplnk = dataplnk['GAUFLUX'] * 1e-3 # [Jy]
                lgalpntsplnk = dataplnk['GLON'] # [deg]
                lgalpntsplnk = (lgalpntsplnk - 180.) % 360. - 180.
                bgalpntsplnk = dataplnk['GLAT'] # [deg]
                fwhmpntsplnk = dataplnk['GAU_FWHM_EFF'] / 60. # [deg]
                stdvpntsplnk = fwhmpntsplnk / 2. / sqrt(2. * log(2.)) # [deg]
                numbpntsplnk = fluxpntsplnk.size
                print 'numbpntsplnk'
                print numbpntsplnk

                ## calculate PS map using PCCS
                pathmapspntsplnk = pathdata + 'mapspntsplnk%s.fits' % strgmapsplnk[k]
                if os.path.isfile(pathmapspntsplnk):
                    print 'Reading Planck PS map for %d GHz from %s' % (float(strgmapsplnk[k]), pathmapspntsplnk)
                    mapspntsplnk = pf.getdata(pathmapspntsplnk)
                else:
                    
                    print 'numbside'
                    print numbside
                    print
                    mapspntsplnk = tdpy.util.retr_mapspnts(lgalpntsplnk, bgalpntsplnk, stdvpntsplnk, fluxpntsplnk, verbtype=1, numbside=numbside)
                    tdpy.util.plot_maps(pathimag + 'mapspntsplnk%s.pdf' % strgmapsplnk[k], mapspntsplnk, satu=True)
                    pf.writeto(pathmapspntsplnk, mapspntsplnk, clobber=True)
        
                ## plot the PCSC map
                tdpy.util.plot_maps(pathimag + 'mapspntsplnk%s.pdf' % strgmapsplnk[k], mapspntsplnk, satu=True)

                mapsorigplnk = mapsplnk - mapspntsplnk
            else:
                mapsorigplnk = mapsplnk

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
        path = pathimag + 'tran_%s.pdf' % rtag
        plt.savefig(path)
        plt.close(figr)
    
    # get Planck PS mask
    path = pathdata + 'HFI_Mask_PointSrc_2048_R2.00.fits'
    mapsmask = pf.open(path)[1].data['F353']
    mapsmask = hp.reorder(mapsmask, n2r=True)
    indxpixlmask = where(mapsmask == 1)
    tdpy.util.plot_maps(pathimag + 'mapsmask_%s.pdf' % rtag, mapsmask, satu=True)
    
    ## multipole
    maxmmpol = 3. * numbside - 1.
    numbalmc = int(maxmmpol * (maxmmpol + 1.) / 2. + maxmmpol + 1)
    numbmpol = int(maxmmpol) + 1
    mpol = arange(numbmpol)
    mpolgrid, temp = hp.Alm.getlm(lmax=maxmmpol)

    # get input maps
    ## Planck map
    tdpy.util.plot_maps(pathimag + 'mapsplnkorig_%s.pdf' % rtag, mapsplnkorig, satu=True)
    
    path = pathdata + 'mapsplnk_%s.fits' % rtag
    if os.path.isfile(path):
        mapsplnk = pf.getdata(path)
    else:
        mapsplnkorig -= mean(mapsplnkorig)
        mapsplnkorig /= std(mapsplnkorig)
        mapsplnk, mapsalmc, mpolsmthplnk, wghtsmthplnk = smth(mapsplnkorig, mpolsmth, mpol=True, retrfull=True, indxpixlmask=indxpixlmask, numbsideoutp=numbside)
        tdpy.util.plot_maps(pathimag + 'mapsplnk_%s.pdf' % rtag, mapsplnk)
        #mapsplnk = hp.ud_grade(mapsplnksmth, numbside)
        pf.writeto(path, mapsplnk, clobber=True)

    almcplnktemp = hp.map2alm(mapsplnk)

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

    # Gaussian noise map
    mapsgaus = 0.1 * randn(numbpixl)

    ## Fermi Diffuse Model
    # temp
    mapsfdfmorig = tdpy.util.retr_fdfm(binsener, numbside=numbside)
    mapsfdfmorig -= mean(mapsfdfmorig, 1)[:, None]
    mapsfdfmorig /= std(mapsfdfmorig, 1)[:, None]
    mapsfdfm = empty_like(mapsfdfmorig)
    almcfdfmorig = empty((numbener, numbalmc), dtype=complex) 
    almcfdfm = empty((numbener, numbalmc), dtype=complex) 
    for i in arange(numbener):
        almcfdfmorig[i, :] = hp.map2alm(mapsfdfmorig[i, :])
        tdpy.util.plot_maps(pathimag + 'mapsfdfmorig%04d_%s.pdf' % (i, rtag), mapsfdfmorig[i, :])
        mapsfdfm[i, :], almcfdfm[i, :], mpolsmthfdfm, wghtsmthfdfm = smth(mapsfdfmorig[i,:], mpolsmth, mpol=True, retrfull=True)

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

    axis.axvline(numbside, ls='--', color='black', alpha=alph, label='$N_{side}$')
    axis.axvline(mpolmerg, ls='-.', color='black', alpha=alph, label='$l_{merg}$')
    axis.set_ylabel('$w_l$')
    axis.set_xlabel('$l$')
    axis.set_ylim([1e-4, 1.])
    axis.set_xlim([amin(mpol), amax(mpol)])
    axis.legend(loc=2)
    plt.tight_layout()
    path = pathimag + 'wght_%s.pdf' % rtag
    plt.savefig(path)
    plt.close(figr)
   
    # compute the power spectra
    factmpol = mpol * (2. * mpol + 1.) / 4. / pi
    psecplnktemp = factmpol * hp.anafast(mapsplnk)
    psecgaus = factmpol * hp.anafast(mapsgaus)
    psecfdfm = empty((numbener, numbmpol))
    psecplnk = empty((numbener, numbmpol))
    indxmpoltemp = where(mpol < 10.)[0]
    almcplnk = empty((numbener, numbalmc), dtype=complex)
    for i in arange(numbener):
        psecfdfm[i, :] = factmpol * hp.anafast(mapsfdfm[i, :])
        ## correct the Planck variance
        # temp
        #factcorr = sum(psecfdfm[i, indxmpoltemp]) / sum(psecplnktemp[indxmpoltemp])
        factcorr = 1.
        psecplnk[i, :] = factcorr * psecplnktemp
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
        
        axis.loglog(mpol, psecfdfm[i, :], label='FDM')
        axis.loglog(mpol, psecplnk[i, :], label='Planck')
        axis.loglog(mpol, psecmerg[i, :], label='Merged')
        axis.loglog(mpol, psecgaus, label='Uncorrelated Noise', alpha=alph, ls='--')
        
        axis.axvline(numbside, ls='--', color='black', alpha=alph, label='$N_{side}$')
        axis.axvline(mpolmerg, ls='-.', color='black', alpha=alph, label='$l_{merg}$')
        axis.axvline(mpolsmth, color='black', alpha=alph, label='$l_{smth}$')
        axis.set_ylabel('$l(2l+1)C_l/4\pi$')
        axis.set_xlabel('$l$')
        axis.set_ylim([1e-3, 1.])
        axis.set_xlim([amin(mpol), amax(mpol)])
        axis.legend(loc=3, ncol=2)
        plt.tight_layout()
        path = pathimag + 'psec%04d_%s.pdf' % (i, rtag)
        plt.savefig(path)
        plt.close(figr)

        for plotigal in [False, True]:

            if plotigal:
                minmlgal = -20.
                maxmlgal = 20.
                minmbgal = -20.
                maxmbgal = 20.
                strg = 'igal'
            else:
                minmlgal = -180.
                maxmlgal = 180.
                minmbgal = -90.
                maxmbgal = 90.
                strg = 'ngal'
               
            path = pathimag + '%s/mapsfdfm%04d_%s.pdf' % (strg, i, rtag)
            tdpy.util.plot_maps(path, mapsfdfm[i, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
            
            if i == 0:
                path = pathimag + '%s/mapsplnk_%s.pdf' % (strg, rtag)
                tdpy.util.plot_maps(path, mapsplnk, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
            
            path = pathimag + '%s/mapsmerg%04d_%s.pdf' % (strg, i, rtag)
            tdpy.util.plot_maps(path, mapsmerg[i, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
            
            path = pathimag + '%s/mapsresifdfm%04d_%s.pdf' % (strg, i, rtag)
            tdpy.util.plot_maps(path, mapsmerg[i, :] - mapsfdfm[i, :], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, resi=True, satu=True)
            
            path = pathimag + '%s/mapsresiplnk%04d_%s.pdf' % (strg, i, rtag)
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


def retr_rtag(gdat):
    
    gdat.rtag = '%s' % (gdat.datatype)


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


def retr_ener(gdat):

    # temp
    gdat.binsenerfull = array([0.1, 0.3, 1., 3., 10., 100.])
    gdat.meanenerfull = sqrt(gdat.binsenerfull[1:] * gdat.binsenerfull[:-1])
    gdat.numbenerfull = gdat.binsenerfull.size - 1
    gdat.indxenerfull = arange(gdat.numbenerfull)
    gdat.binsener = gdat.binsenerfull[gdat.indxenerincl[0]:gdat.indxenerincl[-1] + 2]
    gdat.meanener = sqrt(gdat.binsener[1:] * gdat.binsener[:-1])
    gdat.diffener = gdat.binsener[1:] - gdat.binsener[:-1]
    gdat.numbener = gdat.meanener.size
    gdat.indxener = arange(gdat.numbener)
    gdat.strgbinsener = ['%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]) for i in gdat.indxener]
    

def regrback( \
              numbproc=1, \
              numbswep=None, \
              datatype='inpt', \
              verbtype=1, \
              makeplot=False, \
              strgexpr='fermflux_cmp0_igal.fits', \
              strgexpo='fermexpo_cmp0_igal.fits', \
              strgback=['isotflux', 'fdfmflux', 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits', 'wssa_sample_1024.fits', 'lambda_sfd_ebv.fits', 'darktemp'], \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(3, 4), \
              maxmgang=deg2rad(20.)
             ):

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

    gdat.numbback = len(strgback)
    gdat.indxback = arange(gdat.numbback)

    retr_ener(gdat)

    ## event type
    gdat.evttfull = array([4, 8, 16, 32])
    gdat.numbevttfull = gdat.evttfull.size
    gdat.indxevttfull = arange(gdat.numbevttfull)
    
    gdat.evtt = gdat.evttfull[gdat.indxevttincl]
    gdat.numbevtt = gdat.evtt.size
    gdat.indxevtt = arange(gdat.numbevtt)

    boolsmth = False

    ## pixelization
    gdat.numbside = 256
    gdat.numbpixlfull = gdat.numbside**2 * 12
    gdat.lgalheal, gdat.bgalheal, gdat.numbpixl, gdat.apix = tdpy.util.retr_healgrid(gdat.numbside)
    gdat.indxpixlrofi = where((abs(gdat.lgalheal) < gdat.maxmgang) & (abs(gdat.bgalheal) < gdat.maxmgang))[0]
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.lgalgrid = gdat.lgalheal[gdat.indxpixlrofi]
    gdat.bgalgrid = gdat.bgalheal[gdat.indxpixlrofi]
       
    minmlgal = -gdat.maxmgang * 180. / pi
    maxmlgal = gdat.maxmgang * 180. / pi
    minmbgal = -gdat.maxmgang * 180. / pi
    maxmbgal = gdat.maxmgang * 180. / pi

    boolplotpsec = True

    indxdatacubefilt = meshgrid(gdat.indxenerincl, gdat.indxpixlrofi, gdat.indxevttincl, indexing='ij')
    
    # get data structure
    datapara = retr_datapara(gdat)
    
    if gdat.numbswep == None:
        gdat.numbswep = 4 * gdat.numbpara * 10000

    # setup
    retr_rtag(gdat)
    
    ## data
    if gdat.datatype == 'inpt':
        
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + gdat.strgexpr
        gdat.exprflux = pf.getdata(path)
        
        ### filter
        gdat.dataflux = gdat.exprflux[indxdatacubefilt]

    ## exposure
    path = os.environ["PCAT_DATA_PATH"] + '/' + gdat.strgexpo
    gdat.expo = pf.getdata(path)
    gdat.expo = gdat.expo[indxdatacubefilt]
   
    gdat.datacnts = gdat.dataflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]

    # power spectrum calculation
    gdat.numbmapsplot = gdat.numbback + 1
    gdat.mapsplot = empty((gdat.numbmapsplot, gdat.numbpixlfull))
    
    gdat.mapsplot[0, gdat.indxpixlrofi] = sum(gdat.datacnts[1, :, :], 1)

    ## templates
    gdat.fluxback = empty((gdat.numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:

        if boolsmth:
            strg = '_smth'
        else:
            strg = ''
        if c == 0:
            strg = 'isotflux'
        if c == 1:
            strg = 'fdfmflux' + strg
        if c == 2:
            strg = 'plnkdust' + strg
        if c == 3:
            strg = 'wisestar' + strg
        if c == 4:
            strg = 'finkdust' + strg
        if c == 5:
            strg = 'darktemp' + strg

        # temp -- ROI should be fixed at 40 X 40 degree^2
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/' + strg + '.fits'
        if os.path.isfile(path) or not boolplotpsec:
            fluxback = pf.getdata(path)
        else:
            
            if c == 0:
                fluxbackorig = tdpy.util.retr_isot(gdat.binsenerfull)
            if c == 1:
                fluxbackorig = tdpy.util.retr_fdfm(gdat.binsenerfull) 
            if c == 2:
                pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + gdat.strgback[c]
                fluxbackorig = pf.getdata(pathtemp, 1)['RADIANCE']
                fluxbackorig = hp.ud_grade(fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 3:
                pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + gdat.strgback[c]
                fluxbackorig = pf.getdata(pathtemp, 0)
                fluxbackorig = hp.ud_grade(fluxbackorig, gdat.numbside, order_in='RING', order_out='RING')
            if c == 4:
                pathtemp = os.environ["FERM_IGAL_DATA_PATH"] + gdat.strgback[c]
                fluxbackorig = pf.getdata(pathtemp)['TEMPERATURE']
                fluxbackorig = hp.ud_grade(fluxbackorig, gdat.numbside, order_in='NESTED', order_out='RING')
            if c == 5:
                fluxbackorig = tdpy.util.retr_nfwp(1., gdat.numbside)

            # normalize the templates
            if c != 0 and c != 1:
                fluxbackorig /= mean(fluxbackorig[gdat.indxpixlrofi])
            
            # smooth the templates
            fluxback = empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            for m in gdat.indxevttfull:
                if c == 0 or c == 1:
                    fluxback[:, :, m] = fluxbackorig
                else:
                    for i in gdat.indxenerfull:
                        fluxback[i, :, m] = fluxbackorig
            
            if boolsmth:
                fluxback = tdpy.util.smth_ferm(fluxback, gdat.meanenerfull, gdat.indxevttfull)
            # temp
            fluxback[where(fluxback < 0.)] = 0.

            pf.writeto(path, fluxback, clobber=True)

            # load the map to the array whose power spectrum will be calculated
            if c == 0:
                gdat.mapsplot[c+1, :] = fluxbackorig[1, :]
            elif c == 1:
                gdat.mapsplot[c+1, :] = fluxbackorig[1, :]
            else:
                gdat.mapsplot[c+1, :] = fluxbackorig
    
        # take only the energy bins, spatial pixels and event types of interest
        fluxback = fluxback[indxdatacubefilt]
        indxdatacubetemp = meshgrid(array([c]), gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        gdat.fluxback[indxdatacubetemp] = fluxback

    gdat.pathbase = os.environ["FERM_IGAL_DATA_PATH"]
    gdat.pathimag = gdat.pathbase + '/imag/regrback/%s/' % gdat.rtag
    cmnd = 'mkdir -p ' + gdat.pathimag
    os.system(cmnd)

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
    chan = tdpy.mcmc.init(retr_llik, datapara, numbproc=gdat.numbproc, numbswep=gdat.numbswep, initsamp=initsamp, gdatextr=gdat, \
                verbtype=gdat.verbtype, pathbase=gdat.pathbase, rtag=gdat.rtag, numbplotside=numbplotside)
    
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


def ferm_expr_igal(strgexpr='fermflux_cmp0_igal.fits', strgexpo='fermexpo_cmp0_igal.fits'):
    
    pcat.main.init( \
              psfntype='doubking', \
              numbswep=2000000, \
              randinit=False, \
              maxmgang=deg2rad(20.), \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(2, 4), \
              minmflux=3e-11, \
              maxmflux=3e-6, \
              pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
              regitype='igal', \
              strgback=['isotflux.fits', 'fdfmflux.fits'], \
              strgexpo=strgexpo, \
              datatype='inpt', \
              strgexpr=strgexpr, \
             )
    
    
def ferm_mock_igal_brok():
     
    listmockfluxdistbrek = array([1e-10, 3e-10, 1e-9, 3e-9, 1e-8])
    listmockfluxdistsloplowr = array([1.9, 2.2, 2.8, 3.1, 3.4])

    for fluxdistbrek in listfluxdistbrek:
    
        pcat.main.init( \
                       pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                       numbswep=10000, \
                       randinit=False, \
                       indxevttincl=arange(2, 4), \
                       indxenerincl=arange(1, 4), \
                       strgexpo='fermexpo_cmp0_igal.fits', \
                       strgback=['isotflux.fits'], \
                       regitype='igal', \
                       
                       maxmgang=deg2rad(20.), \
                       fluxdisttype=['brok'], \
                       psfntype='doubking', \
                       
                       maxmnumbpnts=array([400]), \
                       minmflux=3e-11, \
                       maxmflux=3e-7, \

                       datatype='mock', \
                       mocknumbpnts=array([400]), \
                       
                       mockfluxdistbrek=array([1e-9]), \
                       mockfluxdistsloplowr=array([2.6]), \
                       mockfluxdistslopuppr=array([1.6]), \
                       
                       mocksinddiststdv=array([.5]), \
                       mocksinddistmean=array([2.]), \
                      )


def ferm_mock_igal_popl():
     
    pcat.main.init( \
                   pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                   numbswep=10000, \
                   randinit=False, \
                   indxevttincl=arange(2, 4), \
                   indxenerincl=arange(1, 4), \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits'], \
                   regitype='igal', \
                   
                   maxmgang=deg2rad(20.), \
                   fluxdisttype=['brok'], \
                   psfntype='doubking', \
                   
                   maxmnumbpnts=array([400]), \
                   minmflux=3e-11, \
                   maxmflux=3e-7, \

                   datatype='mock', \
                   mocknumbpnts=array([400]), \
                   
                   mockspatdisttype=['unif', 'disc', 'gang'], \
                   
                   mocksinddiststdv=array([.5]), \
                   mocksinddistmean=array([2.]), \
                  )


def ferm_mock_igal():
     
    pcat.main.init( \
                   psfntype='doubking', \
                   numbswep=1000000, \
                   randinit=False, \
                   indxevttincl=arange(2, 4), \
                   indxenerincl=arange(1, 4), \
                   strgexpo='fermexpo_cmp0_igal.fits', \
                   strgback=['isotflux.fits', 'fdfmflux.fits'], \
                   pathdata=os.environ["FERM_IGAL_DATA_PATH"], \
                   regitype='igal', \
                   
                   #maxmnumbpnts=array([2, 2, 4]), \
                   maxmnumbpnts=array([200, 200, 400]), \
                   maxmgang=deg2rad(20.), \
                   minmflux=3e-11, \
                   maxmflux=3e-7, \
                   
                   sinddiststdv=array([.5, .5, .5]), \
                   sinddistmean=array([2., 2., 2.]), \
                   
                   datatype='mock', \
                   #mocknumbpnts=array([1, 1, 2]), \
                   mocknumbpnts=array([100, 100, 200]), \
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

