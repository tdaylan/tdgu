from __init__ import *

def read_pstr(lgal, bgal):

    numbpnts = lgal.size
    indxpixl = empty(numbpnts, dtype=int)
    for k in range(numbpnts):
        indxpixl[k] = ang2pix(32., pi / 2. - deg2rad(bgal), deg2rad(lgal)) 
    indxpixl = unique(indxpixl)

    for k in range(indxpixl):
        path = '/n/fink2/dfink/decam-ucal-qz-time/chunks-qz-star-v3/ps1-%05d.fits' % indxpixl[k]
        print 'Reading %s' % path
        catl = pf.getdata(path)

    return pstrcatl


def comp():

    strgproc = os.uname()[1]
    if strgproc == 'fink1':
        pathimag = '/n/pan/www/tansu/imag/gaia_init/'
        pathdata = '/n/fink2/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/fits/'
        strg = 'GaiaSource_000-000-000.fits'
        #pathdata = '/n/fink2/gaia/gaia-sweep.fits'
    else:
        pathbase = os.environ["TDGU_DATA_PATH"]
        pathdata = pathbase + '/data/gaia_init/'
        pathimag = pathbase + '/imag/gaia_init/'
        os.system('mkdir -p %s' % pathimag)
        os.system('mkdir -p %s' % pathdata)
    
    rasc = rand(100) * 360.
    decl = (rand(100) - 0.5) * 180.
    
    #tdpy.util.read_fits(pathdata + strg, pathimag=pathimag)
    
    
    strg = 'chunk-00000.fits'
    rascgaia = pf.getdata(pathdata + strg, 1)['ra']
    declgaia = pf.getdata(pathdata + strg, 1)['dec']
    lgalgaia, bgalgaia = tdpy.util.retr_galcfromequc(rascgaia, declgaia)
    fluxgaia = pf.getdata(pathdata + strg, 1)['phot_g_mean_mag']
    
    #numbtest = 100
    #lgalgaia = lgalgaia[:numbtest]
    #bgalgaia = bgalgaia[:numbtest]
    #fluxgaia = fluxgaia[:numbtest]
    
    strg = 'ps1-00000.fits'
    rascpstr = pf.getdata(pathdata + strg, 1)['ra']
    declpstr = pf.getdata(pathdata + strg, 1)['dec']
    lgalpstrorig, bgalpstrorig = tdpy.util.retr_galcfromequc(rascpstr, declpstr)
    specpstr = pf.getdata(pathdata + strg, 1)['MEAN'].T
    numbband = specpstr.shape[0]
    
    strgband = ['g', 'r', 'i', 'z', 'y']
    
    # plot histograms of fluxes and positions
    figr, axis = plt.subplots()
    axis.set_ylabel('$N$')
    axis.set_xlabel('$l$ [deg]')
    axis.hist(lgalpstrorig, label='PS', alpha=0.5)
    axis.hist(lgalgaia, label='Gaia', alpha=0.5)
    plt.tight_layout()
    path = pathimag + 'histbgal.pdf'
    plt.savefig(path)
    plt.close(figr)
    
    figr, axis = plt.subplots()
    axis.set_ylabel('$N$')
    axis.set_xlabel('$b$ [deg]')
    axis.hist(bgalpstrorig, label='PS', alpha=0.5)
    axis.hist(bgalgaia, label='Gaia', alpha=0.5)
    plt.tight_layout()
    path = pathimag + 'histbgal.pdf'
    plt.savefig(path)
    plt.close(figr)
    
    figr, axis = plt.subplots()
    axis.set_ylabel('$N$')
    axis.set_xlabel('$f$ [mag]')
    minm = min(amin(fluxgaia), amin(specpstr))
    maxm = max(amax(fluxgaia), amax(specpstr))
    bins = linspace(minm, maxm, 50)
    for i in range(numbband):
        axis.hist(specpstr[i, :], label='PS, %s band' % strgband[i], lw=3, bins=bins, histtype='step')
    axis.hist(fluxgaia, label='Gaia', bins=bins, histtype='step', lw=3)
    axis.legend(loc=2)
    plt.tight_layout()
    path = pathimag + 'histflux.pdf'
    plt.savefig(path)
    plt.close(figr)
        
    for i in range(numbband):
    
        indxpstrgood = where(specpstr[i, :] < 29.)[0]
        lgalpstr = lgalpstrorig[indxpstrgood]
        bgalpstr = bgalpstrorig[indxpstrgood]
        fluxpstr = specpstr[i, indxpstrgood]
    
        indxpstr = tdpy.util.corr_catl(lgalpstr, bgalpstr, fluxpstr, lgalgaia, bgalgaia, fluxgaia, verbtype=2)
        indxindxpstr = where(indxpstr > 0)[0]
        lgalpstr = lgalpstr[indxpstr[indxindxpstr]]
        bgalpstr = bgalpstr[indxpstr[indxindxpstr]]
        fluxpstr = fluxpstr[indxpstr[indxindxpstr]]
        lgalgaia = lgalgaia[indxindxpstr]
        bgalgaia = bgalgaia[indxindxpstr]
        fluxgaia = fluxgaia[indxindxpstr]
        
        figr, axis = plt.subplots()
        axis.set_ylabel('$l_{P1}$ [deg]')
        axis.set_xlabel('$l_{G}$ [deg]')
        axis.scatter(lgalgaia, lgalpstr, s=0.4)
        coef, func, strg = tdpy.util.regr(lgalgaia, lgalpstr, 1)
        plt.plot(lgalgaia, func(lgalgaia), lw=1)
        axis.text(0.05, 0.9, strg, transform=axis.transAxes)
        axis.ticklabel_format(useOffset=False)
        plt.tight_layout()
        path = pathimag + 'corrlgal_%s.pdf' % strgband[i]
        plt.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots()
        axis.set_ylabel('$l_{P1}/l_{G} - 1$ [%]')
        axis.set_xlabel('$l_{G}$ [deg]')
        axis.scatter(lgalgaia, 100. * lgalpstr / lgalgaia - 100., s=0.4)
        plt.tight_layout()
        path = pathimag + 'difflgal_%s.pdf' % strgband[i]
        plt.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots()
        axis.set_ylabel('$b_{P1}/b_{G} - 1$ [%]')
        axis.set_xlabel('$b_{G}$ [deg]')
        axis.scatter(bgalgaia, 100. * bgalpstr / lgalgaia - 100., s=0.4)
        plt.tight_layout()
        path = pathimag + 'diffbgal_%s.pdf' % strgband[i]
        plt.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots()
        axis.set_ylabel('$f_{P1}/f_{G} - 1$ [%]')
        axis.set_xlabel('$f_{G}$ [mag]')
        axis.scatter(fluxgaia, 100. * fluxpstr / fluxgaia - 100., s=0.4)
        plt.tight_layout()
        path = pathimag + 'diffflux_%s.pdf' % strgband[i]
        plt.savefig(path)
        plt.close(figr)
    

def plot_gums():
    magv = arange(10, 22)
    prlxerrr = array([4., 4., 4.2, 6.0, 9.1, 14.3, 23.1, 38.8, 69.7, 138., 312., 1786.])
    
    fig, ax = plt.subplots()
    ax.plot(magv, prlxerrr)
    ax.set_ylabel(r'$\sigma_\pi [\mu$as]')
    ax.set_xlabel('V [mag]')
    ax.set_yscale('log')
    plt.show()
    
    #RESOURCE=yCat_6137
    #Name: VI/137
    #Title: GaiaSimu Universe Model Snapshot (A.C.Robin + 2012)
    #Coosys	J2000_2010.000:	eq_FK5 J2000
    #Coosys	G:	galactic
    #Table	VI_137_gum_mw:
    #Name: VI/137/gum_mw
    #Title: Gaia Universe Model Snapshot (GUMS): Milky Way stars (among 2,143,475,885 stars)
    #Column	_Glon	(F8.4)	Galactic longitude at Epoch=J2010, proper motions taken into account  (computed by VizieR, not part of the original data)	[ucd=pos.galactic.lon]
    #Column	_Glat	(F8.4)	Galactic latitude at Epoch=J2010, proper motions taken into account  (computed by VizieR, not part of the original data)	[ucd=pos.galactic.lat]
    #Column	VMAG	(F6.3)	V-band absolute magnitude (meanAbsoluteV)	[ucd=phot.mag]
    #Column	Mbol	(F6.3)	Bolometric absolute magnitude (mbol)	[ucd=phys.magAbs.bol]
    #Column	Gmag	(F6.3)	Gaia G-band magnitude (350-1050nm) (magG)	[ucd=phot.mag;em.opt]
    #Column	GBmag	(F6.3)	Gaia Gbp band magnitude (350-770nm) (magGBp)	[ucd=phot.mag;em.opt.B]
    #Column	GRmag	(F6.3)	Gaia Grp band magnitude (650-1050nm) (magGRp)	[ucd=phot.mag;em.opt.R]
    #Column	RAJ2000	(F14.10)	Right ascension (ICRS, Epoch=J2010) (alpha)	[ucd=pos.eq.ra;meta.main]
    #Column	DEJ2000	(F14.10)	Declination (ICRS, Epoch=J2010) (delta)	[ucd=pos.eq.dec;meta.main]
    #Column	r	(F7.1)	Barycentric distance (distance)	[ucd=pos.distance;pos.heliocentric]
    #Column	V-I	(F6.3)	Intrinsic color V-I (colorVminusI)	[ucd=phot.color;em.opt.V;em.opt.I]
    #Column	Av	(F6.3)	Absorption (Av)	[ucd=phot.mag]
    #Column	Mass	(F9.4)	Stellar mass (mass)	[ucd=phys.mass]
    #Column	[Fe/H]	(F5.2)	Metallicity [Fe/H] (feH)	[ucd=phys.abund.Z]
    #Column	Teff	(I6)	Effective temperature (teff)	[ucd=phys.temperature.effective]
    #Column	logg	(F6.3)	Gravity (log) (logg)	[ucd=phys.gravity]
    #Column	fI	(I1)	[0,1] 1 (true) if interacting with a companion (flagInteracting)	[ucd=meta.code.multip]
    #_Glon|_Glat|VMAG|Mbol|Gmag|GBmag|GRmag|RAJ2000|DEJ2000|r|V-I|Av|Mass|[Fe/H]|Teff|logg|fI
    #deg|deg|mag|mag|mag|mag|mag|deg|deg|pc|mag|mag|Sun|[Sun]|K|[cm/s2]|
    
    
    path = os.environ["DUST_PRLX_PATH"] + '/dat/gums.dat'
    gums = loadtxt(path, skiprows=72)
    
    nstar = gums.shape[0]
    nattr = gums.shape[1]
    
    attrstrg = ['$l$', '$b$', '$V$', '$M_{bol}$', '$G$', '$G_b$', '$G_r$', 'RA',             'DEC', '$r$', '$V-I$', '$A_v$', 'M', 'Fe/H', '$T_{eff}$', '$\log g$', 'f_I']
    
    jattr = [0, 1, 2, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15]
    fig, axgrd = plt.subplots(6, 2, figsize=(12, 40))
    for a, axrow in enumerate(axgrd):
        for b, ax in enumerate(axrow):
            k = 2 * a + b
            ax.hist(gums[:, jattr[k]], 20)
            ax.set_xlabel(attrstrg[jattr[k]])
            plt.show
    
    fig, ax = plt.subplots()
    ax.scatter(gums[:, 0], gums[:, 1])
    ax.set_xlim([amin(gums[:, 0]), amax(gums[:, 0])])
    ax.set_ylim([amin(gums[:, 1]), amax(gums[:, 1])])
    ax.set_aspect('equal')
    plt.show()
    
    fig, ax = plt.subplots()
    ax.scatter(gums[:, 10], gums[:, 2], edgecolor='none', s=5)
    plt.show()
    

def writ_tgasdata():

    #path = os.environ["TDGU_DATA_PATH"] + '/gaia_init/data/tgashalometa.fits'
    path = os.environ["TDGU_DATA_PATH"] + '/gaia_init/data/tgashalokine.fits'
    
    ekin = pf.getdata(path)['Ek']
    amom = pf.getdata(path)['Lz']
    #amom = pf.getdata(path)['L'][:, 2]
    numbbins = 201
    minmamom = percentile(amom, 10.)
    maxmamom = percentile(amom, 90.)
    maxmamom = max(abs(minmamom), maxmamom)
    minmamom = -maxmamom
    
    minmekin = percentile(ekin, 0.)
    maxmekin = percentile(ekin, 90.)
    
    binsamom = linspace(minmamom, maxmamom, numbbins)
    binsekin = linspace(minmekin, maxmekin, numbbins)
    
    datacntstemp = histogram2d(ekin, amom, bins=(binsekin, binsamom))[0]
    
    numbside = int(sqrt(datacntstemp.size))
    datacnts = zeros((1, numbside, numbside, 1))
    datacnts[0, :, :, 0] = datacntstemp.T
    datacnts *= numbbins**2 / 4.
    
    path = os.environ["PCAT_DATA_PATH"] + '/data/inpt/tgas.fits'
    pf.writeto(path, datacnts, clobber=True)
    
    backcnts = copy(datacnts)
    backcnts[0, :, :, 0] = sp.ndimage.filters.gaussian_filter(backcnts[0, :, :, 0], sigma=15)
    backcnts[where(backcnts <= 0.)] = 1e-20
   
    set_printoptions(precision=1)

    path = os.environ["PCAT_DATA_PATH"] + '/data/inpt/tgasback.fits'
    pf.writeto(path, backcnts, clobber=True)


def pcat_tgas_mock():
    
    pcat.main.init( \
         numbswep=20000, \
         exprtype='sdyn', \
         psfninfoprio=False, \
         checprio=False, \
         strgexpo=1., \
         elemtype='clus', \
         fittback=['tgasback.fits'], \
         truenumbpnts=array([20]), \
         fittmaxmnumbpnts=array([40]), \
        )


def pcat_tgas_inpt():
    
    pcat.main.init( \
         numbswep=20000, \
         exprtype='sdyn', \
         psfninfoprio=False, \
         checprio=False, \
         optihess=False, \
         strgexpo=1., \
         elemtype='clus', \
         fittback=['tgasback.fits'], \
         strgexprflux='tgas.fits', \
         inittype='rand', \
         fittmaxmnumbpnts=array([40]), \
        )

globals().get(sys.argv[1])()


