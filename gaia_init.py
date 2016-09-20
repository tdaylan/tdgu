from tdgu.__init__ import *
import tdpy.util

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


#pathdata = '/n/fink2/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/fits/'

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
    
    
