from astroML.datasets import fetch_wmap_temperatures
from __init__ import *

numbside = 512
numbpixl = numbside**2 * 12
numbmpol = 3 * numbside
mpol = arange(numbmpol)
factmpol = mpol * (2. * mpol - 1.) / 4. / pi

mapswmap = fetch_wmap_temperatures(masked=False)
mapswmapmask = fetch_wmap_temperatures(masked=True)
mapsgaus = np.ma.asarray(np.random.normal(0, 0.062, mapswmapmask.shape))
mapsgaus = 0.062 * randn(numbpixl)

print 'mapswmap'
print mapswmap
print 'mapswmapmask'
print mapswmapmask

tdpy.util.plot_maps('./wmap.pdf', mapswmap, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps('./wmapmask.pdf', mapswmapmask, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps('./wmapmaskflld.pdf', mapswmapmask.filled(), numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps('./wmapmaskmask.pdf', mapswmapmask.mask, numbsidelgal=1000, numbsidebgal=1000, satu=True)

psecwmapflld = hp.anafast(mapswmapmask.filled())
psecwmap = hp.anafast(mapswmapmask)
psecgaus = hp.anafast(mapsgaus)

figr, axis = plt.subplots()
axis.loglog(mpol, factmpol * psecwmap, label='WMAP')
axis.loglog(mpol, factmpol * 2. * psecwmapflld, label='WMAP Filled')
axis.loglog(mpol, factmpol * psecgaus, label='Uncorrelated Noise')
axis.set_xlabel('$l$')
axis.set_ylabel(r'$l(2l+1)/4\pi C_l$')
axis.legend()
axis.set_xlim([amin(mpol), amax(mpol)])
plt.tight_layout()
path = './psec.pdf'
plt.savefig(path)
plt.close(figr)

