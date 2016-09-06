from __init__ import *
from astroML.datasets import fetch_wmap_temperatures

pathfold = os.environ["TDGU_DATA_PATH"] + '/imag/psec/'
os.system('mkdir -p ' + pathfold)

numbside = 512
numbpixl = numbside**2 * 12
numbmpol = 3 * numbside
mpol = arange(numbmpol)
factmpol = mpol * (2. * mpol - 1.) / 4. / pi

mapswmap = fetch_wmap_temperatures(masked=False)
mapswmapmask = fetch_wmap_temperatures(masked=True)

mapsgaus = np.ma.asarray(np.random.normal(0, 0.062, mapswmapmask.shape))
mapsgaus = 0.062 * randn(numbpixl)

tdpy.util.plot_maps(pathfold + 'mapswmap.pdf', mapswmap, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps(pathfold + 'mapswmapmask.pdf', mapswmapmask, numbsidelgal=1000, numbsidebgal=1000, satu=True)

psecwmapmask = hp.anafast(mapswmapmask)
psecwmap = hp.anafast(mapswmap)
psecgaus = hp.anafast(mapsgaus)

figr, axis = plt.subplots()
axis.loglog(mpol, factmpol * psecwmapmask, label='WMAP Masked')
axis.loglog(mpol, factmpol * psecwmap, label='WMAP')
axis.loglog(mpol, factmpol * psecgaus, label='Uncorrelated Noise')
axis.set_xlabel('$l$')
axis.set_ylabel(r'$l(2l+1)/4\pi C_l$')
axis.legend()
axis.set_xlim([amin(mpol), amax(mpol)])
plt.tight_layout()
figr.savefig(pathfold + 'psecwmap.pdf')
plt.close(figr)


numbside = 512
lgalheal, bgalheal, numbpixl, apix = tdpy.util.retr_healgrid(numbside)

maps = zeros(numbpixl)
maps = exp(-0.5 * (bgalheal / 5.)**2) + 3. * exp(-0.5 * (sqrt(lgalheal**2 + bgalheal**2) / 10.)**2)

numbpnts = 1000
indxpixl = choice(arange(numbpixl), size=numbpnts)

mapspnts = copy(maps)
mapspnts[indxpixl] = rand(numbpnts)

mapspntsmask = hp.ma(mapspnts)

tdpy.util.plot_maps(pathfold + 'maps.pdf', maps, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps(pathfold + 'mapspnts.pdf', mapspnts, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps(pathfold + 'mapspntsmask.pdf', mapspntsmask, numbsidelgal=1000, numbsidebgal=1000, satu=True)
        
mpol = arange(3 * numbside)
fact = mpol * (2 * mpol + 1.) / 4. / pi

mapspntsmask = copy(mapspnts)
mask = ones(numbpixl)
mask[indxpixl] = 0.
mapspntsmask *= mask
tdpy.util.plot_maps(pathfold + 'mask.pdf', mask, numbsidelgal=1000, numbsidebgal=1000, satu=True)

mapspntsmasksmth = copy(mapspnts)
mask = ones(numbpixl)
mask[indxpixl] = 0.
mask = hp.smoothing(mask, sigma=deg2rad(0.5))
mask -= amin(mask)
mask /= amax(mask)
mapspntsmasksmth *= mask
tdpy.util.plot_maps(pathfold + 'masksmth.pdf', mask, numbsidelgal=1000, numbsidebgal=1000, satu=True)

mapspntsmaskapod = copy(mapspnts)
mask = ones(numbpixl)
maxmdist = deg2rad(5.)
for k in range(numbpnts):
    dir1 = array([lgalheal, bgalheal])
    dir2 = array([lgalheal[indxpixl[k]], bgalheal[indxpixl[k]]])
    dist = angdist(dir1, dir2, lonlat=True)
    indxpixltemp = where(dist < maxmdist)[0]
    mask[indxpixltemp] *= dist[indxpixltemp] / maxmdist
mapspntsmaskapod *= mask
tdpy.util.plot_maps(pathfold + 'maskapod.pdf', mask, numbsidelgal=1000, numbsidebgal=1000, satu=True)

figr, axis = plt.subplots()
axis.loglog(mpol, fact * hp.anafast(maps), label='Diffuse')
axis.loglog(mpol, fact * hp.anafast(mapspnts), label='Diffuse + PS', lw=1)
axis.loglog(mpol, fact * hp.anafast(mapspntsmask), label='Diffuse + PS, (Original Mask)', lw=1)
axis.loglog(mpol, fact * hp.anafast(mapspntsmasksmth), label='Diffuse + PS, (Smoothed Mask)', lw=1)
axis.loglog(mpol, fact * hp.anafast(mapspntsmaskapod), label='Diffuse + PS, (Apodized Mask)', lw=1)
axis.set_ylabel('$l(2l+1)C_l/4\pi$')
axis.set_xlabel('$l$')
axis.set_xlim([amin(mpol), amax(mpol)])
axis.legend(loc=3)
axis.set_ylim([1e-5, 1.])
plt.tight_layout()
plt.savefig(pathfold + 'psectest.pdf')
plt.close(figr)

fig = plt.figure()
hp.mollview(maps)
plt.savefig(pathfold + 'maps.pdf')
plt.show()    

fig = plt.figure()
hp.mollview(mapspnts)
plt.savefig(pathfold + 'mapspnts.pdf')
plt.show()    

fig = plt.figure()
hp.mollview(mapspntsmask)
plt.savefig(pathfold + 'mapspntsmask.pdf')
plt.show()    

fig = plt.figure()
hp.mollview(mapspntsmasksmth)
plt.savefig(pathfold + 'mapspntsmasksmth.pdf')
plt.show()    

fig = plt.figure()
hp.mollview(mapspntsmaskapod)
plt.savefig(pathfold + 'mapspntsmaskapod.pdf')
plt.show()    

mapsarry = ones((3, numbpixl))
mapsarry[1, choice(arange(numbpixl), size=100)] = rand(100)
mapsarry[2, choice(arange(numbpixl), size=10000)] = rand(10000)

figr, axis = plt.subplots()
axis.loglog(mpol, fact * hp.anafast(mapsarry[0, :]), label='0')
axis.loglog(mpol, fact * hp.anafast(mapsarry[1, :]), label='100')
axis.loglog(mpol, fact * hp.anafast(mapsarry[2, :]), label='10000')
axis.set_ylabel('$l(2l+1)C_l/4\pi$')
axis.set_xlabel('$l$')
axis.set_xlim([amin(mpol), amax(mpol)])
axis.legend(loc=3)
axis.set_ylim([1e-5, 1.])
plt.tight_layout()
plt.savefig(pathfold + 'psectest.pdf')
plt.close(figr)

