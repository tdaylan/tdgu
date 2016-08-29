from __init__ import *

numbside = 256    
lgalheal, bgalheal, numbpixl, apix = tdpy.util.retr_healgrid(numbside)

maps = zeros(numbpixl)
maps = exp(-0.5 * (bgalheal / 5.)**2) + 3. * exp(-0.5 * (sqrt(lgalheal**2 + bgalheal**2) / 10.)**2)

numbpnts = 100
indxpixl = choice(arange(numbpixl), size=numbpnts)
#mapsmask = zeros(numbpixl, dtype=bool)
#mapsmask[indxpixl] = False

mapspnts = copy(maps)
mapspnts[indxpixl] = rand(numbpnts)

mapspntsmask = hp.ma(mapspnts)

tdpy.util.plot_maps('./maps.pdf', maps, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps('./mapspnts.pdf', mapspnts, numbsidelgal=1000, numbsidebgal=1000, satu=True)
tdpy.util.plot_maps('./mapspntsmask.pdf', mapspntsmask, numbsidelgal=1000, numbsidebgal=1000, satu=True)
        
print 'mapspntsmask.mask'
print amin(mapspntsmask.mask)
print amax(mapspntsmask.mask)
print 'mapspntsmask.filled'
print amin(mapspntsmask.filled)
print amax(mapspntsmask.filled)
print 'mapspntsmask'
print amin(mapspntsmask)
print amax(mapspntsmask)

figr, axis = plt.subplots()

mpol = arange(3 * numbside)
fact = mpol * (2 * mpol + 1.) / 4. / pi
axis.loglog(mpol, fact * hp.anafast(maps), label='Diffuse')
axis.loglog(mpol, fact * 2. * hp.anafast(mapspnts), label='Diffuse + PS')

#mapsmask = zeros(numbpixl)
#mapsmask = zeros(numbpixl, dtype=bool)
#mapspntsmask.mask = mapsmask
mapspntsmask = copy(mapspnts)
mapspntsmask = hp.ma(mapspntsmask)
axis.loglog(mpol, fact * 3. * hp.anafast(mapspntsmask), label='Diffuse + PS, Masked (None)')

#mapsmask = zeros(numbpixl)
#mapsmask = zeros(numbpixl, dtype=bool)
#mapsmask[indxpixl[:numbpnts / 2]] = True
#mapspntsmask.mask = mapsmask
mapspntsmask = copy(mapspnts)
mapspntsmask[indxpixl[:numbpnts / 2]] = hp.UNSEEN
mapspntsmask = hp.ma(mapspntsmask)
axis.loglog(mpol, fact * 4. * hp.anafast(mapspntsmask), label='Diffuse + PS, Masked (Half)')

#mapsmask = zeros(numbpixl)
#mapsmask = zeros(numbpixl, dtype=bool)
#mapsmask[indxpixl] = True
#mapspntsmask.mask = mapsmask
mapspntsmask = copy(mapspnts)
masppntsmask[indxpixl] = hp.UNSEEN
mapspntsmask = hp.ma(mapspntsmask)
axis.loglog(mpol, fact * 5. * hp.anafast(mapspntsmask), label='Diffuse + PS, Masked (Full)')

axis.set_ylabel('$l(2l+1)C_l/4\pi$')
axis.set_xlabel('$l$')
axis.set_xlim([amin(mpol), amax(mpol)])
axis.legend(loc=3, ncol=2)
plt.tight_layout()
path = './psec.pdf'
plt.savefig(path)
plt.close(figr)

    
