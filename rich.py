from astropy.io import fits 
from astropy.table import hstack
from __init__ import *

pathdata = os.environ['TDGU_DATA_PATH'] + '/xray_back/data'

hstow = fits.open(pathdata + '/acis_DE_01236_stowed_evt_110210.fits')
datastow = hstow[1].data
hstow.close()

lgal = datastow['x']
bgal = datastow['y']
e_values = datastow['energy']
cropx_x = []
cropx_y = []
cropx_e = []
cropy_x = []
cropy_y = []
listenerstow = []

numbbins = 200

# counts
for i in range(len(lgal)):
    a = lgal[i]
    if a > 3009:
        if a < 5099:
            cropx_x.append(a)
            cropx_y.append(bgal[i])
            cropx_e.append(e_values[i])
for j in range(0, len(cropx_y)):
    b = cropx_y[j]
    if b > 2109:
        if b < 4200:
            cropy_x.append(cropx_x[j])
            cropy_y.append(b)
            listenerstow.append(cropx_e[j])
cnts = histogram(array(listenerstow) * 1e-3, binsener)[0]

# effective area
with open(pathdata + '/effective_area_cycle3.txt', 'r') as p:
    enereffa = []
    effatemp = []
    for line in p:
        line = line.split()
        enereffa.append(line[0])
        effatemp.append(line[1])
enereffa = asarray(enereffa, dtype='float') #in keV
effatemp = asarray(effatemp, dtype='float') #cm^2
effa = empty(numbener)
for i in indxener:
    indx = where((enereaff > binsener[i]) & (enereaff < binsener[i]))[0]
    effa[i] = trapz(eafftemp[indx], enereaff[indx])

# exposure time
expotime = 6.03e5

# energy axis
minmener = 0.5 # [keV]
maxmener = 7 # [keV]
binsener = array([0.5, 2., 7.])
deltener = diff(binsener)
meanener = sqrt(binsener[1:] * binsener[:-1])

# pixel area
apix = (0.492 / 3600. * pi / 180.)**2

# flux
flux = cnts / expotime / deltener / apix / effa

# plots
## effective area
plt.plot(enereffa, effatemp)
plt.xscale('log')
plt.xlabel('Energy [keV]')
plt.ylabel('Effective Area [cm$^2$])')

## stowed spectra
plt.figure(figsize=(6, 6))
plt.loglog(meanener, flux)
plt.xlabel('Energy [keV]')
plt.ylabel('$E^2dI/dE$ [keV/s/cm$^2$/sr]')
plt.savefig('specstow.pdf')

## stowed image
img_zero, yedges, xedges = histogram2d(cropy_x, cropy_y, numbbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.figure(figsize=(10,10))
plt.imshow(img_zero, extent=extent, interpolation='nearest', cmap='gist_yarg', origin='lower')
plt.xlabel('Chip Y Coordinates (pixel)')
plt.ylabel('Chip X Coordinates (pixel)')
plt.savefig('stowimag.pdf')
plt.close()


