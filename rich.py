import numpy as np 
import matplotlib.pyplot as plt 
import astropy.stats
from astropy.io import fits 
from astropy.table import hstack
from __init__ import *

pathdata = os.environ['TDGU_DATA_PATH'] + '/xray_back/data'

hstow = fits.open(pathdata + '/acis_DE_01236_stowed_evt_110210.fits')
hstow.info()
datastow = hstow[1].data
hstow.close()

x_coords = datastow['x']
y_coords = datastow['y']
e_values = datastow['energy']
cropx_x = []
cropx_y = []
cropx_e = []
cropy_x = [] #these three are final cropped image for DE
cropy_y = []
cropy_e = []

NBINS = 200 #arbitrary binning

for i in range(0,len(x_coords)):
    a = x_coords[i]
    if a > 3009:
        if a < 5099:
            cropx_x.append(a)
            cropx_y.append(y_coords[i])
            cropx_e.append(e_values[i])
for j in range(0, len(cropx_y)):
    b = cropx_y[j]
    if b > 2109:
        if b < 4200:
            cropy_x.append(cropx_x[j])
            cropy_y.append(b)
            cropy_e.append(cropx_e[j])

img_zero, yedges, xedges = np.histogram2d(cropy_x, cropy_y, NBINS)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.figure(figsize=(10,10))
plt.imshow(img_zero, extent=extent, interpolation='nearest', cmap='gist_yarg', origin='lower')
plt.xlabel('Chip Y Coordinates (pixel)')
plt.ylabel('Chip X Coordinates (pixel)')
plt.savefig('stowimag.pdf')
plt.close()

# In[15]:

eV_to_keV = 0.001
e_low = 0.5 #keV
e_high = 7 #keV
exp_time = 6.03e5
div = (np.min(cropy_e)+np.max(cropy_e))/NBINS*eV_to_keV

traditional_bins = np.array([500, 2000, 7000]) #for 0.5-2 and 2-7 keV
traditional_norm = np.diff(traditional_bins)*eV_to_keV #will divide to get per keV
bin_com = np.sqrt(traditional_bins[1:]*traditional_bins[:-1]) #geometric means
bin_av = np.array([1250,4500]) #arithmetic means

print(traditional_norm, bin_com)


energyhist = np.histogram(cropy_e, traditional_bins)
# stowed_energy = np.asarray(energyhist[1], dtype='float')*eV_to_keV
counts = energyhist[0]
print(bin_com, counts)

plt.figure(figsize=(10,8))
plt.plot(bin_com*eV_to_keV, counts/exp_time/traditional_norm, color='black', linewidth=1)
plt.plot(bin_com*eV_to_keV, counts/exp_time/traditional_norm, 'ro')
plt.xscale('log')
plt.yscale('log')
plt.xlim(e_low,e_high)
plt.ylim(5e-2,5e-1)
plt.xlabel('Energy (keV)', fontsize=14)
plt.ylabel('S (Counts/s/keV)', fontsize=14)
plt.title('X-ray Spectrum of Chandra Stowed Background (DE): ' + str(e_low) + ' - ' + str(e_high) + 'keV', fontsize=15)
plt.savefig('stowed_spectra_300_10000_eV.pdf')


# In[34]:

with open(pathdata + '/effective_area_cycle3.txt', 'r') as p:
    energy = []
    ea = []
    for line in p:
        line = line.split()
        energy.append(line[0])
        ea.append(line[1])

energy_array = np.asarray(energy, dtype='float') #in keV
ea = np.asarray(ea, dtype='float') #cm^2

plt.plot(energy_array, ea)
plt.xscale('log')
plt.xlim(0.3,11)
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm^2)')
plt.title('Effective Area of ACIS-I Detector')

bin1 = []
bin2 = []
for a in range(0, len(energy_array)):
    p = energy_array[a]
#     print(p)
    if p > 0.5:
        if p <= 2:
            bin1.append(ea[a])
        elif p <= 7:
            bin2.append(ea[a])

ea_bin1_weighted = np.sum(bin1)/len(bin1) # average value over range in cm^2
ea_bin2_weighted = np.sum(bin2)/len(bin2)

print(ea_bin1_weighted, ea_bin2_weighted)

ea_weighted = [ea_bin1_weighted, ea_bin2_weighted]

total_pixel_area = (2048*0.492)**2 #in arcsecond^2
arcsec_to_sterad = 4.25e10 
total_steradians = total_pixel_area/arcsec_to_sterad
print(total_steradians)


# In[35]:

plt.figure(figsize=(10,8))
plt.plot(bin_com*eV_to_keV, counts/exp_time/traditional_norm/total_steradians/ea_weighted, color='black', linewidth=1)
plt.plot(bin_com*eV_to_keV, counts/exp_time/traditional_norm/total_steradians/ea_weighted, 'ro')
plt.xscale('log')
plt.yscale('log')
plt.xlim(e_low,e_high)
# plt.ylim(5e-2,5e-1)
plt.xlabel('Energy (keV)', fontsize=14)
plt.ylabel('S (Counts/s/keV/sr/cm^2)', fontsize=14)
plt.title('X-ray Spectrum of Chandra Stowed Background (DE): ' + str(e_low) + ' - ' + str(e_high) + 'keV', fontsize=15)
plt.savefig('stowed_spectra_300_10000_eV.pdf')



