
# coding: utf-8

# In[1]:

import json, requests
import h5py
from numpy import *
import healpy as hp
import tdpy_util
import matplotlib.pyplot as plt
import os
import pyfits as pf
# plotting
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
sns.set(context='poster', style='dark', color_codes=True)


# In[2]:

def query(lon, lat, coordsys='gal', mode='full'):
    
    url = 'http://argonaut.skymaps.info/gal-lb-query-light'
    
    payload = {'mode': mode}
    
    if coordsys.lower() in ['gal', 'g']:
        payload['l'] = lon
        payload['b'] = lat
    elif coordsys.lower() in ['equ', 'e']:
        payload['ra'] = lon
        payload['dec'] = lat
    else:
        raise ValueError("coordsys '{0}' not understood.".format(coordsys))
    
    headers = {'content-type': 'application/json'}
    
    r = requests.post(url, data=json.dumps(payload), headers=headers)
    
    try:
        r.raise_for_status()
    except requests.exceptions.HTTPError as e:
        print('Response received from Argonaut:')
        print(r.text)
        raise e
    
    return json.loads(r.text)


# In[3]:

get_ipython().magic(u'matplotlib inline')
magv = arange(10, 22)
prlxerrr = array([4., 4., 4.2, 6.0, 9.1, 14.3, 23.1, 38.8, 69.7, 138., 312., 1786.])

fig, ax = plt.subplots()
ax.plot(magv, prlxerrr)
ax.set_ylabel(r'$\sigma_\pi [\mu$as]')
ax.set_xlabel('V [mag]')
ax.set_yscale('log')
plt.show()


# In[17]:

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



# In[18]:

nside = 256
maxmgang = 10.
lghp, bghp, nside, npixl, apix = tdpy_util.get_heal(nside) 
mpixl = where((abs(lghp) < maxmgang) & (abs(90. - bghp) < maxmgang))[0]

#qresult = query(list(lghp[mpixl]), list(bghp[mpixl]))
for key in qresult.keys():
    fig, ax = plt.subplots()
    ax.hist(qresult[key], 100)
    ax.set_title(key)
    plt.show()

#qresult = query(list(lghp[mpixl]), list(bghp[mpixl]), mode='sfd')
for key in qresult.keys():
    fig, ax = plt.subplots()
    ax.hist(qresult[key])
    ax.set_title(key)
    plt.show()
    
#qresult = query(list(lghp[mpixl]), list(bghp[mpixl]), mode='lite')
for key in qresult.keys():
    fig, ax = plt.subplots()
    ax.hist(qresult[key])
    ax.set_title(key)
    plt.show()


# In[ ]:



