import os, sys

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from tqdm import tqdm

import numpy as np

import time
import pandas as pd

import sim_pipeline
import sim_pipeline.image_simulation
from sim_pipeline.gg_lens import GGLens
from sim_pipeline.Pipelines.skypy_pipeline import SkyPyPipeline
import sim_pipeline.Observations.roman_speclite
import speclite
from astropy.cosmology import FlatLambdaCDM
from astropy.units import Quantity

import tdpy
import chalcedon
import pergamon
import nicomedia
from tdpy import summgene


def retr_areatest(velodisp):
    '''
    Calculate the test area around the lensing galaxy
    '''
    radieinstinft = chalcedon.retr_radieinstinft(velodisp)
    
    areatest = np.pi * (1. * radieinstinft)**2
    
    return areatest


def retr_s2nr(timeexpo):
    '''
    Roman F087 SNR at 23 AB magnitude
    https://roman.gsfc.nasa.gov/science/apttables2021/table-signaltonoise.html
    '''
    s2nr = 21. * np.sqrt(timeexpo / 55.)
    
    return s2nr


timeinitpipe = time.time()

pathbase = os.environ['TDGU_DATA_PATH'] + '/Roman_Lensing_DM_Substructure/'
pathdata = pathbase + 'data/'
pathvisu = pathbase + 'visuals/'
os.system('mkdir -p %s' % pathdata)
os.system('mkdir -p %s' % pathvisu)

area = np.linspace(1., 3000., 100)
timeexpo = np.linspace(100., 110000., 50)

areamesh, timeexpomesh = np.meshgrid(area, timeexpo)

s2nr = retr_s2nr(timeexpo)
s2nrmesh = retr_s2nr(timeexpomesh)

areadhls = 100.
areahlws = 1700.

path = pathvisu + 'AreaExposure.png'

figr, axis = plt.subplots(figsize=(8, 3.5))
listlsty = ['--', '-', '-.']
listmonttotl = np.array([12.]) # [month]
listnumbexpo = np.array([3., 13.]) # [month]
numbmonttotl = len(listmonttotl)
indxmonttotl = np.arange(numbmonttotl)
for k in indxmonttotl:
    labl = '%d-month' % listmonttotl[k]
    timetotl = listmonttotl[k] * 30. * 24. * 3600. # [sec]
    areacove = 0.28 * timetotl / timeexpo# / 2.# / listnumbexpo[k]
    axis.plot(timeexpo, areacove, color='green', ls=listlsty[k], label=labl)
axis.set_xlim([min(timeexpo), max(timeexpo)])
axis.set_ylim([min(area), max(area)])
axis.set_yscale('log')
axis.plot(2250, areahlws, ls='', color='k', marker='o', ms=10, label='HLWAS')
axis.plot(21600, 100, ls='', color='k', marker='*', ms=15, label='HLDS')
axis.plot(1e5, 5., ls='', color='k', marker='D', ms=10, label='HLTDS')

axis.legend()

axistwin = axis.twiny()
listlabl = []
for s2nrtemp in s2nr:
    listlabl.append('%.3g' % s2nrtemp)
axistwin.set_xticklabels(listlabl)
axistwin.set_xlabel('F087 SNR at 23 AB magnitude')
#figr.colorbar(c, ax=axis)
axis.set_xlabel('Exposure time [seconds]')
axis.set_ylabel('Survey area [degree$^2$]')
print('Writing to %s...' % path)
#plt.subplots_adjust(top=0.9, right=0.8)#tight_layout()
plt.savefig(path)
plt.close()


# construct global object
gdat = tdpy.gdatstrt()
    
# import default Roman Space Telescope configuration
## available filters: 
#listnameband = ['F062', 'F087', 'F106', 'F129', 'F158', 'F184', 'F146', 'F213']
namebandfilt = 'F087'
listnameband = [namebandfilt]
sim_pipeline.Observations.roman_speclite.configure_roman_filters()
path = os.path.dirname(sim_pipeline.__file__)
module_path, _ = os.path.split(path)
skypy_config = os.path.join(module_path, 'data/SkyPy/roman-like.yml')
roman_filters = sim_pipeline.Observations.roman_speclite.filter_names()
speclite.filters.load_filters(roman_filters[0], roman_filters[1], roman_filters[2], roman_filters[3],
                              roman_filters[4], roman_filters[5], roman_filters[6], roman_filters[7])

# cosmology
gdat.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# survey area
gdat.areasurv = 10.
gdat.sky_area = Quantity(value=gdat.areasurv, unit='deg2')

listnamepara = ['velodisp', 'massstel', 'angleins', 'redslens', 'redssour', 'magtsour', 'magtlens']
listlablpara = [r'$\sigma_v$', r'$\log(M_{*})$', r'$\theta_E$', r'$z_{\rm l}$', r'$z_{\rm s}$', r'$m_{\rm source}$', r'$m_{\rm lens}$']
gdat.strgheadpara = ''
for k, labl in enumerate(listlablpara):
    if k != 0:
        gdat.strgheadpara += ', '
    gdat.strgheadpara += labl

dictparaggln = dict()

path = pathdata + 'dictparaggln_Area%08.4g.csv' % gdat.areasurv
if not os.path.exists(path):
    print('Defining galaxies using SkyPyPipeline...')
    pipeline = SkyPyPipeline(skypy_config=skypy_config, sky_area=gdat.sky_area)
    
    print('Lens galaxies...')
    kwargs_deflector_cut = {'band': namebandfilt, 'band_max':23, 'z_min': 0.01, 'z_max': 2.}
    kwargs_mass2light = {}
    lens_type = 'early-type'
    if lens_type == 'early-type':
        from sim_pipeline.Lenses.early_type_lens_galaxies import EarlyTypeLensGalaxies
        gdat.lens_galaxies = EarlyTypeLensGalaxies(pipeline.red_galaxies, kwargs_cut=kwargs_deflector_cut,
                                                    kwargs_mass2light=kwargs_mass2light, cosmo=gdat.cosmo,
                                                    sky_area=gdat.sky_area)
    elif lens_type == 'all-galaxies':
        from sim_pipeline.Lenses.all_lens_galaxies import AllLensGalaxies
        gdat.lens_galaxies = AllLensGalaxies(pipeline.red_galaxies, pipeline.blue_galaxies,
                                              kwargs_cut=kwargs_deflector_cut, kwargs_mass2light=kwargs_mass2light,
                                              cosmo=gdat.cosmo, sky_area=gdat.sky_area)
    else:
        raise ValueError('lens_type %s is not supported' % lens_type)
    numbdelf = gdat.lens_galaxies.deflector_number()
    print('Number of deflectors: %d, scaled to HLWAS (%d sq deg): %d' % (numbdelf, areahlws, int(numbdelf * areahlws / gdat.areasurv)))
    
    print('Source galaxies...')
    kwargs_source_cut = {'band': namebandfilt, 'band_max':24, 'z_min': 0.01, 'z_max': 4.}
    source_type = 'galaxies'
    if source_type == 'galaxies':
        from sim_pipeline.Sources.galaxies import Galaxies
        gdat.source_galaxies = Galaxies(pipeline.blue_galaxies, kwargs_cut=kwargs_source_cut, cosmo=gdat.cosmo)
    else:
        raise ValueError('source_type %s is not supported' % source_type)
    
    print('Sampling GG lenses...')
   
    # Draw a population of galaxy-galaxy lenses within the area.
    # Initialize an empty list to store the GGLens instances
    listobjtgglntotl = []
    listparalens = []
    
    # total number of lens galaxies
    numbgalxdefl = gdat.lens_galaxies.deflector_number()

    # total number of source galaxies
    numbgalxsour = gdat.source_galaxies.galaxies_number()
        
    for kk in range(numbgalxdefl):
        
        deflector = gdat.lens_galaxies.draw_deflector()
        areatest = retr_areatest(deflector['vel_disp'])
        
        # Poisson mean of the number of source galaxies in the test area
        numbgalxsourtestmean = areatest * numbgalxsour / gdat.areasurv / 12960000.
        
        # number of source galaxies in the test area
        numbgalxsourtest = np.random.poisson(lam=numbgalxsourtestmean)
        for n in range(numbgalxsourtest):
            source = gdat.source_galaxies.draw_galaxy()
            if source['z'] > deflector['z']:
                objtggln = GGLens(deflector_dict=deflector, source_dict=source, cosmo=gdat.cosmo, test_area=areatest)
                listobjtgglntotl.append(objtggln)
    
    
    dictparaggln['Candidate'] = dict()
    listnamepara = ['velodisp', 'massstel', 'angleins', 'redssour', 'redslens', 'xposlens', 'yposlens', 'xpossour', 'ypossour', 'numbimag', 'magnsour', 'maxmdistimag']
    for nameband in listnameband:
        listnamepara += ['magtlens%s' % nameband]
        listnamepara += ['magtsour%s' % nameband]
        listnamepara += ['magtsourMagnified%s' % nameband]
    
    for namepara in listnamepara:
        dictparaggln['Candidate'][namepara] = np.empty(len(listobjtgglntotl))
    
    from tqdm import tqdm

    numbobjtgglntotl = len(listobjtgglntotl)
    for k in tqdm(range(numbobjtgglntotl)):
        
        if listobjtgglntotl[k].einstein_radius != listobjtgglntotl[k]._theta_E:
            raise Exception('')

        dictparaggln['Candidate']['velodisp'][k] = listobjtgglntotl[k].deflector_velocity_dispersion()
        dictparaggln['Candidate']['massstel'][k] = listobjtgglntotl[k].deflector_stellar_mass() * 1e-12
        for nameband in listnameband:
            dictparaggln['Candidate']['magtsour%s' % nameband][k] = listobjtgglntotl[k].source_magnitude(band=nameband)
            dictparaggln['Candidate']['magtsourMagnified%s' % nameband][k] = listobjtgglntotl[k].source_magnitude(band=nameband, lensed=True)
            dictparaggln['Candidate']['magtlens%s' % nameband][k] = listobjtgglntotl[k].deflector_magnitude(band=nameband)
        print('dictparaggln[Candidate][magtsourMagnifiedF087]')
        print(dictparaggln['Candidate']['magtsourMagnifiedF087'])
        dictparaggln['Candidate']['angleins'][k] = listobjtgglntotl[k].einstein_radius
        dictparaggln['Candidate']['redssour'][k] = listobjtgglntotl[k].source_redshift
        dictparaggln['Candidate']['redslens'][k] = listobjtgglntotl[k].lens_redshift
        dictparaggln['Candidate']['magnsour'][k] = listobjtgglntotl[k].host_magnification()
        
        # central positions of the lens and the source
        posilens, posisour = listobjtgglntotl[k].position_alignment()
        dictparaggln['Candidate']['xposlens'][k] = posilens[0]
        dictparaggln['Candidate']['yposlens'][k] = posilens[1]
        dictparaggln['Candidate']['xpossour'][k] = posisour[0]
        dictparaggln['Candidate']['ypossour'][k] = posisour[1]
        
        # image positions
        posiimag = listobjtgglntotl[k].get_image_positions()
        dictparaggln['Candidate']['numbimag'][k] = posiimag[0].size
        
        if posiimag[0].size > 1:
            maxmdistimag = np.amax(np.sqrt((posiimag[0][:, None] - posiimag[0][None, :])**2 + (posiimag[1][:, None] - posiimag[1][None, :])**2))
        else:
            maxmdistimag = np.nan
        dictparaggln['Candidate']['maxmdistimag'][k] = maxmdistimag

        #if k == 0:
        #    for nn in range(dictparaggln['Candidate']['numbimag'][k]):
        #        namepara = 'xposimag%d' % nn
        #        dictparaggln['Candidate'][namepara] = np.empty(len(listobjtgglntotl))
        #        listnamepara += [namepara]
        #
        #for nn in range(dictparaggln['Candidate']['numbimag'][k]):
        #    dictparaggln['Candidate']['ypossourimag%d' % nn][k] = posisour[0][0]
    
    print('Writing to %s..' % path)
    pd.DataFrame.from_dict(dictparaggln['Candidate']).to_csv(path, index=False)


    numbyaxi = 2
    numbxaxi = 3
    numbfram = numbyaxi * numbxaxi
    print('Plotting sample images...')
    num_pix = 50
    observatory = 'Roman'
    rgb_band_list = ['F184', 'F129', 'F062']
    add_noise = True
    
    figr, axes = plt.subplots(numbyaxi, numbxaxi, figsize=(numbxaxi * 3, numbyaxi * 3))
    
    k = 0
    kk = 0
    while kk < numbfram and k < dictparaggln['Candidate']['maxmdistimag'].size:
        
        if dictparaggln['Candidate']['maxmdistimag'][k] > 0.4 and dictparaggln['Candidate']['numbimag'][k] > 1 and dictparaggln['Candidate']['angleins'][k] > 0.5 and \
                                            dictparaggln['Candidate']['magtsourMagnified%s' % namebandfilt][k] < 21.:
                                            #dictparaggln['Candidate']['magtlens%s' % namebandfilt][k] > 21. and \
            
            print('dictparaggln[Candidate][maxmdistimag][k]')
            print(dictparaggln['Candidate']['maxmdistimag'][k])
            print('dictparaggln[Candidate][numbimag][k]')
            print(dictparaggln['Candidate']['numbimag'][k])
            print('dictparaggln[Candidate][magtsourMagnified%s % namebandfilt][k]')
            print(dictparaggln['Candidate']['magtsourMagnified%s' % namebandfilt][k])
            print('kk')
            print(kk)
            i = kk // numbyaxi
            j = kk % numbyaxi
            ax = axes[j, i]
            
            #titl = 'maxmdistimag: %g, numbimag'
            #figr, ax = plt.subplots(numbyaxi, numbxaxi, figsize=(numbxaxi * 3, numbyaxi * 3))
            
            # simulate images
            image_r = sim_pipeline.image_simulation.simulate_image(lens_class=listobjtgglntotl[k], band=rgb_band_list[0], num_pix=num_pix, add_noise=add_noise, observatory=observatory)
            image_g = sim_pipeline.image_simulation.simulate_image(lens_class=listobjtgglntotl[k], band=rgb_band_list[1], num_pix=num_pix, add_noise=add_noise, observatory=observatory)
            image_b = sim_pipeline.image_simulation.simulate_image(lens_class=listobjtgglntotl[k], band=rgb_band_list[2], num_pix=num_pix, add_noise=add_noise, observatory=observatory)
        
            image = np.empty(list(image_b.shape) + [3])

            image[:, :, 0] = image_r
            image[:, :, 1] = image_g
            image[:, :, 2] = image_b
            image = np.arcsinh(image)
            image -= np.amin(image)
            image /= np.amax(image)
        
            ax.imshow(image, aspect='equal', origin='lower')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.autoscale(False)
            
            #path = pathvisu + 'grid_%d.png' % k
            #print('Writing to %s...' % path)
            #figr.savefig(path)
            
            kk += 1
        
        k += 1

    figr.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    path = pathvisu + 'grid.png'
    print('Writing to %s...' % path)
    figr.savefig(path)
    








else:
    print('Reading from %s...' % path)
    
    dictparaggln['Candidate'] = pd.read_csv(path).to_dict(orient='list')
    for name in dictparaggln['Candidate'].keys():
        dictparaggln['Candidate'][name] = np.array(dictparaggln['Candidate'][name])

numbggln = dictparaggln['Candidate']['velodisp'].shape[0]

print('Filtering GG lenses with the additional criteria below...')
dictfilt = dict()
dictfilt['Roman-Detectable'] = {'min_image_separation': 0.5, 'mag_arc_limit': {namebandfilt: 20.}}
dictfilt['Roman-Substructure'] = {'min_image_separation': 0.8, 'mag_arc_limit': {namebandfilt: 19}}

# conditions
dictbool = dict()
dictbool['Candidate'] = dictparaggln['Candidate']['numbimag'] > -100
dictbool['All'] = dictparaggln['Candidate']['numbimag'] > 1

#dictbool['sepalar0'] = dictparaggln['Candidate']['angleins'] > dictfiltdete['min_image_separation']
#dictbool['SeparationWithinFOV'] = dictparaggln['Candidate']['angleins'] < 0.5 * dictfiltdete['max_image_separation']

dictbool['DetectableSeparation'] = dictparaggln['Candidate']['maxmdistimag'] > dictfilt['Roman-Detectable']['min_image_separation']
dictbool['SubstructureSeparation'] = dictparaggln['Candidate']['maxmdistimag'] > dictfilt['Roman-Substructure']['min_image_separation']

dictbool['SeparationWithinFOV'] = dictparaggln['Candidate']['maxmdistimag'] < 10.

dictbool['Roman-DetectableBright'] = False
for nameband in dictfilt['Roman-Detectable']['mag_arc_limit'].keys():
    dictbool['Roman-DetectableBright'] = dictparaggln['Candidate']['magtsourMagnified%s' % nameband] < dictfilt['Roman-Detectable']['mag_arc_limit'][nameband]

dictbool['Roman-SubstructureBright'] = False
for nameband in dictfilt['Roman-Substructure']['mag_arc_limit'].keys():
    dictbool['Roman-SubstructureBright'] = dictparaggln['Candidate']['magtsourMagnified%s' % nameband] < dictfilt['Roman-Substructure']['mag_arc_limit'][nameband]

dictbool['Roman-Detectable'] = dictbool['All'] & dictbool['DetectableSeparation'] & dictbool['SeparationWithinFOV'] & dictbool['Roman-DetectableBright']
dictbool['Roman-Substructure'] = dictbool['All'] & dictbool['SubstructureSeparation'] & dictbool['SeparationWithinFOV'] & dictbool['Roman-SubstructureBright']


listnamecond = list(dictbool.keys())


dictindxggln = dict()
for namecond in listnamecond:
    dictindxggln[namecond] = np.where(dictbool[namecond])[0]


## TODO: test for SN ratio in surface brightness


numbgglnobsv = dictindxggln['Roman-Detectable'].size

print('Plotting the lens parameters...')
    
dictnumbsamp = {'Candidate': numbggln}
dictindxsamp = dict()
dictindxsamp['Candidate'] = dict()

nicomedia.retr_subp(dictparaggln, dictnumbsamp, dictindxsamp, 'Candidate', 'All', dictindxggln['All'])
nicomedia.retr_subp(dictparaggln, dictnumbsamp, dictindxsamp, 'Candidate', 'Roman-Detectable', dictindxggln['Roman-Detectable'])
nicomedia.retr_subp(dictparaggln, dictnumbsamp, dictindxsamp, 'Candidate', 'Roman-Substructure', dictindxggln['Roman-Substructure'])
for name in ['Candidate', 'All', 'Roman-Detectable', 'Roman-Substructure']:
    print('Number of %s lenses: %d, scaled to HLWAS (%d deg2): %d, scaled to HLDS (%d deg2): %d' % (name, dictindxggln[name].size, \
                                                                                                    areahlws, int(dictindxggln[name].size * areahlws / gdat.areasurv), \
                                                                                                    areadhls, int(dictindxggln[name].size * areadhls / gdat.areasurv)))

listdictlablcolrpopl = []
listboolcompexcl = []
listtitlcomp = []

if dictindxggln['Roman-Substructure'].size > 0:
    listdictlablcolrpopl.append(dict())
    listdictlablcolrpopl[-1]['All'] = ['All', 'red', 1]
    listdictlablcolrpopl[-1]['Roman-Detectable'] = ['Roman-Detectable', 'blue', 3]
    listdictlablcolrpopl[-1]['Roman-Substructure'] = ['Roman-Characterizable Substructure', 'green', 6]
    listtitlcomp.append('Galaxy-galaxy-type strong lenses')
    listboolcompexcl.append(True)

lablsampgene = 'G-G lens'
pergamon.init('GalaxyGalaxyLensPopulation', \
              dictpopl=dictparaggln, \
              booldiag=True, \
              pathbase=pathbase, \
              listtitlcomp=listtitlcomp, \
              lablsampgene=lablsampgene, \
              listboolcompexcl=listboolcompexcl, \
              listdictlablcolrpopl=listdictlablcolrpopl, \
             )


print('The entire pipeline has run in %g seconds...' % (time.time() - timeinitpipe))

