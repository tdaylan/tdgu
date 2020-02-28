from tdpy.util import summgene
import tdpy.util
import tdpy.mcmc
import tesstarg.util

import pandas as pd

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt

import os

def retr_modl(gdat, para, inpt):
    
    
    angl = para[0]
    gamm = para[1]

    #slop = para[0]
    #intc = para[1]
    
    slop = -1. / np.tan(angl)
    intc = gamm / np.sin(angl)
    
    line = inpt * slop + intc
    
    return line, []


def retr_llik(gdat, para):
    
    angl = para[0]
    gamm = para[1]

    #slop = para[0]
    #intc = para[1]
    
    #angl = np.arctan(-1. / slop)
    #gamm = np.sin(angl) * intc

    dist = np.cos(angl) * gdat.tempinpt + np.sin(angl) * gdat.tempresu - gamm
    vari = np.cos(angl)**2 * gdat.tempinptstdv**2 + np.sin(angl)**2 * gdat.tempresustdv**2

    llik = -0.5 * np.sum(dist**2 / vari)

    return llik

# paths
pathbase = os.environ['TDGU_DATA_PATH'] + '/pcurcorr/'
pathdata = pathbase + 'data/'
pathimag = pathbase + 'imag/'

# plotting
strgplotextn = 'png'

# construct global object
gdat = tdpy.util.gdatstrt()

lablalbe = '$A_g$'
for strgtype in ['all', 'NoWASP-18b']:

    listdatainpt = np.array([\
                        ['WASP-36 b'  , 0.23, 0.19],\
                        ['WASP-43 b'  , 0.14, 0.07],\
                        ['WASP-46 b'  , 0.12, 0.16],\
                        ['WASP-64 b'  , 0.39, 0.24],\
                        ['WASP-77 A b', 0.06, 0.04],\
                        ['WASP-78 b'  , 0.22, 0.22],\
                        ['WASP-100 b' , 0.26, 0.07],\
                        ['WASP-19 b'  , 0.19, 0.08],\
                        ['WASP-121 b' , 0.27, 0.04],\
                        ])
    if strgtype == 'all':
        listdatainpt = np.vstack([listdatainpt, np.array(['WASP-18 b', 0., 0.02])])
    listalbe = listdatainpt[:, 1].astype(float)
    listalbestdv = listdatainpt[:, 2].astype(float)
        
    listlablinpt = ['$g$', '$T_p$', '[Fe/H]']
    listnameinpt = ['logg', 'ptmp', 'meta']
    
    numbplan = 0
    path = pathdata + 'TESSYear1Parameters.txt'
    print('Reading from %s...' % path)
    objtfile = open(path, 'r')
    for k, line in enumerate(objtfile):
        numbplan += 1
    if strgtype != 'all':
        numbplan -= 1
    
    arryinpt = np.empty((numbplan, 6))
    arryinptstdv = np.empty((numbplan, 6))
    liststrgplaninpt = np.empty(numbplan, dtype=object)
    
    objtfile = open(path, 'r')
    for k, line in enumerate(objtfile):
        if k >= numbplan:
            continue
        linesplt = line.split(',')
        # smax [AU], Planet R [R_J], P M [M_J], St Tempterature [K], St Rad [R_S], St Fe/H
        arryinpt[k, 0] = float(linesplt[0])
        arryinpt[k, 1] = float(linesplt[2])
        arryinpt[k, 2] = float(linesplt[4])
        arryinpt[k, 3] = float(linesplt[6])
        arryinpt[k, 4] = float(linesplt[8])
        arryinpt[k, 5] = float(linesplt[10])
        arryinptstdv[k, 0] = float(linesplt[1])
        arryinptstdv[k, 1] = float(linesplt[3])
        arryinptstdv[k, 2] = float(linesplt[5])
        arryinptstdv[k, 3] = float(linesplt[7])
        arryinptstdv[k, 4] = float(linesplt[9])
        arryinptstdv[k, 5] = float(linesplt[11])
    
    numbcompinpt = 3
    
    indxcompinpt = np.arange(numbcompinpt)
    arryinpttemp = np.empty((numbplan, numbcompinpt))
    arryinptstdvtemp = np.empty((numbplan, numbcompinpt))
    # log of the plenatery surface gravity
    grav = 2479 * arryinpt[:, 2] / arryinpt[:, 1]**2
    stdvgrav = 2479 * np.sqrt(arryinptstdv[:, 2]**2 + 4. * arryinptstdv[:, 1]**2)
    arryinpttemp[:, 0] = np.log10(grav)
    arryinptstdvtemp[:, 0] = 0.434 * stdvgrav / grav
    # plenatery temperature
    arryinpttemp[:, 1] = arryinpt[:, 3] * np.sqrt(arryinpt[:, 4] / arryinpt[:, 0] / 215.)
    arryinptstdvtemp[:, 1] = np.sqrt(arryinptstdv[:, 3]**2 + 0.25 * arryinptstdv[:, 4]**2 + 0.25 * (215. * arryinptstdv[:, 0])*2)
    # stellar metallicity
    arryinpttemp[:, 2] = arryinpt[:, 5]
    arryinptstdvtemp[:, 2] = arryinptstdv[:, 5]
    
    numbsampwalk = 10000
    numbsampburnwalk = 20000
    listlablpara = [[r'$\alpha$', ''], [r'$\rho$', '']]
    listscalpara = ['self', 'self']
    listmeangauspara = None
    liststdvgauspara = None
    listminmpara = [0, -1e1]
    listmaxmpara = [np.pi, 1e1]
    
    numbeval = 2
    numbsampfeww = 1000
    listmodl = np.empty((numbcompinpt, numbsampfeww, numbeval))
    for k in indxcompinpt: 
        gdat.tempinpt = arryinpttemp[:, k]
        gdat.tempresu = listalbe
        gdat.tempinptstdv = arryinptstdvtemp[:, k]
        gdat.tempresustdv = listalbestdv
        
        print('k')
        print(k)
        print('gdat.tempinpt')
        print(gdat.tempinpt)
        print('gdat.tempresu')
        print(gdat.tempresu)
        print('gdat.tempinptstdv')
        print(gdat.tempinptstdv)
        print('gdat.tempresustdv')
        print(gdat.tempresustdv)
        print('')
        numbdata = 2 * gdat.tempinpt.size
        strgextn = strgtype + '_' + listnameinpt[k]
        
        minmxpos = np.amin(gdat.tempinpt) / 1.1
        maxmxpos = np.amax(gdat.tempinpt) * 1.1
        minmypos = np.amin(gdat.tempresu) / 1.1
        maxmypos = np.amax(gdat.tempresu) * 1.1
        
        inpteval = np.linspace(minmxpos, maxmxpos, numbeval)
        parapost = tdpy.mcmc.samp(gdat, pathimag, numbsampwalk, numbsampburnwalk, retr_llik, \
                                        listlablpara, listscalpara, listminmpara, listmaxmpara, listmeangauspara, liststdvgauspara, \
                                            numbdata, strgextn=strgextn, strgplotextn=strgplotextn)
        
        numbsamp = parapost.shape[0]
        indxsamp = np.arange(numbsamp)
        indxsampplot = np.random.choice(indxsamp, size=numbsampfeww)
        for i in np.arange(numbsampfeww):
            listmodl[k, i, :], _ = retr_modl(gdat, parapost[i, :], inpteval)
                
        # plot those correlations with the smallest p values
        figr, axis = plt.subplots(figsize=(4, 4))
        yerr = listalbestdv
        xerr = arryinptstdvtemp[:, k]
        for i in np.arange(numbsampfeww):
            axis.plot(inpteval, listmodl[k, i, :], color='b', alpha=0.03)
        
        axis.errorbar(arryinpttemp[:, k], listalbe, yerr=yerr, xerr=xerr, fmt='o', color='k')
        axis.set_ylim([0., 1.1 * np.amax(gdat.tempresu + gdat.tempresustdv)])
        axis.set_xlabel(listlablinpt[k])
        axis.set_ylabel('$A_g$')


        postslop = -1. / np.tan(parapost[:, 0])
        medislop = np.median(postslop)
        lowrslop = np.median(postslop) - np.percentile(postslop, 16)
        upprslop = np.percentile(postslop, 84) - np.median(postslop)

        axis.set_title('Slope: %g +%g -%g' % (medislop, upprslop, lowrslop))
        plt.tight_layout()
        path = pathimag + 'scat_%s_%s.%s' % (strgtype, listnameinpt[k], strgplotextn)
        print('Writing to %s...' % path)
        print('')
        plt.savefig(path)
        plt.close()
    
