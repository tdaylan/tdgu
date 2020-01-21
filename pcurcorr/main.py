from tdpy.util import summgene

import pandas as pd

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt

import os

pathbase = os.environ['TDGU_DATA_PATH'] + '/pcurcorr/'
pathdata = pathbase + 'data/'
pathimag = pathbase + 'imag/'

listnameresu = [r'$D_{TESS}$ [ppm]', '$D_{3.6}$ [ppm]', '$D_{4.5}$ [ppm]', '$T_D$ [K]', '$A_g$']

listdata = np.array([\
            #['WASP-36 b'  , 120, 80 , 100 , 914 , 578 , 1953, 544 , 1420, 160, 180,   0., 0.52],\
            ['WASP-43 b'  , 170, 70 , 80  , 3230, 60  , 3830, 80  , 1571, 44 , 43 , 0.13, 0.07],\
            #['WASP-46 b'  , 110, 80 , 110 , 1360, 701 , 4446, 589 , 1870, 130, 130,   0., 0.40],\
            #['WASP-64 b'  , 230, 110, 120 , 2859, 270 , 2071, 471 , 1960, 110, 110,   0., 0.76],\
            ['WASP-77 A b', 55 , 18 , 20  , 2016, 94  , 2487, 127 , 1739, 30 , 31, 0.07, 0.03],\
            #['WASP-78 b'  , 220, 85 , 84  , 2001, 218 , 2013, 351 , 2540, 150, 170,   0., 0.49],\
            ['WASP-100 b' , 101, 13 , 13  , 1267, 98  , 1720, 119 , 2451, 89 , 83 , 0.20, 0.08],\
            #['WASP-18 b'  , 341, 18 , 17  , 3037, 62  , 4033, 97  , 3037, 36 , 36 , 0., 0.04],\
            ['WASP-19 b'  , 101, 106, 131 , 5015,175  , 5343., 318 , 2174, 55 , 55 , 0.19, 0.08],\
            ['WASP-121 b' , 101, 43 , 42  , 3685, 114 , 4684, 121 , 2577, 59 , 63 , 0.24, 0.04],\
            ])
arryresu = listdata[:, np.array([1, 4, 6, 8, 11])].astype(float)

# look into the Exoplanet Archive to get the mass and radii of planets
path = pathdata + 'compositepars_2019.12.12_04.51.12.csv'
objtarch = pd.read_csv(path, skiprows=124)
arryinpt = objtarch.to_numpy()

# finding the indices of the known planets that match the planets of interest
listindx = []
for nameplan in listdata[:, 0]:
    indx = np.where(objtarch['fpl_name'] == nameplan)[0]
    if indx.size == 0:
        print('nameplan')
        print(nameplan)
        raise Exception('')
    listindx.append(indx[0])
listindx = np.array(listindx)
arryinpt = arryinpt[listindx, :]

numbcompresu = arryresu.shape[1]
indxcompresu = np.arange(numbcompresu)
numbplan = arryresu.shape[0]

# find the columns of the known planets that are floats
indxcolsflot = []
for k in range(arryinpt.shape[1]):
    booltemp = False
    for n in range(arryinpt.shape[0]):
        if isinstance(arryinpt[n, k], str):
            booltemp = True
    if not booltemp:
        indxcolsflot.append(k)
indxcolsflot = np.array(indxcolsflot)
arryinpt = arryinpt[:, indxcolsflot]
listnameinpt = objtarch.columns[indxcolsflot]

numbcompinpt = arryinpt.shape[1]
indxcompinpt = np.arange(numbcompinpt)

numbsamp = 1
indxsamp = np.arange(numbsamp)

coef = np.empty((numbsamp, numbcompinpt, numbcompresu))
pval = np.empty((numbsamp, numbcompinpt, numbcompresu))

# find the p value
numbtest = numbcompinpt * numbcompresu
for i in indxsamp:
    arryinpttemp = arryinpt * (np.random.randn() * 0.1 + 1.)
    arryresutemp = arryresu * (np.random.randn() * 0.1 + 1.)
    for n in indxcompinpt: 
        for k in indxcompresu:
            stdvinpt = np.std(arryinpttemp[:, n])
            stdvresu = np.std(arryresutemp[:, k])
            if stdvinpt < 1e-6 or stdvresu < 1e-6:
                coef[i, n, k] = 0.
                pval[i, n, k] = 2.
            else:
                tempinpt = (arryinpttemp[:, n] - np.mean(arryinpttemp[:, n])) / stdvinpt
                tempresu = (arryresutemp[:, k] - np.mean(arryresutemp[:, k])) / stdvresu
                coef[i, n, k], pval[i, n, k] = scipy.stats.pearsonr(tempinpt, tempresu)
pval = np.mean(pval, 0)
coef = np.mean(coef, 0)

# sort with respect to pvalue
indxsort = np.argsort(pval.flatten())
indxtest = np.arange(numbtest)
indxcompinptmesh, indxcompresumesh = np.meshgrid(indxcompinpt, indxcompresu, indexing='ij')
indxinptsort = indxcompinptmesh.flatten()[indxsort]
indxresusort = indxcompresumesh.flatten()[indxsort]

# plot those correlations with the smallest p values
for a in range(indxsort.size):
    titl = 'p = %g, C = %g, %s %s' % (pval.flatten()[indxsort[a]], coef.flatten()[indxsort[a]], listnameinpt[indxinptsort[a]], listnameresu[indxresusort[a]])
    indxxpos = indxinptsort[a]
    indxypos = indxresusort[a]
    if pval.flatten()[indxsort[a]] != 2. and (np.std(arryinpt[:, indxinptsort[a]]) < 1e-6 or np.std(arryresu[:, indxresusort[a]]) < 1e-6):
        raise Exception('')
    figr, axis = plt.subplots(figsize=(6, 6))
    if listnameresu[indxresusort[a]] == '$A_g$':
        yerr = listdata[:, 12].astype(float)
        xerr = None
    else:
        yerr = None
        xerr = None
    axis.errorbar(arryinpt[:, indxinptsort[a]], arryresu[:, indxresusort[a]], yerr=yerr, xerr=xerr, fmt='o')
    axis.set_xlabel(listnameinpt[indxinptsort[a]])
    axis.set_ylabel(listnameresu[indxresusort[a]])
    axis.set_title(titl)
    plt.tight_layout()
    path = pathimag + 'scat_%02d.png' % (a)
    print('Writing to %s...' % path)
    print
    plt.savefig(path)
    plt.close()
    
