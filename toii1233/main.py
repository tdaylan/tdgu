import pexo.main
import numpy as np
import os
import matplotlib.pyplot as plt
import astropy
import tesstarg.util
import tdpy.util
from tdpy.util import summgene
import scipy.signal
import allesfitter

SMALL_SIZE = 10
MEDIUM_SIZE = 10
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


strgtarg = 'TOI1233'
ticitarg = 260647166
strgmast = 'TIC ' + str(ticitarg)
labltarg = 'HD 108236'

inclprio = np.array([90., 90., 90., 90.])
strgtoii = ['1233.01', '1233.02', '1233.03', '1233.04']

smaxprio = np.array([0.1107, 0.1374, 0.0638, 0.04599]) * 2093 # [R_J]
liststrgplan = ['d', 'e', 'c', 'b']
listlabltoii = ['.01', '.02', '.03', '.04']
indxplan = np.array([3, 2, 0, 1])
gdat = tdpy.util.gdatstrt()
gdat.factmsmj = 1048.
factrsrj = 9.95
gdat.massstar = 0.97 * gdat.factmsmj # [M_J]
radistar = 0.888 * factrsrj # [R_J]

pathbase = os.environ['PEXO_DATA_PATH'] + '/TOI1233/'
pathdata = pathbase + 'data/'
pathimag = pathbase + 'imag/'
pathimagtlss = pathimag + 'tlss/'
os.system('mkdir -p %s' % pathimagtlss)
pathimaginit = pathimag + 'init/'
os.system('mkdir -p %s' % pathimaginit)
pathimagradv = pathimag + 'radv/'
os.system('mkdir -p %s' % pathimagradv)
pathinptalle = pathdata + 'inptalle/'
os.system('mkdir -p %s' % pathinptalle)

listcolrplan = ['red', 'green', 'blue', 'magenta', 'yellow', 'cyan']

strgplotextn = 'pdf'

figrsize = [4., 3.]
figrsizeydob = [8., 3.]

numbplan = 4

# .01 
# .02
# .03
# .04
data = np.array([ \
# ExoFOP
[2458571.335571, 0.001771, 14.175671, 0.000956], \
[2458586.566895, 0.001802, 19.593409, 0.002461], \
[2458572.398315, 0.003184,  6.203183, 0.000652], \
[2458572.111694, 0.002823,  3.795304, 0.000381], \
# ExoFAST
[2458571.3389, 0.0034, 14.1743,  0.0016], \
[2458586.5653, 0.0039, 19.5977,  0.0052], \
[2458572.3972, 0.0023, 6.20368, 0.00057], \
[2458572.1147, 0.0053, 3.79546, 0.00063], \
])
epoc = data[:, 0]
epocstdv = data[:, 1]
peri = data[:, 2]
peristdv = data[:, 3]

# TTV
path = pathdata + '/posterior_toi1233_ttv.csv'
print('Reading from %s...' % path)
objtfile = open(path, 'r')
timettvr = [[] for k in range(4)]
ttvr = [[] for k in range(4)]
stdvttvr = [[] for k in range(4)]
for line in objtfile:
    linesplt = line.split(',')
    if linesplt[0].startswith('tts'):
        indxplantemp = int(linesplt[0][4])
        timettvr[indxplantemp].append(float(linesplt[1]))
    if linesplt[0].startswith('ttvs'):
        indxplantemp = int(linesplt[0][5])
        ttvr[indxplantemp].append(float(linesplt[1]) * 24. * 60.)
        stdvttvr[indxplantemp].append(float(linesplt[2]) * 24. * 60.)

figr, axis = plt.subplots(3, 1, figsize=(4, 5), sharex=True)
for jj, j in enumerate(indxplan[:-1]):
    axis[jj].errorbar(timettvr[jj], ttvr[jj], yerr=stdvttvr[jj], color=listcolrplan[j], marker='o', ls='')
    axis[jj].text(0.1, 0.8, liststrgplan[j], transform=axis[jj].transAxes, color=listcolrplan[j])
axis[-1].set_xlabel('Time [BJD - 2457000]')
axis[1].set_ylabel('Transit Timing Variation [minute]')
#axis[1].yaxis.set_label_coords(-0.1, -0.5)
plt.subplots_adjust(left=0.2, bottom=0.2, hspace=0)
path = pathimaginit + 'ttvr.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()





# WASP
path = pathdata + 'WASP/TOI-1233_WASP_200_ORFG_TAMTFA.lc'
arrylcur = np.loadtxt(path, skiprows=109)
arrylcur[:, 0] += 2450000

# correct HJD to BJD
## location of SAAO where WASP-South data were taken
objtloca = astropy.coordinates.EarthLocation.from_geodetic(20.8105, -32.3783)
objttime = astropy.time.Time(arrylcur[:, 0], format='jd', location=objtloca)

## location of the target
rasc = 186.574063
decl = -51.363052
objtcoor = astropy.coordinates.SkyCoord(rasc, decl, unit="deg")

## convert the HJD to BJD
offsbary = np.array(objttime.light_travel_time(objtcoor, kind='barycentric').jd)
offsheli = np.array(objttime.light_travel_time(objtcoor, kind='heliocentric').jd)
offsdelt = offsbary - offsheli
print('offsdelt')
print(type(offsdelt))
print('arrylcur')
print(type(arrylcur))
arrylcur[:, 0] += offsdelt

print('offsbary 24 * 60')
summgene(offsbary * 24. * 60.)
print('offsheli 24 * 60')
summgene(offsheli * 24. * 60.)
print('offsdelt * 24. * 60.')
summgene(offsdelt * 24. * 60.)

# convert differential mag to relative flux
relf, stdvrelf = tesstarg.util.retr_reflfromdmag(arrylcur[:, 1], arrylcur[:, 2])
arrylcur[:, 1] = relf
arrylcur[:, 2] = stdvrelf

indxsort = np.argsort(arrylcur[:, 0])
arrylcur = arrylcur[indxsort, :]
# temp
diff = (arrylcur[1:, 0] - arrylcur[:-1, 0])
print('diff')
summgene(diff)

# save the light curve
path = pathinptalle + 'WASP.csv'
print('Writing to %s...' % path)
np.savetxt(path, arrylcur, delimiter=',')

# plot the light curve
figr, axis = plt.subplots(figsize=figrsizeydob)
axis.plot(arrylcur[:, 0], arrylcur[:, 1])
axis.set_xlabel('Time [HJD - 2458925]')
axis.set_ylabel('Relative Flux')
#axis.set_ylim([0.997, 1.003])
plt.subplots_adjust(left=0.2, bottom=0.2)
path = pathimaginit + 'lcurwasp.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()

print('arrylcur[:, 0]')
summgene(arrylcur[:, 0])
print('arrylcur[:, 1]')
summgene(arrylcur[:, 1])
print('arrylcur[:, 2]')
summgene(arrylcur[:, 2])
dicttlss = tesstarg.util.exec_tlss(arrylcur, pathimagtlss)
time = arrylcur[:, 0]
flux = arrylcur[:, 1]
stdvflux = arrylcur[:, 2]
logname = 'ejknfksd'
#allesfitter.transit_search.injection_recovery.inject_and_tls_search(time, flux, stdvflux, listperi, listrradiplan, logfname, SNR_threshold=5.)


# RV data

def retr_llik_radv(gdat, para):
    
    offs = para[:3]
    llik = 0.
    massplan = para[3] / gdat.factmjme # [M_J]
    epoc = para[4]
    peri = para[5]
    for k, strginst in enumerate(gdat.liststrginst):
        radvmodl = tesstarg.util.retr_radv(gdat.listdata[k][:, 0], epoc, peri, massplan, gdat.massstar/gdat.factmsmj, 90., 0., 0.)
        llik += -0.5 * np.sum((gdat.listdata[k][:, 1] - radvmodl - para[k])**2 / gdat.listdata[k][:, 2]**2)
    
    return llik


#gdat.liststrginst = ['CHIRON/CHIRON.csv', 'PFS', 'NRES']
gdat.liststrginst = ['CHIRON', 'PFS', 'NRES']
gdat.listdata = []
numbdata = 0
for a, strginst in enumerate(gdat.liststrginst):
    path = pathdata + 'TFOP/' + strginst + '/' + strginst + '.csv'
    data = np.loadtxt(path, skiprows=1, delimiter=',')
    if not strginst == 'PFS':
        data[:, 1:3] *= 1e3
    numbdata += data.shape[0]
    rmsq = np.std(data[:, 1])
    
    print('strginst')
    print(strginst)
    print('rmsq')
    print(rmsq)
    print('')

    gdat.listdata.append(data)

    # save the PFS data for allesfitter
    if strginst == 'PFS':
        path = pathinptalle + 'PFS.csv'
        print('Writing to %s...' % path)
        data[:, 1:3] *= 1e-3 # [km/s]
        np.savetxt(path, data, delimiter=',')

maxmtime = -1e12
minmtime = 1e12
for a, strg in enumerate(gdat.liststrginst):
    minmtime = min(np.amin(gdat.listdata[a][:, 0]), minmtime)
    maxmtime = max(np.amax(gdat.listdata[a][:, 0]), maxmtime)

listlablpara = [['$O_{C}$', 'm s$^{-1}$'], ['$O_{P}$', 'm s$^{-1}$'], ['$O_{N}$', 'm s$^{-1}$']]
listminmpara = np.array([ 1e4,-1e1, 1e4])
listmaxmpara = np.array([ 2e4, 1e1, 2e4])
listlablpara += [['$M$', '$M_E$'], ['$T_{0}$', 'BJD'], ['$P$', 'days']]
listminmpara = np.concatenate([listminmpara, np.array([ 0., minmtime,  0.1])])
listmaxmpara = np.concatenate([listmaxmpara, np.array([1e4, maxmtime, 400.])])
listmeangauspara = None
liststdvgauspara = None
numbpara = len(listlablpara)
indxpara = np.arange(numbpara)
listscalpara = ['self' for k in indxpara]

numbsampwalk = 500
numbsampburnwalk = 0
numbsampburnwalkseco = 400
strgextn = 'radv'

listcolrinst = ['cyan', 'brown', 'olive']
## plot data only
figr, axis = plt.subplots(figsize=figrsizeydob)
for a, strg in enumerate(gdat.liststrginst):
    axis.errorbar(gdat.listdata[a][:, 0] - 2458000, gdat.listdata[a][:, 1], yerr=gdat.listdata[a][:, 2], ls='', marker='o', color=listcolrinst[a], label=strg)
axis.set_xlabel('Time [BJD - 2458000]')
axis.set_ylabel('Radial velocity [m s$^{-1}$]')
axis.legend()
plt.subplots_adjust(left=0.2, bottom=0.2)
path = pathimagradv + 'radvdata.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()

## plot Lomb Scargle periodogram
figr, axis = plt.subplots(figsize=(12, 4))
axis.set_ylabel('Power')
axis.set_xlabel('Frequency [1/day]')
arryfreq = np.linspace(0.1, 10., 2000)
for a in range(3):
    if a == 1:
        ydat = scipy.signal.lombscargle(gdat.listdata[a][:, 0], gdat.listdata[a][:, 1], arryfreq)
        axis.plot(arryfreq * 2. * np.pi, ydat, color='black')
#for a in range(4):
#    axis.axvline(a / peri, ls='--', color='black')
plt.tight_layout()
path = pathimagradv + 'lspd.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()


gdat.factmjme = 317.907
listpara = tdpy.mcmc.samp(gdat, pathimagradv, numbsampwalk, numbsampburnwalk, numbsampburnwalkseco, retr_llik_radv, \
             listlablpara, listscalpara, listminmpara, listmaxmpara, listmeangauspara, liststdvgauspara, numbdata, strgextn=strgextn)

numbsamp = listpara.shape[0]
indxsamp = np.arange(numbsamp)


def plot_compsampresi(xdat, listydatdata, lablydat, lablxdat, pathimag, strgplot, listlabldata=None, strgplotextn='pdf'):
    
    numbdata = len(listydatdata)
    indxdata = np.arange(numbdata) 
    figr, axis = plt.subplots(figsize=figrsizeydob)
    for a in range(3):
        # upper data and samples
        
        ## data
        for b in indxdata:
            axis[0].errorbar(xdat, listydatdata[b], yerr=gdat.listdata[a][:, 2], ls='', marker='o', color=listcolrinst[a])
        
        ## 
        # upper data and samples
        
        axis[a].set_ylabel(lablydat)
    axis[2].set_xlabel(lablxdat)
    #plt.subplots_adjust(left=0.2, bottom=0.2)
    path = pathimag + '%s.%s' % (strgplot, strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()


figr, axis = plt.subplots(figsize=figrsizeydob)
for b, strg in enumerate(gdat.liststrginst):
    axis.errorbar(gdat.listdata[b][:, 0] - 2458000, gdat.listdata[b][:, 1] - np.median(listpara[:, b]), yerr=gdat.listdata[b][:, 2], \
                                                                                                    ls='', marker='o', color=listcolrinst[b])

timefine = np.linspace(minmtime, maxmtime, 1000)
### sample model phas
numbsampplot = 1
indxsampplot = np.random.choice(indxsamp, numbsampplot, replace=False)

for nn, n in enumerate(indxsampplot):
    massplan = listpara[n, 0] / gdat.factmjme # [M_J]
    epocplan = listpara[n, 1]
    periplan = listpara[n, 2]
    radv = tesstarg.util.retr_radv(timefine, epocplan, periplan, massplan, gdat.massstar/gdat.factmsmj, 90., 0., 0.)
    axis.plot(timefine - 2458000, radv, alpha=0.1, color='b')

axins = axis.inset_axes([0.5, 0.5, 0.47, 0.47])
for b, strg in enumerate(gdat.liststrginst):
    axis.errorbar(gdat.listdata[b][:, 0] - 2458000, gdat.listdata[b][:, 1] - np.median(listpara[:, b]), yerr=gdat.listdata[b][:, 2], \
                                                                                                    ls='', marker='o', color=listcolrinst[b])
# sub region of the original image
x1, x2, y1, y2 = 630, 730, -10, 10
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels('')
axins.set_yticklabels('')
axis.indicate_inset_zoom(axins)

axis.set_xlabel('Time [BJD - 2458000]')
axis.set_ylabel('Radial velocity [m s$^{-1}$]')
plt.subplots_adjust(left=0.2, bottom=0.2)
path = pathimagradv + 'radvsamp_%d.%s' % (a, strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()


# transit timee uncertainties
numbtran = np.empty(numbplan, dtype=int)
timefutu = 2459000
for j in indxplan:
    numbtran[j] = (timefutu - epoc[j]) / peri[j]
time = [np.empty(numbtran[j]) for j in indxplan]
labltime = [np.empty(numbtran[j], dtype=object) for j in indxplan]
stdvtimetran = [np.empty(numbtran[j]) for j in indxplan]
timeexof = [np.empty(numbtran[j]) for j in indxplan]
for a in range(2):
    for j in indxplan:
        indxtran = np.arange(numbtran[j])
    
        timeexof[j] = epoc[4*a+j] + indxtran * peri[4*a+j]
        
        objttime = astropy.time.Time(timeexof[j], format='jd', scale='utc', out_subfmt='date_hm')
        
        labltime[j] = objttime.iso
        stdvtimetran[j] = np.sqrt(epocstdv[4*a+j]**2 + (indxtran * peristdv[4*a+j])**2) * 24 # [hr]
    
    figr, axis = plt.subplots(figsize=figrsize)
    for j, strgplan in enumerate(liststrgplan):
        axis.plot(timeexof[j], stdvtimetran[j], color=listcolrplan[j], label=listlabltoii[j])
        xlim = axis.get_xlim()
        timeaxis = np.linspace(xlim[0], xlim[1], 7)
        
        objttime = astropy.time.Time(timeaxis, format='jd', scale='utc', out_subfmt='date_hm')
        labltimeaxis = [labltimeaxistem[:10] for labltimeaxistem in objttime.iso]
        axis.set_xticks(timeaxis)
        axis.set_xticklabels(labltimeaxis, rotation=20)
    axis.set_ylabel('$\sigma$ [hour]')
    axis.legend()
    plt.subplots_adjust()
    path = pathimag + 'stdvtimetran_%d.%s' % (a, strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()


# plot Stephen Kane's output
#[]Time (years)         a          e          i          mass          long       node         M    
listlabl = ['Time [year]', '$a$ [AU]', 'Eccentricity', 'Inclination [degree]', 'Mass [$M_E$]', 'Longitude of ascending node [degree]', 'node', 'Mean anomaly']
liststrg = ['time', 'smax', 'ecce', 'incl', 'mass', 'long', 'node', 'M']
numbcols = len(liststrg)
indxcols = np.arange(numbcols)

numbplan = len(strgtoii)
for j, strgplan in enumerate(liststrgplan):
    path = pathdata + 'stability/planet%d.aei' % (j + 1)
    data = np.loadtxt(path, skiprows=4)
    
    if j == 0:
        numbtime = data.shape[0]
        listecce = np.empty((numbtime, numbplan))
    
    #time = data[:, 0]
    #for k in indxcols[1:]:
    #    figr, axis = plt.subplots(figsize=figrsize)
    #    axis.plot(time, data[:, k])
    #    axis.set_xlabel('Time [year]')
    #    axis.set_ylabel(listlabl[k])
    #    plt.subplots_adjust()
    #    path = pathimaginit + '%s_%d.%s' % (liststrg[k], j, strgplotextn)
    #    print('Writing to %s...' % path)
    #    plt.savefig(path)
    #    plt.close()
    
    listecce[:, j] = data[:, 2]

figr, axis = plt.subplots(1, numbplan, figsize=figrsizeydob, sharey=True)
bins = np.linspace(0., np.amax(listecce), 100)
for jj, j in enumerate(indxplan):
    axis[jj].hist(listecce[:, j], color=listcolrplan[j])
    axis[jj].set_xlabel('$e_{%s}$' % liststrgplan[j])
axis[0].set_ylabel('Number of samples')
plt.subplots_adjust(bottom=0.2)
path = pathimaginit + 'histecce.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()


# plot Keivan Stassun's model
path = pathdata + 'stellar_characterization/specstardata.dat'
arry = np.loadtxt(path, skiprows=2)
wlendata = arry[:, 0]
stdvwlendata = arry[:, 1]
fluxdata = arry[:, 2]
stdvfluxdata = arry[:, 3]
path = pathdata + 'stellar_characterization/specstarmodl.dat'
arry = np.loadtxt(path, skiprows=2)
wlenmodl = arry[:, 0]
fluxmodl = 10**arry[:, 1]
figr, axis = plt.subplots(figsize=figrsize)
axis.plot(wlenmodl, fluxmodl, color='b', lw=0.5)
axis.errorbar(wlendata, fluxdata, xerr=stdvwlendata, yerr=stdvfluxdata, ls='', marker='o', color='k', ms=1, lw=0.5)
axis.set_xlabel('Wavelength [$\mu$m]')
axis.set_ylabel('$\lambda F_\lambda$ [erg s$^{-1}$ cm$^{-2}$]')
axis.set_yscale('log')
axis.set_xscale('log')
axis.set_ylim([1e-13, None])
plt.subplots_adjust(left=0.2, bottom=0.2)
path = pathimaginit + 'specstar.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()

# plot LCO/NRES
path = pathdata + 'TFOP/NRES/lscnrs01-fa09-20190612-0037-e91/lscnrs01-fa09-20190612-0037-e91.fits.fz'
listhdun = astropy.io.fits.open(path)
listhdun.info()
for k in range(len(listhdun)):
    #print(repr(listhdun[k].header))
    data = listhdun[k].data
    continue
    if data is None:
        continue
    figr, axis = plt.subplots(figsize=figrsize)
    if data.ndim == 2:
        print('HEEEY')
        print('HEEEY')
        print('HEEEY')
        print('HEEEY')
        print('HEEEY')
        summgene(data)
        axis.imshow(data)
    else:
        pass
        #axis.errorbar(wlendata, fluxdata, xerr=stdvwlendata, yerr=stdvfluxdata, ls='', marker='o', color='k', ms=1)
    #axis.plot(wlenmodl, fluxmodl, color='b')
    #axis.set_xlabel('Wavelenth [$\mu$m]')
    #axis.set_ylabel('Flux [erg s$^{-1}$ cm$^{-2}$)]')
    #axis.set_yscale('log')
    #axis.set_ylim([1e-13, None])
    #plt.subplots_adjust(left=0.2, bottom=0.2)
    path = pathimaginit + 'specnres%d.%s' % (k, strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()

# plot HD 108236 d follow-up detections
dictarry = {}

## MEarth
timeoffs = 2458925
for strg in ['11', '12', '13', '15', '16', '17']:
    strgarry = 'tl' + strg
    path = pathdata + 'TFOP/20200317_MEarth-South/TIC260647166-01_20200317_MEarth-South_defocus_lcs/tel%s/T1_ap2.txt' % strg
    print('Reading from %s...' % path)
    data = np.loadtxt(path, skiprows=17)
    numbtime = data.shape[0]
    dictarry[strgarry] = np.empty((numbtime, 3))
    dictarry[strgarry][:, 0] = data[:, 0]
    dictarry[strgarry][:, 1] = data[:, 18]
    dictarry[strgarry][:, 2] = data[:, 2]
    
    relf, stdvrelf = tesstarg.util.retr_reflfromdmag(dictarry[strgarry][:, 1], dictarry[strgarry][:, 2])
    dictarry[strgarry][:, 1] = relf
    dictarry[strgarry][:, 2] = stdvrelf

# merge data from different defocused telescopes
dictarry['tl00'] = []
for strg in ['11', '12', '13', '15', '16', '17']:
    strgarry = 'tl' + strg
    dictarry['tl00'].append(dictarry[strgarry])
dictarry['tl00'] = np.concatenate(dictarry['tl00'], 0)
## sort the data in time
indx = np.argsort(dictarry['tl00'][:, 0])
dictarry['tl00'] = dictarry['tl00'][indx, :]

# save the MEarth-South light curve
path = pathinptalle + 'MEARTHS.csv'
print('Writing to %s...' % path)
np.savetxt(path, dictarry['tl00'], delimiter=',')

## LCO-CTIO
path = pathdata + 'TFOP/20200317_LCO-CTIO/TIC260647166-01_20200317_LCO-CTIO-1m0_zs_bjd-flux-err-fwhm-humid.dat'
print('Reading from %s...' % path)
data = np.loadtxt(path, skiprows=1)
numbtime = data.shape[0]
dictarry['lcocraww'] = np.empty((numbtime, 3))
dictarry['lcocdetr'] = np.empty((numbtime, 3))
dictarry['lcocraww'][:, 0] = data[:, 0]
dictarry['lcocraww'][:, 1] = data[:, 1]
dictarry['lcocraww'][:, 2] = data[:, 2]
dictarry['lcocdetr'][:, 0] = data[:, 0]
dictarry['lcocdetr'][:, 1] = data[:, 3]
dictarry['lcocdetr'][:, 2] = data[:, 4]

# save the light curve
path = pathinptalle + 'LCO.csv'
print('Writing to %s...' % path)
np.savetxt(path, dictarry['lcocdetr'], delimiter=',')

liststrg = list(dictarry.keys())

numbbins = 40
for strgarry in liststrg:
    strgarrybind = strgarry + 'bind'
    dictarry[strgarrybind] = tesstarg.util.rebn_lcur(dictarry[strgarry], numbbins)
    
    figr, axis = plt.subplots(figsize=figrsize)
    axis.errorbar(dictarry[strgarry][:, 0] - timeoffs, dictarry[strgarry][:, 1], ls='', marker='o', color='grey', ms=1)
    axis.errorbar(dictarry[strgarrybind][:, 0] - timeoffs, dictarry[strgarrybind][:, 1], \
                                                            yerr=dictarry[strgarrybind][:, 2], ls='', marker='o', color='red', ms=3)
    axis.set_xlabel('Time [BJD - 2458925]')
    axis.set_ylabel('Relative Flux')
    plt.subplots_adjust(left=0.2, bottom=0.2)
    path = pathimaginit + 'lcurtfop_%s.%s' % (strgarry, strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()

# plot the two light curves together
figr, axis = plt.subplots(2, 1, figsize=figrsizeydob, sharex=True)
axis[0].errorbar(dictarry['lcocdetr'][:, 0] - timeoffs, dictarry['lcocdetr'][:, 1], ls='', marker='o', color='grey', ms=1)
axis[0].errorbar(dictarry['lcocdetrbind'][:, 0] - timeoffs, dictarry['lcocdetrbind'][:, 1], \
                                                                    yerr=dictarry['lcocdetrbind'][:, 2], ls='', marker='o', color='red', ms=3)
axis[1].errorbar(dictarry['tl00'][:, 0] - timeoffs, dictarry['tl00'][:, 1], ls='', marker='o', color='grey', ms=1)
axis[1].errorbar(dictarry['tl00bind'][:, 0] - timeoffs, dictarry['tl00bind'][:, 1], yerr=dictarry['tl00bind'][:, 2], ls='', marker='o', color='red', ms=3)
indx = np.argmin(abs(timeexof[0] - 2458925))
axis[0].axvline(timeexof[0][indx] - timeoffs, ls='--', alpha=0.8, color='k')
axis[1].axvline(timeexof[0][indx] - timeoffs, ls='--', alpha=0.8, color='k')
axis[1].set_xlabel('Time [BJD - 2458925]')
axis[0].set_ylabel('Relative Flux')
axis[0].yaxis.set_label_coords(-0.07, -0.2)
axis[0].set_ylim([0.997, 1.003])
axis[1].set_ylim([0.997, 1.003])
plt.subplots_adjust(left=0.2, bottom=0.2)
path = pathimaginit + 'lcurtfop.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()


# injection recovery


pexo.main.main( \
                   strgtarg=strgtarg, \
                   labltarg=labltarg, \
                   
                   strgmast=strgmast, \
                   #ticitarg=ticitarg, \
                
                   booltlss=False, \
                   smaxprio=smaxprio, \
                   
                   boolmaskqual=False, \

                   toiitarg=1233, \
                   priotype='exof', \

                   liststrgplan=liststrgplan, \

                   inclprio=inclprio, \

                   massstar=gdat.massstar, \
                   radistar=radistar, \
                    
                  )

