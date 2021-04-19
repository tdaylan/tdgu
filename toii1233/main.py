import pexo.main
import numpy as np
import os
import matplotlib.pyplot as plt
import astropy
import ephesus.util
import tdpy.util
from tdpy.util import summgene
import scipy.signal
import allesfitter
from allesfitter.v2.classes import allesclass2
from allesfitter.v2.translator import translate




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

# Boolean flag to make plots before pexo
boolplotinit = False

strgtarg = 'TOI1233'
labltarg = 'HD 108236'
liststrgplan = ['d', 'e', 'c', 'b']
    
gdat = tdpy.util.gdatstrt()

gdat.factrsrj, gdat.factmsmj, gdat.factrjre, gdat.factmjme, gdat.factaurj = ephesus.util.retr_factconv()
    
# impose the priors from Keivan's SED fit
gdat.radistar = 0.888 * gdat.factrsrj # [R_J]
gdat.stdvradistar = 0.017 * gdat.factrsrj # [R_J]
gdat.massstar = 0.97 * gdat.factmsmj # [M_J]
gdat.stdvmassstar = 0.06 * gdat.factmsmj # [M_J]
gdat.tmptstar = 5730. # [K]
gdat.stdvtmptstar = 50. # [K]

#rvsaprio = np.array([1.97, 2.36, 1.93, 1.30]) * 1e-3 # [km/s]

numbpred = 3
indxpred = np.arange(numbpred)
    
# line styles for different ephemerides
liststylline = ['-', '--', '-.']
data = [ \
# ExoFOP
np.array([ \
 [2458571.335571, 0.001771, 14.175671, 0.000956, np.sqrt(917.111206e-6)], \
 [2458586.566895, 0.001802, 19.593409, 0.002461, np.sqrt(1214.447632e-6)], \
 [2458572.398315, 0.003184,  6.203183, 0.000652, np.sqrt(481.47464e-6)], \
 [2458572.111694, 0.002823,  3.795304, 0.000381, np.sqrt(321.974854e-6)]]), \
# ExoFAST
np.array([ \
 [2458571.3389, 0.0034, 14.1743,  0.0016], \
 [2458586.5653, 0.0039, 19.5977,  0.0052], \
 [2458572.3972, 0.0023, 6.20368, 0.00057], \
 [2458572.1147, 0.0053, 3.79546, 0.00063]]), \
# mine (eccentric)
np.array([ \
 [2458571.3368, 0.0015, 14.17555, 0.0011], \
 [2458586.5677, 0.0014,  19.5917, 0.0022], \
 [2458572.3949, 0.0025,  6.20370, 0.00064], \
 [2458572.1128, 0.0036,  3.79523, 0.00047]]), \
]
rrat = [[] for a in indxpred]
epoc = [[] for a in indxpred]
peri = [[] for a in indxpred]
epocstdv = [[] for a in indxpred]
peristdv = [[] for a in indxpred]
for a in indxpred:
    epoc[a] = data[a][:, 0]
    epocstdv[a] = data[a][:, 1]
    peri[a] = data[a][:, 2]
    peristdv[a] = data[a][:, 3]
    if a == 0:
        rrat[a] = data[a][:, 4]

indxsorttcee = np.argsort(peri[0])
print('indxsorttcee')
summgene(indxsorttcee)
    
if boolplotinit:
    
    strgtoii = ['1233.01', '1233.02', '1233.03', '1233.04']
    
    listlabltoii = ['.01', '.02', '.03', '.04']
    indxplan = np.array([3, 2, 0, 1])
    
    pathbase = os.environ['PEXO_DATA_PATH'] + '/TOI1233/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    pathimagtlsq = pathimag + 'tlsq/'
    os.system('mkdir -p %s' % pathimagtlsq)
    pathimaginit = pathimag + 'init/'
    os.system('mkdir -p %s' % pathimaginit)
    pathimagrvel = pathimaginit + 'rvel/'
    os.system('mkdir -p %s' % pathimagrvel)
    pathinptalle = pathdata + 'inptalle/'
    os.system('mkdir -p %s' % pathinptalle)
    
    listcolrplan = ['red', 'green', 'orange', 'magenta', 'yellow', 'cyan']
    
    strgplotextn = 'pdf'
    
    figrsize = [4., 3.]
    figrsizeydob = [8., 3.]
    
    numbplan = 4
    
    # .01 
    # .02
    # .03
    # .04

    # get transit model based on TESS ephemerides
    #rr = 0.1
    #rsuma = 0.1
    #epoch = 0
    #period = 1
    #cosi = 0
    #ld_law = 'quad'                     #limb darkening law
    #q1 = 0.5                            #transformed ld q1
    #q2 = 0.2                            #transformed ld q2
    #r_host = 1                          #in Rsun
    #m_host = 1                          #in Msun
    #time = np.linspace(-0.5,0.5,1001)   #in days
    #alles = allesclass2()
    #alles.settings = {'companions_phot':['b'], 'inst_phot':['telescope'], 'host_ld_law_telescope':ld_law}
    #alles.params = {'b_rr':rr, 'b_rsuma':rsuma, 'b_epoch':epoch, 'b_period':period, 'b_cosi':cosi, 'host_ldc_q1_telescope':q1, 'host_ldc_q2_telescope':q2}
    #alles.params_host = {'R_host':r_host, 'M_host':m_host}
    #alles.fill()
    #model_flux = alles.generate_model(time, inst='telescope', key='flux')
    
    
    
    # plot TTV
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
    arrylcur[:, 0] += offsdelt
    
    # convert differential mag to relative flux
    rflx, stdvrflx = ephesus.util.retr_rflxfromdmag(arrylcur[:, 1], arrylcur[:, 2])
    arrylcur[:, 1] = rflx
    arrylcur[:, 2] = stdvrflx
    
    indxsort = np.argsort(arrylcur[:, 0])
    arrylcur = arrylcur[indxsort, :]
    
    # save the light curve
    path = pathinptalle + 'WASP.csv'
    print('Writing to %s...' % path)
    np.savetxt(path, arrylcur, delimiter=',')
    
    timeoffstess = 2457000
    timeoffsfolw = 2458000
    timeoffslcoc = 2458925
    timeoffswasp = 2450000
    
    # plot the light curve
    figr, axis = plt.subplots(figsize=figrsizeydob)
    axis.plot(arrylcur[:, 0] - timeoffswasp, arrylcur[:, 1])
    axis.set_xlabel('Time [HJD - %d]' % timeoffswasp)
    axis.set_ylabel('Relative Flux')
    #axis.set_ylim([0.997, 1.003])
    plt.subplots_adjust(left=0.2, bottom=0.2)
    path = pathimaginit + 'lcurwasp.%s' % (strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()
    
    dicttlsq = ephesus.util.exec_tlsq(arrylcur, pathimagtlsq)
    time = arrylcur[:, 0]
    flux = arrylcur[:, 1]
    stdvflux = arrylcur[:, 2]
    logname = 'ejknfksd'
    #allesfitter.transit_search.injection_recovery.inject_and_tls_search(time, flux, stdvflux, listperi, listrradiplan, logfname, SNR_threshold=5.)
    
    
    # RV data
    
    # convert PFS data
    path = pathdata + 'TFOP/PFS/HD108236_PFS.vels'
    np.genfromtxt(path, delimiter=",")
    
    offsbary = np.array(objttime.light_travel_time(objtcoor, kind='barycentric').jd)
    objtloca = astropy.coordinates.EarthLocation.from_geodetic(20.8105, -32.3783)
    objttime = astropy.time.Time(arrylcur[:, 0], format='jd', location=objtloca)
    
    
    
    def retr_llik_rvelsimp(para, gdat):
        
        massplan = para[0] / gdat.factmjme # [M_J]
        sema = retr_rvelsema(1000., massplan, gdat.massstar/gdat.factmsmj, 90., 0.)
        llik = -0.5 * np.sum((sema / gdat.stdvsema)**2)
        
        return llik
    
    
    def retr_llik_rvel(para, gdat):
        
        offs = para[:3]
        llik = 0.
        massplan = para[3] / gdat.factmjme # [M_J]
        epoctemp = para[4]
        peritemp = para[5]
        for k, strginst in enumerate(gdat.liststrginst):
            if strginst == 'PFS':
                rvelmodl = ephesus.util.retr_rvel(gdat.listdata[k][:, 0], epoctemp, peritemp, massplan, gdat.massstar/gdat.factmsmj, 90., 0., 0.)
            else:
                rvelmodl = 0.
            llik += -0.5 * np.sum((gdat.listdata[k][:, 1] - rvelmodl - para[k])**2 / gdat.listdata[k][:, 2]**2)
        
        return llik
    
    
    gdat.liststrginst = ['CHIRON', 'PFS', 'NRES']
    gdat.listdata = []
    numbdata = 0
    rveloffs = []
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
            datatemp = np.copy(data)
            datatemp[:, 1:3] *= 1e-3 # [km/s]
            np.savetxt(path, datatemp, delimiter=',')
    
        rveloffs.append(data[:, 1] - np.mean(data[:, 1]))
    rveloffs = np.concatenate(rveloffs)
    
    gdat.stdvsema = np.std(rveloffs)
    print('gdat.stdvsema')
    print(gdat.stdvsema)
    massplan = 1. / gdat.factmjme # [M_J]
    sema = ephesus.util.retr_rvelsema(1000., massplan, gdat.massstar/gdat.factmsmj, 90., 0.)
    massplan1sig = gdat.stdvsema / sema
    print('massplan1sig')
    print(massplan1sig)
    
    maxmtime = -1e12
    minmtime = 1e12
    for a, strg in enumerate(gdat.liststrginst):
        if strg == 'PFS':
            minmtime = min(np.amin(gdat.listdata[a][:, 0]), minmtime)
            maxmtime = max(np.amax(gdat.listdata[a][:, 0]), maxmtime)
    
    listlablpara = [['$O_{C}$', 'm s$^{-1}$'], ['$O_{P}$', 'm s$^{-1}$'], ['$O_{N}$', 'm s$^{-1}$']]
    listminmpara = np.array([ 14e3,-1e0, 16e3])
    listmaxmpara = np.array([ 16e3, 1e0, 18e3])
    listlablpara += [['$M$', '$M_E$'], ['$T_{0}$', 'BJD'], ['$P$', 'days']]
    listminmpara = np.concatenate([listminmpara, np.array([ 10., minmtime,  50.])])
    listmaxmpara = np.concatenate([listmaxmpara, np.array([1e4, maxmtime, 200.])])
    listmeangauspara = None
    liststdvgauspara = None
    numbpara = len(listlablpara)
    indxpara = np.arange(numbpara)
    listscalpara = ['self' for k in indxpara]
    
    numbsampwalk = 30000
    numbsampburnwalk = 0
    numbsampburnwalkseco = 4000
    strgextn = 'rvel'
    
    listcolrinst = ['cyan', 'brown', 'olive']
    ## plot data only
    figr, axis = plt.subplots(figsize=figrsizeydob)
    for a, strg in enumerate(gdat.liststrginst):
        axis.errorbar(gdat.listdata[a][:, 0] - timeoffsfolw, gdat.listdata[a][:, 1], yerr=gdat.listdata[a][:, 2], ms=1, \
                                                                        ls='', marker='o', color=listcolrinst[a], label=strg)
    axis.set_xlabel('Time [BJD - %d]' % timeoffsfolw)
    axis.set_ylabel('Differential radial velocity [m s$^{-1}$]')
    axis.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    path = pathimagrvel + 'rveldata.%s' % (strgplotextn)
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
    path = pathimagrvel + 'lspd.%s' % (strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()
    
    
    listpara, _ = tdpy.mcmc.samp(gdat, pathimagrvel, numbsampwalk, numbsampburnwalk, numbsampburnwalkseco, retr_llik_rvel, \
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
    
    axisinse = axis.inset_axes([0.5, 0.13, 0.4, 0.4])
    
    pathalle = pathbase + 'allesfits/allesfit_orbt_rvel/'
    print('Reading from %s...' % pathalle)
    objtallervel = allesfitter.allesclass(pathalle)
    timefine = np.linspace(minmtime, maxmtime, 1000)
    ### sample model phas
    numbsampplot = 100
    indxsampplot = np.random.choice(indxsamp, numbsampplot, replace=False)
    for axistemp in [axis, axisinse]:
        for b, strg in enumerate(gdat.liststrginst):
            axistemp.errorbar(gdat.listdata[b][:, 0] - timeoffsfolw, gdat.listdata[b][:, 1] - np.median(listpara[:, b]), yerr=gdat.listdata[b][:, 2], \
                                                                                              ms=2, label=strg, ls='', marker='o', color=listcolrinst[b])
        
        #for nn, n in enumerate(indxsampplot):
        #    massplan = listpara[n, 3] / gdat.factmjme # [M_J]
        #    epocplan = listpara[n, 4]
        #    periplan = listpara[n, 5]
        #    
        #    rvelmodl = objtallervel.get_posterior_median_model('PFS', 'rv', xx=timefine)
        #    #rvel = ephesus.util.retr_rvel(timefine, epocplan, periplan, massplan, gdat.massstar/gdat.factmsmj, 90., 0., 0.)
        #    axistemp.plot(timefine - timeoffsfolw, rvelmodl, alpha=0.1, color='b')
    
    axis.legend(loc='upper center', bbox_to_anchor=(0.6, 0.95), ncol=3, fancybox=True, shadow=True)
        
    # sub region of the original image
    x1, x2, y1, y2 = 170, 220, -10, 10
    axisinse.set_xlim(x1, x2)
    axisinse.set_ylim(y1, y2)
    axis.indicate_inset_zoom(axisinse)
    
    axis.set_xlabel('Time [BJD - %d]' % timeoffsfolw)
    axis.set_ylabel('Differential radial velocity [m s$^{-1}$]')
    plt.subplots_adjust(left=0.2, bottom=0.2)
    path = pathimagrvel + 'rvel.%s' % (strgplotextn)
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()
    
    # transit timee uncertainties
    #time = [np.empty(numbtran[j]) for j in indxplan]
    #labltime = [np.empty(numbtran[j], dtype=object) for j in indxplan]
    #stdvtimetran = [np.empty(numbtran[j]) for j in indxplan]
    #for a in indxpred:
    #    for j in indxplan:
    #        objttime = astropy.time.Time(timetranpred[a][j], format='jd', scale='utc', out_subfmt='date_hm')
    #        labltime[j] = objttime.iso
    #        stdvtimetran[j] = np.sqrt(epocstdv[4*a+j]**2 + (indxtran * peristdv[4*a+j])**2) * 24 # [hr]
    #    
    #    figr, axis = plt.subplots(figsize=figrsize)
    #    for j, strgplan in enumerate(liststrgplan):
    #        axis.plot(timetranpred[a][j], stdvtimetran[j], color=listcolrplan[j], label=listlabltoii[j])
    #        xlim = axis.get_xlim()
    #        timeaxis = np.linspace(xlim[0], xlim[1], 7)
    #        
    #        objttime = astropy.time.Time(timeaxis, format='jd', scale='utc', out_subfmt='date_hm')
    #        labltimeaxis = [labltimeaxistem[:10] for labltimeaxistem in objttime.iso]
    #        axis.set_xticks(timeaxis)
    #        axis.set_xticklabels(labltimeaxis, rotation=20)
    #    axis.set_ylabel('$\sigma$ [hour]')
    #    axis.legend()
    #    plt.subplots_adjust()
    #    path = pathimag + 'stdvtimetran_%d.%s' % (a, strgplotextn)
    #    print('Writing to %s...' % path)
    #    plt.savefig(path)
    #    plt.close()
    
    
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
    
    # plot follow-up light curves
    dictarry = {}
    
    ## MEarth
    for strgdate in  ['0317', '0303', '0311']:
        strgarrybase = 'mear' + strgdate
        strgtelemerg = 'mear%stl00' % strgdate
        dictarry[strgtelemerg] = []
        
        if strgdate == '0311':
            strgplanmear = '02'
        else:
            strgplanmear = '01'
        for strg in ['11', '12', '13', '15', '16', '17']:
            strgarrytele = strgarrybase + 'tl' + strg
            path = pathdata + 'TFOP/MEarth-South/TIC260647166-%s_2020%s_MEarth-South_defocus_lcs/tel%s/T1_ap2.txt' % (strgplanmear, strgdate, strg)
            print('Reading from %s...' % path)
            data = np.loadtxt(path, skiprows=17)
            numbtime = data.shape[0]
            dictarry[strgarrytele] = np.empty((numbtime, 3))
            dictarry[strgarrytele][:, 0] = data[:, 0]
            dictarry[strgarrytele][:, 1] = data[:, 18]
            dictarry[strgarrytele][:, 2] = data[:, 2]
            
            rflx, stdvrflx = ephesus.util.retr_rflxfromdmag(dictarry[strgarrytele][:, 1], dictarry[strgarrytele][:, 2])
            dictarry[strgarrytele][:, 1] = rflx
            dictarry[strgarrytele][:, 2] = stdvrflx
    
            dictarry[strgtelemerg].append(dictarry[strgarrytele])
        # merge data from different defocused telescopes
        dictarry[strgtelemerg] = np.concatenate(dictarry[strgtelemerg], 0)
        
        ## sort the data in time
        indx = np.argsort(dictarry[strgtelemerg][:, 0])
        dictarry[strgtelemerg] = dictarry[strgtelemerg][indx, :]
    
        # save the MEarth-South light curve
        path = pathinptalle + 'MEARTHS.csv'
        print('Writing to %s...' % path)
        np.savetxt(path, dictarry[strgtelemerg], delimiter=',')
    
    ## LCOGT
    pathmeas = pathdata + 'TFOP/measurements/'
    listextnlcog = [ \
                    ['01', '20200302', 'SAAO'], \
                    ['01', '20200317', 'CTIO'], \
                    ['02', '20200111', 'CTIO'], \
                    ['02', '20200131', 'SAAO'], \
                    ['02', '20200311', 'CTIO'], \
                    ['03', '20200202', 'SAAO'], \
                    ['03', '20200311', 'CTIO'], \
                    ['04', '20200111', 'SAAO'], \
                ]
    
    for k in range(len(listextnlcog)):
        strgplanlcog = listextnlcog[k][0]
        strgdatelcog = listextnlcog[k][1]
        strglocalcog = listextnlcog[k][2]
        strgdatelcogextn = 'lcog' + strgdatelcog
        path = pathmeas + 'TIC260647166-%s_%s_LCO-%s-1m_measurements.tbl' % (strgplanlcog, strgdatelcog, strglocalcog)
        print('Reading from %s...' % path)
        objtfile = open(path, 'r')
        listtime = []
        liststdvrflx = []
        listrflx = []
        for n, line in enumerate(objtfile):
            linesplt = line.split('\t') 
            if n == 0:
                cols = np.array(linesplt)
                indxtime = np.where(cols == 'BJD_TDB')[0][0]
                indxrflx = np.where(cols == 'rel_flux_T1')[0][0]
                indxstdvrflx = np.where(cols == 'rel_flux_err_T1')[0][0]
            else:
                listtime.append(float(linesplt[indxtime]))
                liststdvrflx.append(float(linesplt[indxstdvrflx]))
                listrflx.append(float(linesplt[indxrflx]))
        numbtime = len(listtime)
        dictarry[strgdatelcogextn] = np.empty((numbtime, 3))
        dictarry[strgdatelcogextn][:, 0] = np.array(listtime)
        dictarry[strgdatelcogextn][:, 1] = np.array(listrflx)
        dictarry[strgdatelcogextn][:, 2] = np.array(liststdvrflx)
    
    # plot LCO-CTIO detection
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
    
    timefutu = 2459000
    numbtran = np.empty((numbpred, numbplan), dtype=int)
    for a in indxpred:
        for j in indxplan:
            numbtran[a, j] = (timefutu - epoc[a][j]) / peri[a][j]
    timetranpred = [[np.empty(numbtran[a, j]) for j in indxplan] for a in indxpred]
    for a in indxpred:
        for j in indxplan:
            indxtran = np.arange(numbtran[a, j])
            timetranpred[a][j] = epoc[a][j] + indxtran * peri[a][j]
    
    numbbins = 40
    for strgarry in liststrg:
        
        if strgarry.startswith('lcoc'):
            colr = 'red'
            timeoffs = timeoffslcoc
        else:
            colr = 'k'
            timeoffs = timeoffsfolw
        
        strgarrybind = strgarry + 'bind'
        dictarry[strgarrybind] = ephesus.util.rebn_lcur(dictarry[strgarry], numbbins)
        
        figr, axis = plt.subplots(figsize=(6, 4))
        axis.errorbar(dictarry[strgarry][:, 0] - timeoffs, dictarry[strgarry][:, 1], ls='', marker='o', color='grey', ms=1)
        
        axis.errorbar(dictarry[strgarrybind][:, 0] - timeoffs, dictarry[strgarrybind][:, 1], \
                                                                yerr=dictarry[strgarrybind][:, 2], ls='', marker='o', color=colr, ms=3)
        
        # overplot the model for LCO light curve of planet d
        if strgarry == 'lcocdetr':
            pathalle = pathbase + 'allesfits/allesfit_orbt_folw/'
            print('Reading from %s...' % pathalle)
            alles = allesfitter.allesclass(pathalle)
            lcurmodl = alles.get_posterior_median_model('LCO', 'flux', xx=dictarry[strgarry][:, 0])
            axis.plot(dictarry[strgarry][:, 0] - timeoffs, lcurmodl, color='b')
        
        # overplot the TESS-predicted mid-transit time
        for a in indxpred:
            # only overplot my posterior
            if a != 2:
                continue
            for j in indxplan:
                indxtranplan = np.argmin(abs(timetranpred[a][j] - np.mean(dictarry[strgarry][:, 0])))
                if abs(timetranpred[a][j][indxtranplan] - np.mean(dictarry[strgarry][:, 0])) < 0.2:
                    axis.axvline(timetranpred[a][j][indxtranplan] - timeoffs, ls=liststylline[a], alpha=0.8, color='k')
        axis.set_xlabel('Time [BJD - %d]' % timeoffs)
        axis.set_ylabel('Relative Flux')
        plt.subplots_adjust(left=0.3, bottom=0.2)
        path = pathimaginit + 'lcurtfop_%s.%s' % (strgarry, strgplotextn)
        print('Writing to %s...' % path)
        plt.savefig(path)
        plt.close()
    
    # injection recovery

numbruns = 2
indxruns = np.arange(numbruns)

for a in indxruns:
    
    if a == 0:
        priotype = 'inpt'
        
        rratprio = np.concatenate([rrat[0][indxsorttcee], np.array([np.sqrt(0.23e-3)])])
        epocprio = np.concatenate([epoc[0][indxsorttcee], np.array([2458570.6781])])
        periprio = np.concatenate([peri[0][indxsorttcee], np.array([10.9113])])
        listtypeallemodl = ['pla5']
        dictdictallesett = {'pla5': {}}
        dictdictallesett['pla5']['use_host_density_prior'] = 'False'
        dictdictallesett['pla5']['mcmc_total_steps'] = '6000'
        dictdictallesett['pla5']['mcmc_burn_steps'] = '5000'
        dictdictallepara = {'pla5': {}}
        stdvperiprio = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        stdvepocprio = np.array([0.04, 0.04, 0.04, 0.04, 0.04])
        indxplanalle = np.array([4])
        liststrgplan = ['b', 'c', 'd', 'e', 'f']
    else:
        periprio = None
        epocprio = None
        rratprio = None
        priotype = 'exof'
        listtypeallemodl = ['orbt']
        indxplanalle = None
        liststrgplan = None
        stdvperiprio = None
        stdvepocprio = None
        dictdictallesett = None
        dictdictallepara = None
    
    if a == 0:
        continue

    for a in range(2):
    
        if a == 1:
            radiplan += [2.0]
            rsma += [0.88 / (215. * 0.1758)]
            epoc += [2458793.2786]
            peri += [29.54115]
            cosi += [0.]
                
        pexo.main.init( \
                       toiitarg=1233, \
                       strgtarg=strgtarg, \
                       labltarg=labltarg, \
                       
                       priotype=priotype, \
                       
                       #boolmaskqual=False, \
                        
                       # plot exoplanet properties
                       boolplotprop=True, \
                        
                       #booldatatser=False, \
                       boolprocprio=False, \
                       boolexecalle=True, \
                       
                       boolexof=False, \

                       listtypeallemodl=listtypeallemodl, \
                       dictdictallesett=dictdictallesett, \
                       dictdictallepara=dictdictallepara, \
                        
                       stdvperiprio=stdvperiprio, \
                       periprio=periprio, \
                       epocprio=epocprio, \
                       rratprio=rratprio, \
                       stdvepocprio=stdvepocprio, \

                       liststrgplan=liststrgplan, \
                       
                       indxplanalle=indxplanalle, \

                       # RM proposal
                       #listlablinst=[['TESS'], ['PFS']], \
                       #listdatatype=[['real'], ['mock']], \
                       #vsiistarprio=4.7, \
                       #stdvvsiistarprio=0.5, \
                       #lambstarprio=45., \
                       #stdvlambstarprio=10., \
    
                       radistar=gdat.radistar, \
                       stdvradistar=gdat.stdvradistar, \
                       massstar=gdat.massstar, \
                       stdvmassstar=gdat.stdvmassstar, \
                       tmptstar=gdat.tmptstar, \
                       stdvtmptstar=gdat.stdvtmptstar, \
                       
                       #rvsaprio=rvsaprio, \
    
                      )
    
