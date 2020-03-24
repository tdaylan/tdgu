import pexo.main
import numpy as np
import os
import matplotlib.pyplot as plt
import astropy
from tdpy.util import summgene

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

factmsmj = 1048.
factrsrj = 9.95
massstar = 0.97 * factmsmj # [M_J]
radistar = 0.888 * factrsrj # [R_J]

pathbase = os.environ['PEXO_DATA_PATH'] + '/TOI1233/'
pathdata = pathbase + 'data/'
pathimag = pathbase + 'imag/'
pathimagstab = pathimag + 'stab/'
os.system('mkdir -p %s' % pathimagstab)

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

indxplan = np.arange(numbplan, dtype=int)

numbtran = np.empty(numbplan, dtype=int)
timefutu = 2459000
for j in indxplan:
    numbtran[j] = (timefutu - epoc[j]) / peri[j]
time = [np.empty(numbtran[j]) for j in indxplan]
labltime = [np.empty(numbtran[j], dtype=object) for j in indxplan]
stdvtimetran = [np.empty(numbtran[j]) for j in indxplan]
for a in range(2):
    for j in indxplan:
        indxtran = np.arange(numbtran[j])
    
        timeexof = epoc[4*a+j] + indxtran * peri[4*a+j]
        
        objttime = astropy.time.Time(timeexof, format='jd', scale='utc', out_subfmt='date_hm')
        
        time[j] = objttime.jd
        labltime[j] = objttime.iso
        stdvtimetran[j] = np.sqrt(epocstdv[4*a+j]**2 + (indxtran * peristdv[4*a+j])**2) * 24 # [hr]
    
    figr, axis = plt.subplots(figsize=figrsize)
    for j, strgplan in enumerate(liststrgplan):
        axis.plot(time[j], stdvtimetran[j], color=listcolrplan[j], label=listlabltoii[j])
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


# plot Stephen Kane's outpu
#[]Time (years)         a          e          i          mass          long       node         M    
listlabl = ['Time [year]', '$a$ [AU]', 'Eccentricity', 'Inclination [degree]', 'Mass [$M_E$]', 'Longitude of ascending node [degree]', 'node', 'Mean anomaly']
liststrg = ['time', 'smax', 'ecce', 'incl', 'mass', 'long', 'node', 'M']
numbcols = len(liststrg)
indxcols = np.arange(numbcols)

numbplan = len(strgtoii)
for j, strgplan in enumerate(liststrgplan):
    path = pathdata + 'stabilityfiles/planet%d.aei' % (j + 1)
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
    #    path = pathimagstab + '%s_%d.%s' % (liststrg[k], j, strgplotextn)
    #    print('Writing to %s...' % path)
    #    plt.savefig(path)
    #    plt.close()
    
    listecce[:, j] = data[:, 2]

figr, axis = plt.subplots(1, numbplan, figsize=figrsizeydob, sharey=True)
bins = np.linspace(0., np.amax(listecce), 100)
for j, strgplan in enumerate(liststrgplan):
    axis[j].hist(listecce[:, j], color=listcolrplan[j])
    axis[j].set_xlabel('$\epsilon_{%s}$' % liststrgplan[j])
axis[0].set_ylabel('Number of samples')
plt.subplots_adjust(bottom=0.2)
path = pathimagstab + 'histecce.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()


# plot Keivan Stassun's model
path = pathdata + 'specstardata.dat'
arry = np.loadtxt(path, skiprows=2)
wlendata = arry[:, 0]
stdvwlendata = arry[:, 1]
fluxdata = arry[:, 2]
stdvfluxdata = arry[:, 3]
path = pathdata + 'specstarmodl.dat'
arry = np.loadtxt(path, skiprows=2)
wlenmodl = arry[:, 0]
fluxmodl = 10**arry[:, 1]
figr, axis = plt.subplots(figsize=figrsize)
axis.errorbar(wlendata, fluxdata, xerr=stdvwlendata, yerr=stdvfluxdata, ls='', marker='o', color='k', ms=1)
axis.plot(wlenmodl, fluxmodl, color='b')
axis.set_xlabel('Wavelenth [$\mu$m]')
axis.set_ylabel('Flux [erg s$^{-1}$ cm$^{-2}$)]')
axis.set_yscale('log')
axis.set_ylim([1e-13, None])
plt.subplots_adjust(left=0.2, bottom=0.2)
path = pathimagstab + 'specstar.%s' % (strgplotextn)
print('Writing to %s...' % path)
plt.savefig(path)
plt.close()

pexo.main.main( \
                   strgtarg=strgtarg, \
                   labltarg=labltarg, \
                   
                   strgmast=strgmast, \
                   #ticitarg=ticitarg, \
                
                   booltlss=False, \
                   smaxprio=smaxprio, \
                   
                   toiitarg=1233, \
                   priotype='exof', \

                   liststrgplan=liststrgplan, \

                   inclprio=inclprio, \

                   massstar=massstar, \
                   radistar=radistar, \
                    
                  )

