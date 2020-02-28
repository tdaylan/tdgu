import tcat.main
import sys
import numpy as np

def cnfg_WASP0121():
    
    strgmast = 'WASP-121b'
    tcat.main.main( \
               strgmast=strgmast, \
              )
        

def cnfg_Luhman16():
    
    tcat.main.main( \
         rasctarg=162.328812, \
         decltarg=-53.319467, \
         labltarg='Luhman 16', \
         weigsplndetr=1e-6, \
         evaltype='ontf', \
         strgtarg='Luhman16', \
        )
        

def cnfg_TOI1233():
    
    tcat.main.main( \
         #boolfittoffs=True, \
         labltarg='TOI 1233', \
         strgtarg='TOI1233', \
         ticitarg=260647166, \
         cntpscaltype='asnh', \
         strgdata='ffim', \
        )
        

def cnfg_WD1856():

    strgtarg = 'WD1856'
    ticitarg = 267574918
    strgmast = 'TIC 267574918'
    labltarg = 'WD 1856'
    
    tcat.main.main( \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         #strgmast=strgmast, \
         ticitarg=ticitarg, \
         boolmakeanim=True, \
         #datatype='sapp', \
        )


def cnfg_GJ299():
    
    tcat.main.main( \
         boolfittoffs=True, \
         labltarg='GJ 299', \
         maxmnumbstar=4, \
         strgtarg='gj299', \
         ticitarg=334415465, \
         epocpmot=2019.3, \
        )
        

def cnfg_spec():
    
    path = os.environ['TCAT_DATA_PATH'] + '/data/List_for_MIT_pilot.txt'
    data = np.loadtxt(path, delimiter='\t', skiprows=1)
    numbtarg = data.shape[0]
    indxtarg = np.arange(numbtarg)
    for k in indxtarg:
        ticitarg = int(data[k, 2])
        tcat.main.main( \
             ticitarg=ticitarg, \
             labltarg='TIC %s' % ticitarg, \
             strgtarg='speculus_%s' % ticitarg, \
            )
        
    
def cnfg_saul():
    
    path = os.environ['TCAT_DATA_PATH'] + '/data/list_saul.txt'
    strgbase = 'saul'

    init_list( \
              path, \
              strgbase, \
             )


def cnfg_maxg122319():
    
    pathdata = os.environ['TCAT_DATA_PATH'] + '/data/'
    path = pathdata + 'Amaury_EBs.csv'
    #path = pathdata + 'SPECULOOS_TIC_WTV.csv'
    listticitarg = []
    data = np.genfromtxt(path, delimiter=',', skip_header=1)
    
    listtici = data[:, 1]
    #listgaia = data[:, 2]
    
    listrasctarg = data[:, 2]
    listdecltarg = data[:, 3]
    
    strgbase = 'maxg122319'
   
    numbtarg = listrasctarg.size
    indxtarg = np.arange(numbtarg)
    listintgresu = np.empty(numbtarg)
    for k in indxtarg:
        
        tici = listtici[k]
        #labltarg = 'GID %s' % gaia
        #strgtarg = 'spec1313_%d' % int(gaia)
        
        labltarg = 'TIC %s' % tici
        strgtarg = '%d' % int(tici)
        
        listintgresu[k] = tcat.main.main( \
                               rasctarg=listrasctarg[k], \
                               decltarg=listdecltarg[k], \
                               labltarg=labltarg, \
                               strgbase=strgbase, \
                               detrtype='mfil', \
                               maxmnumbstar=4, \
                               strgtarg=strgtarg, \
                              )


def chec_runs():
    
    strg = 'spec1313'

    path = os.environ['TCAT_DATA_PATH'] + '/'
    liststrgfile = fnmatch.filter(os.listdir(path), '%s_*' % strg)
    numb = len(liststrgfile)
    listbool = np.zeros(numb, dtype=bool)
    for k, strgfile in enumerate(liststrgfile):
        liststrgextn = fnmatch.filter(os.listdir(path + strgfile + '/data/'), 'rflx_*')
        if len(liststrgextn) == 1:
            listbool[k] = True
    print('numb')
    print(numb)
    print('np.where(listbool).size')
    print(np.where(listbool)[0].size)
    print(float(np.where(listbool)[0].size) / numb)


def cnfg_test347543557():

    ticitarg = 347543557
    labltarg = 'TIC 347543557'
    strgtarg = 'test347543557'
    listlimttimeplot = [[2458428, 2458430]]
    tcat.main.main( \
         ticitarg=ticitarg, \
         labltarg=labltarg, \
         strgtarg=strgtarg, \
         boolplotframtotl=True, \
         listlimttimeplot=listlimttimeplot, \
        )


def cnfg_GRB191016A():
    
    rasctarg = 30.2695
    decltarg = 24.5099
    labltarg = 'GRB191016A'
    strgtarg = 'GRB191016A'
    listpathtescfile = ['/Users/tdaylan/Downloads/tess-s0017-1-4_30.269500_24.509900_11x11_astrocut.fits']
    listlimttimeplot = []
    for timedelt in [1.]:
        listlimttimeplot.append(2458772.67 + np.array([-timedelt, timedelt]))
    listtimeplotline = [2458772.67291666667]
    boolfittoffs = True
    boolcuttqual = False
    tcat.main.main( \
         rasctarg=rasctarg, \
         decltarg=decltarg, \
         labltarg=labltarg, \
         strgtarg=strgtarg, \
         boolcuttqual=boolcuttqual, \
         boolfittoffs=boolfittoffs, \
         numbside=9, \
         booldetr=False, \
         listtimeplotline=listtimeplotline, \
         #listpathtescfile=listpathtescfile, \
         listlimttimeplot=listlimttimeplot, \
        )


globals().get(sys.argv[1])()


