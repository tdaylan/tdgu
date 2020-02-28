import pexo.main
import sys
import numpy as np

def cnfg_WASP0121():

    strgtarg = 'WASP-121'
    strgmast = 'WASP-121'
    labltarg = 'WASP-121'
    ticitarg = 22529346
    strgtoii = ['495.01']
    pexo.main.main( \
         strgmast=strgmast, \
         ticitarg=ticitarg, \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         strgtoii=strgtoii, \
         boolphascurv=True, \
        )


def cnfg_GJ299():
    
    radistar = 0.175016 * 9.95 # [R_J]
    massstar = 0.14502 * 1048. # [M_J]
    pexo.main.main( \
         ticitarg=334415465, \
         strgtarg='gj299', \
         labltarg='GJ 299', \
         weigsplndetr=1e-5, \
         epocpmot=2019.3, \
         radistar=radistar, \
         massstar=massstar, \
        )


def cnfg_josh():
    
    strgtarg = 'HD118203'
    strgmast = 'HD 118203'
    labltarg = 'HD 118203'
    
    listlimttimemask = np.array([ \
                                [0, 1712], \
                                [1724.5, 1725.5], \
                                ])
    listlimttimemask += 2457000
    epocprio = np.array([2458712.662354])
    periprio = np.array([6.134842])
    duraprio = np.array([5.6457]) / 24. # [day]
    rratprio = np.sqrt(np.array([3516.19165]) * 1e-6)
    pexo.main.main( \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         strgmast=strgmast, \
         epocprio=epocprio, \
         periprio=periprio, \
         duraprio=duraprio, \
         rratprio=rratprio, \
         listlimttimemask=listlimttimemask, \
        )


def cnfg_toii1339():
    
    strgtarg = 'TOI1339'
    strgmast = '269701147'
    labltarg = 'TOI 1339'
    
    epocprio = np.array([2458715.354492, 2458726.054199, 2458743.5534])
    periprio = np.array([8.880832, 28.579357, 38.3499])
    duraprio = np.array([3.0864, 4.4457, 5.5336]) / 24. # [day]
    rratprio = np.array([0.0334, 0.0314, 0.0310])
    pexo.main.main( \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         strgmast=strgmast, \
         epocprio=epocprio, \
         periprio=periprio, \
         duraprio=duraprio, \
         rratprio=rratprio, \
        )


def cnfg_HATP19():
    
    ticitarg = 267650535
    strgmast = 'HAT-P-19'
    labltarg = 'HAT-P-19'
    strgtarg = 'hatp0019'
    #strgtoii = '495' 
    
    pexo.main.main( \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         strgmast=strgmast, \
         ticitarg=ticitarg, \
        )


def cnfg_toii0203():

    strgtarg = 'TOI203'
    strgmast = '259962054'
    pexo.main.main( \
         strgtarg=strgtarg, \
         strgmast=strgmast, \
        )


def cnfg_WD1856():

    epocprio = np.array([2458708.978112])
    periprio = np.array([1.4079342])
    duraprio = np.array([1. / 24. / 2.])
    strgtarg = 'WD1856'
    ticitarg = 267574918
    strgmast = 'TIC 267574918'
    labltarg = 'WD-1856'
    
    #strgtoii = '' 
    
    print('HACKING! MAKING UP THE STAR RADIUS')
    pexo.main.main( \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         strgmast=strgmast, \
         ticitarg=ticitarg, \
         boolmakeanim=True, \
         #maxmnumbstartcat=40, \
         #makeprioplot=False, \

         infetype='trap', \
         #dilucorr=0.01, \
         jmag=15.677, \
         #contrati=10, \
         # temp
         tmptstar=5000., \
         datatype='sapp', \
         radistar=10., \
         massstar=1000., \
         booltlss=False, \
         epocprio=epocprio, \
         periprio=periprio, \
         duraprio=duraprio, \
        )


def cnfg_cont():

    strgtarg = 'cont'
    ticitarg = 1717706276
    strgmast = 'TIC 1717706276'
    labltarg = 'Cont'
    
    pexo.main.main( \
         strgtarg=strgtarg, \
         labltarg=labltarg, \
         strgmast=strgmast, \
         ticitarg=ticitarg, \
         boolmakeanim=True, \
         #maxmnumbstartcat=40, \
         contrati=10, \
         datatype='tcat', \
         booltlss=False, \
        )


globals().get(sys.argv[1])(*sys.argv[2:])


