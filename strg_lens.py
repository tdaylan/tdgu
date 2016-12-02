from __init__ import *

def pcat_lens_mock():
    
    minmflux = deg2rad(5e-2 / 3600.)
    maxmflux = deg2rad(5e-1 / 3600.)

    numbiter = 1
    for k in range(numbiter):
        gridchan = pcat.main.init( \
                                  #verbtype=2, \
                                  #makeplot=False, \
                                  numbswep=100000, \
                                  numbswepplot=10000, \
                                  factthin=100, \
                                  diagmode=False, \
                                  #prophypr=False, \
                                  #proppsfp=False, \
                                  #propbacp=False, \
                                  #proplenp=False, \
                                  #optiprop=True, \
                                  #pntstype='lens', \
                                  
                                  stdvprophypr=0.1, \
                                  stdvpropbacp=5e-3, \
                                  stdvproplenp=5e-4, \
                                  #probtran=0., \
                                  #optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  stdvflux=1e-5, \
                                  stdvlbhl=0.01, \
                                  maxmnumbpnts=array([0]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mockbacp=zeros((1, 1)), \
                                  #mocknumbpnts=array([20]), \
                                 )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


