from __init__ import *

def pcat_lens_mock():
    
    minmflux = deg2rad(1e-2 / 3600.)
    maxmflux = deg2rad(1e-1 / 3600.)

    numbiter = 1
    for k in range(numbiter):
        gridchan = pcat.main.init( \
                                  #verbtype=2, \
                                  #makeplot=False, \
                                  numbswep=100000, \
                                  diagmode=False, \
                                  numbswepplot=10000, \
                                  factthin=100, \
                                  #proppsfp=False, \
                                  #propbacp=False, \
                                  #proplenp=False, \
                                  optiprop=True, \
                                  pntstype='lens', \
                                  exprinfo=False, \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  datatype='mock', \
                                  exprtype='hubb', \
                                  stdvflux=0.01, \
                                  stdvlbhl=0.01, \
                                  maxmnumbpnts=array([40]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mockbacp=zeros((1, 1)), \
                                  mocknumbpnts=array([20]), \
                                 )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


