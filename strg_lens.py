from __init__ import *

def pcat_lens_mock():
    
    minmflux = deg2rad(0.5e-20 / 3600.)
    maxmflux = deg2rad(1e-19 / 3600.)

    numbiter = 1
    for k in range(numbiter):
        gridchan = pcat.main.init( \
                                  verbtype=2, \
                                  #makeplot=False, \
                                  numbswep=10, \
                                  diagmode=False, \
                                  numbswepplot=10000, \
                                  factthin=1, \
                                  #prophypr=False, \
                                  proppsfp=False, \
                                  propbacp=False, \
                                  proplenp=False, \
                                  #optiprop=True, \
                                  #pntstype='lens', \
                                  stdvprophypr=0.1, \
                                  stdvpropbacp=5e-3, \
                                  stdvproplenp=5e-4, \
                                  probtran=0., \
                                  #optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  stdvflux=1e-4, \
                                  stdvlbhl=0.01, \
                                  maxmnumbpnts=array([3]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mockbacp=zeros((1, 1)), \
                                  mocknumbpnts=array([2]), \
                                 )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


