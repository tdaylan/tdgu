from __init__ import *

def pcat_lens_mock():
    
    minmflux = deg2rad(1. / 3600.)
    maxmflux = deg2rad(1e2 / 3600.)
    gridchan = pcat.main.init( \
                              numbswep=200000, \
                              numbswepplot=20000, \
                              pntstype='lens', \
                              exprinfo=False, \
                              indxenerincl=arange(1), \
                              indxevttincl=arange(1), \
                              strgback=['zero'], \
                              strgexpo=1e25, \
                              datatype='mock', \
                              exprtype='hubb', \
                              boolpropnormback=False, \
                              maxmnumbpnts=array([100]), \
                              #minmflux=1e-4, \
                              #maxmflux=1e-1, \
                              minmflux=minmflux, \
                              maxmflux=maxmflux, \
                              mocknumbpnts=array([100]), \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


