from __init__ import *

def pcat_lens_mock():
    
    minmflux = deg2rad(1e-3 / 3600.)
    maxmflux = deg2rad(1e-1 / 3600.)
    gridchan = pcat.main.init( \
                              verbtype=2, \
                              #makeplot=False, \
                              numbswep=100, \
                              numbswepplot=500, \
                              factthin=10, \
                              pntstype='lens', \
                              exprinfo=False, \
                              indxenerincl=arange(1), \
                              indxevttincl=arange(1), \
                              strgexpo=1e28, \
                              datatype='mock', \
                              exprtype='hubb', \
                              stdvflux=0.01, \
                              stdvlbhl=0.01, \
                              maxmnumbpnts=array([4]), \
                              #minmflux=1e-4, \
                              #maxmflux=1e-1, \
                              minmflux=minmflux, \
                              maxmflux=maxmflux, \
                              mockbacp=zeros((1, 1)), \
                              #mocknumbpnts=array([100]), \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


