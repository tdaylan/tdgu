from __init__ import *

def pcat_lens_mock():
    
    minmflux = deg2rad(1. / 3600.)
    maxmflux = deg2rad(1e2 / 3600.)
    gridchan = pcat.main.init( \
                              verbtype=2, \
                              numbswep=10000, \
                              #numbproc=20, \
                              numbburn=0, \
                              factthin=10, \
                              numbswepplot=30000, \
                              pntstype='lens', \
                              exprinfo=False, \
                              indxenerincl=arange(1), \
                              indxevttincl=arange(1), \
                              strgback=[1.], \
                              strgexpo=1e28, \
                              datatype='mock', \
                              exprtype='hubb', \
                              scalmaps='asnh', \
                              stdvflux=0.01, \
                              stdvlbhl=0.01, \
                              #optiprop=True, \
                              boolpropnormback=False, \
                              maxmnumbpnts=array([100]), \
                              #minmflux=1e-4, \
                              #maxmflux=1e-1, \
                              minmflux=minmflux, \
                              maxmflux=maxmflux, \
                              mocknormback=zeros((1, 1)), \
                              mocknumbpnts=array([100]), \
                             )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


