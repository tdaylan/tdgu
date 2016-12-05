from __init__ import *


def pcat_lens_mock_arry():

    #liststrgback = linspace(0.1, 3., 3.)
    liststrgback = logspace(-1., 1., 3.)

    minmflux = deg2rad(5e-2 / 3600.)
    maxmflux = deg2rad(5e-1 / 3600.)

    for k in range(len(liststrgback)):
       
        strgback = [liststrgback[k]]
        gridchan = pcat.main.init( \
                                  numbswep=100000, \
                                  numbswepplot=10000, \
                                  factthin=100, \
                                  diagmode=False, \
                                  optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  strgback=strgback, \
                                  maxmnumbpnts=array([10]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mocknumbpnts=array([5]), \
                                 )

    
def pcat_lens_mock():
    
    listminmflux = array([deg2rad(1e-3 / 3600.)])
#    listminmflux = deg2rad(logspace(-3., -2., 2) / 3600.)
    liststrgback = logspace(-2., 0., 3.)
    maxmflux = deg2rad(5e-1 / 3600.)

    numbiter = 1
    for k in range(numbiter):
        seedstat = get_state()

        for k in range(len(listminmflux)):
            for n in range(len(liststrgback)):
                minmflux = listminmflux[k]
                strgback = [liststrgback[n]]
                gridchan = pcat.main.init( \
                                  #verbtype=2, \ 
                                  seedstat=seedstat, \
                                  #makeplot=False, \
                                  numbswep=10000, \
                                  numbswepplot=3000, \
                                  factthin=10, \
                                  diagmode=False, \
                                  #prophypr=False, \
                                  proppsfp=False, \
                                  #propbacp=False, \
                                  #proplenp=False, \
                                  #propcomp=False, \
                                  #stdvprophypr=0.1, \
                                  #stdvpropbacp=5e-3, \
                                  #stdvproplenp=5e-4, \
                                  #probtran=0., \
                                  #optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  strgback=strgback, \
                                  maxmnumbpnts=array([20]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mocknumbpnts=array([10]), \
                                 )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


