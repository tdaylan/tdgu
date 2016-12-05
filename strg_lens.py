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
                                  #propcomp=False, \
                                  stdvprophypr=0.1, \
                                  stdvpropbacp=5e-3, \
                                  stdvproplenp=5e-4, \
                                  #probtran=0., \
                                  optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  stdvflux=1e-5, \
                                  stdvlbhl=0.01, \
                                  strgback=[0.01], \
                                  maxmnumbpnts=array([10]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mocknumbpnts=array([5]), \
                                 )

if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass


