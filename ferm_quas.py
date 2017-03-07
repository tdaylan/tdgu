from __init__ import *

def pcat_ferm_quas_mock():
    
    seedstat = get_state()
    
    numbiter = 3
    for k in range(numbiter):
        if k == 0:
            spatdisttype = ['gaus']
            lgalprio = None
            bgalprio = None
        elif k == 1:
            spatdisttype = ['gaus']
            lgalprio = (rand(20) - 0.5) * 20. * pi / 180.
            bgalprio = (rand(20) - 0.5) * 20. * pi / 180.
        else:
            spatdisttype = None
   
        pcat.main.init( \
                       seedstat=seedstat, \
                       checprio=False, \
                       numbswep=50000, \
                       truenumbpnts=array([50]), \
                       maxmnumbpnts=array([100]), \
                       spatdisttype=spatdisttype, \
                       lgalprio=lgalprio, \
                       bgalprio=bgalprio, \
                       lgalcntr=0., \
                       bgalcntr=pi / 2., \
                       strgexpo='fermexpo_cmp0_ngal.fits', \
                      )

                  # numbswep=2000, \
                  # indxenerincl=arange(1, 4), \
                  # indxevttincl=arange(2, 4), \
                  # maxmgangdata=4./180.*pi, \
                  # lgalcntr=0., \
                  # bgalcntr=pi / 2., \
                  # back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                  # strgexpo='fermexpo_cmp0_ngal.fits', \
                  # numbsideheal=256, \
                  # #trueminmflux=7e-11, \
                  # maxmnumbpnts=array([100]), \
                  # truenumbpnts=array([50]), \
                  #)









globals().get(sys.argv[1])()
