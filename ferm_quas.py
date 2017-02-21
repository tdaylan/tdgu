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
            truespatdisttype = None
            truespatdisttype = None
   
        pcat.main.init( \
                       seedstat=seedstat, \
                       numbswep=10000, \
                       verbtype=2, \
                       spatdisttype=spatdisttype, \
                       lgalprio=lgalprio, \
                       bgalprio=bgalprio, \
                       lgalcntr=0., \
                       bgalcntr=pi / 2., \
                       strgexpo='fermexpo_cmp0_ngal.fits', \
                      )

globals().get(sys.argv[1])()
