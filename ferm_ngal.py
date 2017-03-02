from __init__ import *

def make_maps():

    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec7']
    gdat.enertype = ['pnts']
    gdat.timefrac = [0.1]
    tdpy.util.make_maps_main(gdat, os.environ["FERM_NGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec7', 'pnts', 'ngal', os.environ["FERM_NGAL_DATA_PATH"], 256, 'tim0')


def prep_maps():
    
    tdpy.util.prep_maps('rec7', 'pnts', 'ngal', os.environ["FERM_NGAL_DATA_PATH"], 256, 'tim0')


def plot_spec():

    # Fermi-LAT best-fit components at the NGP
    # temp
    if False:
        if gdat.datatype == 'mock':
            pass
        else:
            if gdat.exprtype == 'ferm':
                listname = ['data', 'pion', 'invc', 'brem', 'pnts', 'isot']
                listmrkr = ['o', 's', 'p', '*', 'D', '^']
                listcolr = ['g', 'g', 'g', 'g', 'g', 'g']
                listlabl = ['Fermi-LAT Data', r'Fermi-LAT $\pi^0$', 'Fermi-LAT ICS', 'Fermi-LAT Brem', 'Fermi-LAT PS', 'Fermi-LAT Iso']
                for k, name in enumerate(listname):
                    path = os.environ["PCAT_DATA_PATH"] + '/fermspec' + name + '.csv'
                    data = loadtxt(path)
                    enertemp = data[:, 0] # [GeV]
                    fluxtemp = data[:, 1] * 1e-3 # [GeV/cm^2/s/sr]
                    fluxtemp = interp(gdat.meanener, enertemp, fluxtemp)
                    #fluxtemp = interpolate.interp1d(enertemp, fluxtemp)(gdat.meanener)
                    axis.plot(gdat.meanener, fluxtemp, marker=listmrkr[k], color=listcolr[k], label=listlabl[k])

    
def pcat_ferm_inpt_ngal_intr( \
                        strgexprflux='fermflux_cmp0_ngal.fits', \
                        strgexpo='fermexpo_cmp0_ngal.fits', \
                       ): 

    karg = {}
    karg['numbswep'] = 3000
    karg['numbburn'] = 2500
    karg['factthin'] = 500
    karg['numbproc'] = 6
    karg['numbswepplot'] = 100000
    #karg['propwithsing'] = True
    karg['inittype'] = 'refr'
    karg['indxenerincl'] = arange(1, 4)
    karg['indxevttincl'] = arange(2, 4)
    karg['proppsfn'] = False
    karg['lgalcntr'] = 0.
    karg['bgalcntr'] = pi / 2.
    karg['back'] = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
    karg['strgexpo'] = strgexpo
    karg['maxmgangdata'] = 20. / 180. * pi
    karg['minmflux'] = 3e-11
    karg['maxmflux'] = 1e-7
    karg['maxmnumbpnts'] = array([200])
    karg['sinddisttype'] = ['atan']
    karg['strgexprflux'] = strgexprflux
    
    return karg


def pcat_ferm_inpt_ngal():
    karg = pcat_ferm_inpt_ngal_intr()
    pcat.main.init(**karg)


def pcat_ferm_inpt_ngal_cmp1():
    karg = pcat_ferm_inpt_ngal_intr(strgexprflux='fermflux_cmp1_ngal.fits', strgexpo='fermexpo_cmp1_ngal.fits')
    pcat.main.init(**karg)


def pcat_ferm_inpt_ngal_cmp2():
    karg = pcat_ferm_inpt_ngal_intr(strgexprflux='fermflux_cmp2_ngal.fits', strgexpo='fermexpo_cmp2_ngal.fits')
    pcat.main.init(**karg)


def pcat_ferm_inpt_ngal_cmp3():
    karg = pcat_ferm_inpt_ngal_intr(strgexprflux='fermflux_cmp3_ngal.fits', strgexpo='fermexpo_cmp3_ngal.fits')
    pcat.main.init(**karg)


def pcat_ferm_inpt_ngal_tim4():
    karg = pcat_ferm_inpt_ngal_intr(strgexprflux='fermflux_rec7_pnts_ngal_0256_tim4.fits', strgexpo='fermexpo_rec7_pnts_ngal_0256_tim4.fits')
    pcat.main.init(**karg)


def pcat_ferm_mock_ngal():
     
    pcat.main.init( \
                   numbswep=2000, \
                   factthin=500, \
                   numbburn=0, \
                   optiprop=True, \
                   #optipropllik=True, \
                   #makeplot=False, \
                   #verbtype=2, \
                   #optiprop=True, \
                   indxenerincl=arange(1, 4), \
                   indxevttincl=arange(2, 4), \
                   maxmgangdata=4./180.*pi, \
                   lgalcntr=0., \
                   bgalcntr=pi / 2., \
                   back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                   strgexpo='fermexpo_cmp0_ngal.fits', \
                   numbsideheal=256, \
                   #trueminmflux=7e-11, \
                   maxmnumbpnts=array([100]), \
                   truenumbpnts=array([50]), \
                  )

globals().get(sys.argv[1])()
