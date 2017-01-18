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
    karg['numbswep'] = 10000
    karg['numbswepplot'] = 3000
    karg['diagmode'] = True
    karg['verbtype'] = 2
    karg['factthin'] = 90
    karg['randinit'] = False
    karg['indxenerincl'] = arange(2, 4)
    karg['indxevttincl'] = arange(3, 4)
    karg['lgalcntr'] = 0.
    karg['bgalcntr'] = pi / 2.
    karg['back'] = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
    karg['strgexpo'] = strgexpo
    karg['minmflux'] = 3e-11
    karg['maxmflux'] = 1e-7
    karg['strgexprflux'] = strgexprflux
    karg['maxmnumbpnts'] = array([20])
    
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
                   numbswep=1000, \
                   numbswepplot=4000, \
                   factthin=90, \
                   proppsfp=False, \
                   indxenerincl=arange(1, 4), \
                   indxevttincl=arange(2, 4), \
                   verbtype=2, \
                   #makeplot=False, \
                   lgalcntr=0., \
                   #propcova=True, \
                   #diagmode=True, \
                   bgalcntr=pi / 2., \
                   back=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                   strgexpo='fermexpo_cmp0_ngal.fits', \
                   minmflux=3e-11, \
                   maxmflux=1e-7, \
                   numbsideheal=256, \
                   maxmnumbpnts=array([10]), \
                   truenumbpnts=array([5]), \
                   truefluxdistslop=array([2.0]), \
                  )

globals().get(sys.argv[1])()
