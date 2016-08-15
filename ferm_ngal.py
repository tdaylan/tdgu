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

    
def pcat_info():
    
    minmflux = array([3e-10, 1e-10, 3e-11, 1e-11])
    numbruns = minmflux.size
    maxmnumbpnts = zeros(numbruns, dtype=int) + 1000
    numbswep = zeros(numbruns, dtype=int) + 2000000
    numbburn = numbswep / 2
    
    numbiter = minmflux.size

    listlevi = zeros(numbiter)
    listinfo = zeros(numbiter)
    
    strgexpo='fermexpo_cmp0_ngal.fits'
    strgexpr='fermflux_cmp0_ngal.fits'

    indxenerincl = arange(2, 3)
    indxevttincl = arange(3, 4)
    numbener = indxenerincl.size

    # temp
    if False:
        for k in range(numbiter):
            gridchan = pcat.main.init( \
                                  psfntype='doubking', \
                                  numbswep=numbswep[k], \
                                  numbburn=numbburn[k], \
                                  probprop=array([0.01, 0.01, 0., 0., 1., 1., 0, 0, 1., 1., 1., 1.], dtype=float), \
                                  randinit=False, \
                                  makeplot=True, \
                                  maxmgang=10., \
                                  maxmnumbpnts=array([maxmnumbpnts[k]]), \
                                  indxenerincl=indxenerincl, \
                                  indxevttincl=indxevttincl, \
                                  minmflux=minmflux[k], \
                                  maxmflux=1e-7, \
                                  regitype='ngal', \
                                  pathdata=os.environ["FERM_NGAL_DATA_PATH"], \
                                  strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                                  strgexpo=strgexpo, \
                                  datatype='inpt', \
                                  strgexpr=strgexpr, \
                                 )
        
        listlevi[k] = gridchan[-2]
        listinfo[k] = gridchan[-1]

    else:
        listlevi = array([-39222.9070569, -39226.1778779, -39300.7982166, -39521.8723332])[::-1]
        listinfo = array([91.0328924911, 98.1275628394, 98.732104824, 88.3453610331])[::-1]

    pcat.visu.plot_minmfluxinfo(minmflux, listinfo, listlevi)


def intr_ferm_expr_ngal( \
                        strgexpr='fermflux_cmp0_ngal.fits', \
                        strgexpo='fermexpo_cmp0_ngal.fits', \
                       ): 

    karg = {}
    karg['psfntype'] = 'doubking'
    karg['numbswep'] = 2000000
    karg['randinit'] = False
    karg['maxmgang'] = 20.
    karg['boolproppsfn'] = False
    karg['initfluxdistslop'] = array([1.9])
    karg['initfluxdistnorm'] = array([300])
    karg['maxmnumbpnts'] = array([500])
    karg['indxenerincl'] = arange(1, 4)
    karg['indxevttincl'] = arange(2, 4)
    karg['minmflux'] = 3e-11
    karg['maxmflux'] = 1e-7
    karg['regitype'] = 'ngal'
    karg['pathdata'] = os.environ["FERM_NGAL_DATA_PATH"]
    karg['strgback'] = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
    karg['strgexpo'] = strgexpo
    karg['datatype'] = 'inpt'
    karg['strgexpr'] = strgexpr

    return karg


def pcat_ferm_expr_ngal():
    karg = intr_ferm_expr_ngal()
    pcat.main.init(**karg)


def pcat_ferm_expr_ngal_cmp1():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_cmp1_ngal.fits', strgexpo='fermexpo_cmp1_ngal.fits')
    pcat.main.init(**karg)


def pcat_ferm_expr_ngal_cmp2():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_cmp2_ngal.fits', strgexpo='fermexpo_cmp2_ngal.fits')
    pcat.main.init(**karg)


def pcat_ferm_expr_ngal_cmp3():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_cmp3_ngal.fits', strgexpo='fermexpo_cmp3_ngal.fits')
    pcat.main.init(**karg)


def pcat_ferm_expr_ngal_tim4():
    karg = intr_ferm_expr_ngal(strgexpr='fermflux_rec7_pnts_ngal_0256_tim4.fits', strgexpo='fermexpo_rec7_pnts_ngal_0256_tim4.fits')
    pcat.main.init(**karg)


def pcat_ferm_mock_ngal():
     
    indxenerincl = arange(1, 4)
    indxevttincl = arange(2, 4)
    numbener = indxenerincl.size

    minmflux = 3e-11
    maxmflux = 2e-8
    
    pcat.main.init(psfntype='doubking', \
                   numbswep=2000000, \
                   randinit=False, \
                   maxmgang=deg2rad(20.), \
                   boolproppsfn=False, \
                   indxenerincl=indxenerincl, \
                   indxevttincl=indxevttincl, \
                   mocknumbpnts=array([300]), \
                   maxmnumbpnts=array([600]), \
                   minmflux=minmflux, \
                   maxmflux=maxmflux, \
                   regitype='ngal', \
                   pathdata=os.environ["FERM_NGAL_DATA_PATH"], \
                   strgback=['fermisotflux.fits', 'fermfdfmflux_ngal.fits'], \
                   strgexpo='fermexpo_cmp0_ngal.fits', \
                   datatype='mock', \
                   numbsideheal=256, \
                   mockfluxdistslop=array([2.0]), \
                  )

    
if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

