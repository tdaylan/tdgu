from __init__ import *

def writ_ferm_ngal():

    gdat = tdpy.util.gdatstrt()
    gdat.recotype = ['rec7']
    gdat.enertype = ['pnts']
    gdat.timefrac = [0.1]
    tdpy.util.make_maps_main(gdat, os.environ["FERM_NGAL_DATA_PATH"])
    tdpy.util.prep_maps('rec7', 'pnts', 'ngal', os.environ["FERM_NGAL_DATA_PATH"], 256, 'tim0')
    tdpy.util.prep_maps('rec7', 'pnts', 'ngal', os.environ["FERM_NGAL_DATA_PATH"], 256, 'tim0')


def test_ferm_quas_mock():
    
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
                       minmflux=6e-11, \
                       spatdisttype=spatdisttype, \
                       lgalprio=lgalprio, \
                       bgalprio=bgalprio, \
                       lgalcntr=0., \
                       bgalcntr=pi / 2., \
                       strgexpo='expofermrec8pntsngal0256.fits', \
                      )


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

    
def pcat_ferm_ngal_inpt(strgcnfgextnexec=None):

    dictargs = {}
    dictargs['lgalcntr'] = 0.
    dictargs['bgalcntr'] = pi / 2.
    dictargs['backtype'] = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
    dictargs['strgexpo'] = 'expofermrec8pntsngal0256.fits'
    dictargs['strgexprsbrt'] = 'sbrtferm.fits'
    
    listnamecnfgextn = ['nomi', 'tim1', 'tim2']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['tim1']['strgexprsbrt'] = 'sbrtferm.fits'
    dictargsvari['tim2']['strgexprsbrt'] = 'sbrtferm.fits'
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_ferm_ngal_mock():
     
    pcat.main.init( \
                   lgalcntr=0., \
                   bgalcntr=pi / 2., \
                   #forccart=True, \
                   #pixltype='cart', \
                   #numbsidecart=100, \
                   numbelempop0reg0=100, \
                   back=['fermisotflux.fits', 'sbrtfermfdfmngal.fits'], \
                   strgexpo='expofermrec8pntsngal0256.fits', \
                   numbsideheal=256, \
                  )

globals().get(sys.argv[1])()
