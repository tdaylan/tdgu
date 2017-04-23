from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(elemtype='lens', defa=True, verbtype=0)

    numbcnfg = 3
    numbiter = 1
   
    varbinpt = tdpy.util.varb(numbcnfg)
    #varbinpt.defn_para('expo', 1.1e3 / gdat.hubbexpofact, 4.4e3 / gdat.hubbexpofact, scal='logt')
    varbinpt.defn_para('truebacpbac0ene0', 3e-8, 3e-7, scal='logt')
    varbinpt.defn_para('truespecsour', 3e-20, 3e-19, scal='logt')

    listnameoutpvarb = ['defsdistsloppop0', 'meanpntspop0', 'dotsassc', 'spechost', 'beinhost']
    numboutpvarb = len(listnameoutpvarb)
    liststrgoutpvarb = []
    listscaloutpvarb = []

    liststrgvarbinpt = [r'$\theta_{E,min}$ [arcsec]', '$A$ [1/cm$^2$/s]', r'$\epsilon$ [cm$^2$ s]', '$f_s$ [erg/cm$^2$/s]', '$f_h$ [erg/cm$^2$/s]']
    grid = empty((4, numboutpvarb, varbinpt.size, numbcnfg))
    
    dictvarb = dict()
    dictvarb['elemtype'] = 'lens'
    
    cntrcnfg = 0
    for k in range(numbiter):
        dictvarb['seedstat'] = get_state()
        for m in range(varbinpt.size):
            for l in range(varbinpt.numb):
                
                if m > 0 and l == varbinpt.numb / 2:
                    grid[:, :, m, l] = grid[:, :, 0, numbcnfg / 2]
                    continue
    
                for p in range(varbinpt.size):
                    if p == m:
                        dictvarb[varbinpt.name[p]] = varbinpt.para[p][l]
                    else:
                        dictvarb[varbinpt.name[p]] = varbinpt.para[p][numbcnfg / 2]
                    
                    if varbinpt.name[p] == 'truebacpbac0ene0':
                        dictvarb[varbinpt.name[p]] = array([dictvarb[varbinpt.name[p]]])
                
                dictvarb['strgcnfg'] = 'pcat_lens_mock_grid_%04d' % cntrcnfg
                
                gdat = pcat.main.init(**dictvarb)
                cntrcnfg += 1

                #for n in range(numboutpvarb):
                #    if listnameoutpvarb[n] == 'dotsassc':
                #        grid[0, n, m, l] = gdat.anglfact * gdat.medidotsassc[0][0]
                #        grid[1:3, n, m, l] = gdat.anglfact * gdat.errrdotsassc[0][:, 0]
                #        grid[3, n, m, l] = gdat.anglfact * gdat.truedots[0][0, 0]

                #        if k == 0 and m == 0 and l == 0:
                #            liststrgoutpvarb.append(r'$\theta_{E,0}$')
                #            listscaloutpvarb.append('logt')
                #    else:
                #        indx = where(gdat.fittnamefixp == listnameoutpvarb[n])[0]
                #        print 'indx'
                #        print indx
                #        print 'gdat.fittfactfixpplot[indx]'
                #        print gdat.fittfactfixpplot[indx]
                #        print 'getattr(gdat, medifixp)[indx]'
                #        print getattr(gdat, 'medifixp')[indx]
                #        print
                #        grid[0, n, m, l] = gdat.fittfactfixpplot[indx] * getattr(gdat, 'medifixp')[indx]
                #        grid[1:3, n, m, l] = gdat.fittfactfixpplot[indx] * getattr(gdat, 'errrfixp')[:, indx].flatten()
                #        grid[3, n, m, l] = gdat.fittfactfixpplot[indx] * getattr(gdat, 'truefixp')[indx]
                #        if k == 0 and m == 0 and l == 0:
                #            liststrgoutpvarb.append(gdat.fittlablfixp[indx][0])
                #            listscaloutpvarb.append(gdat.fittscalfixp[indx])
                            
        path = os.environ["PCAT_DATA_PATH"] + '/imag/%s_lensgrid/' % gdat.strgtimestmp
        os.system('mkdir -p %s' % path)
        for n in range(numboutpvarb):
            for m in range(varbinpt.size):
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.errorbar(varbinpt.para[m], grid[0, n, m, :], yerr=grid[1:3, n, m, :], ls='', marker='o')
                axis.plot(varbinpt.para[m], grid[3, n, m, :], marker='x', color='g')
                axis.set_xlabel(liststrgvarbinpt[m])
                axis.set_ylabel(liststrgoutpvarb[n])
                maxm = amax(varbinpt.para[m])
                minm = amin(varbinpt.para[m])
                if varbinpt.scal[m] == 'logt':
                    axis.set_xscale('log')
                    axis.set_xlim([minm / 2., maxm * 2.])
                else:
                    axis.set_xlim([minm - 1., maxm + 1.])
                plt.tight_layout()
                plt.savefig('%s/%s%d.pdf' % (path, listnameoutpvarb[n], m))
                plt.close(figr)
    

def pcat_lens_mock_perf():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['trueminmdefs'] = [1e-3 / anglfact, 5e-3 / anglfact, 2.5e-2 / anglfact, None,          None,          None]
    dictargsvari['truemaxmdefs'] = [1e-2 / anglfact, 5e-2 / anglfact, 2.5e-1 / anglfact, None,          None,          None]
    dictargsvari['truebacp']     = [None,            None,            None,              array([3e-8]), array([1e-7]), array([3e-7])]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )
    

def pcat_lens_mock_syst():
   
    numbswepnomi = 1000000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
                
    dictargsvari = {}
    dictargsvari['numbswep'] =         [numbswepnomi, numbswepnomi, numbswepnomi, numbswepnomi,       numbswepnomi, numbswepnomi,  5*numbswepnomi]
    dictargsvari['fittminmdefs'] =     [None,         None,         None,         6e-3/3600./180.*pi, None,         None,          None          ]
    dictargsvari['fittminmnumbpnts'] = [None,         None,         None,         None,               array([1]),   array([10]),   None          ]
    dictargsvari['fittmaxmnumbpnts'] = [None,         None,         None,         None,               array([1]),   array([10]),   None          ]
    dictargsvari['priofactdoff'] =     [0.  ,         1.,           1.,           0.,                 0.,           0.,            0.            ]
    dictargsvari['checprio'] =         [True,         False,        False,        False,              False,        False,         False         ]
    dictargsvari['dotnpowr'] =         [0.,           0.,           4.,           0.,                 0.,           0.,            0.            ]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )
    

def pcat_lens_mock_intr():
   
    pcat.main.init( \
                   elemtype='lens', \
                   intreval=True, \
                  )
    

def pcat_lens_mock_sing():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       truenumbpnts=array([1]), \
                       trueminmdefs=5e-2/3600./180.*pi, \
                       truemaxmdefs=1e-1/3600./180.*pi, \
                      )
    

def pcat_lens_mock_spmr():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       numbswep=1000, \
                       factthin=100, \
                       makeplot=False, \
                       probbrde=0., \
                       #checprio=True, \
                       makeplotintr=True, \
                       #shrtfram=True, \
                      )
    

def pcat_lens_mock_dotn():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargsvari = {}
    dictargsvari['dotnpowr'] = [0., 5., 10.]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_doff():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargsvari = {}
    dictargsvari['priofactdoff'] = [0., 0.33, 0.66, 1.]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock():
   
    anglfact = 3600. * 180. / pi
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       numbswep=100000, \
                       trueminmdefs=1e-2/anglfact, \
                       truemaxmdefs=5e-2/anglfact, \
                       fittspatdisttype=['unif'], \
                       priofactdoff=1., \
                       #shrtfram=True, \
                       #makeplot=False, \
                       #checprio=True, \
                       #makeplotintr=True, \
                      )
   

globals().get(sys.argv[1])()
