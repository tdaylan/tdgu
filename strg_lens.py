from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(verbtype=0, elemtype='lens', exprtype='hubb', defa=True)

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
    dictvarb['numbswep'] = 400000
    dictvarb['condcatl'] = False
    dictvarb['elemtype'] = 'lens'
    dictvarb['inittype'] = 'pert'
    dictvarb['exprtype'] = 'hubb'
    
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
   
    liststrgvarboutp = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']
    
    numbswepnomi = 10000000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['inittype'] = 'pert'
    dictargs['condcatl'] = False
    dictargs['priofactdoff'] = 0.6
    
    anglfact = 3600. * 180. / pi
    listlablinpt = ['Nominal', '$N=1$', r'$\alpha_{s,min}$', 'No Penalty', 'Long']
    dictargsvari = {}
    dictargsvari['trueminmdefs'] = [1e-3 / anglfact, 5e-3 / anglfact, 2.5e-2 / anglfact, None,          None,          None]
    dictargsvari['truemaxmdefs'] = [1e-2 / anglfact, 5e-2 / anglfact, 2.5e-1 / anglfact, None,          None,          None]
    dictargsvari['truebacp']     = [None,            None,            None,              array([3e-8]), array([1e-7]), array([3e-7])]

    dictglob = pcat.main.initarry( \
                                  liststrgvarboutp, \
                                  dictargsvari, \
                                  dictargs, \
                                  sameseed=True, \
                                  listlablinpt=listlablinpt, \
                                  makeplotarry=True, \
                                 )
    

def pcat_lens_mock_syst():
   
    liststrgvarboutp = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']
    
    numbswepnomi = 1000000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['inittype'] = 'pert'
    dictargs['condcatl'] = False
                
    listlablinpt = ['Nominal', '$N=1$', r'$\alpha_{s,min}$', 'No Penalty', 'Long']
    dictargsvari = {}
    dictargsvari['numbswep'] =         [numbswepnomi, numbswepnomi, numbswepnomi,       numbswepnomi,       numbswepnomi,  numbswepnomi, 5*numbswepnomi]
    dictargsvari['fittminmdefs'] =     [None,         None,         5e-3/3600./180.*pi, None,               None,          None,         None          ]
    dictargsvari['fittminmnumbpnts'] = [None,         None,         None,               array([1]),         array([10]),   None,         None          ]
    dictargsvari['fittmaxmnumbpnts'] = [None,         None,         None,               array([1]),         array([10]),   None,         None          ]
    dictargsvari['priofactdoff'] =     [0.35,         0.5,          0.5,                0.5,                0.5,           0.5,          0.5           ]
    dictargsvari['checprio'] =         [False,        False,        False,              False,              False,         False,        False         ]

    dictglob = pcat.main.initarry( \
                                  liststrgvarboutp, \
                                  dictargsvari, \
                                  dictargs, \
                                  sameseed=True, \
                                  listlablinpt=listlablinpt, \
                                  makeplotarry=True, \
                                 )
    

def pcat_lens_mock_intr():
   
    pcat.main.init( \
                   intreval=True, \
                   elemtype='lens', \
                   inittype='pert', \
                   exprtype='hubb', \
                  )
    

def pcat_lens_mock_sing():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       numbswep=1000000, \
                       inittype='pert', \
                       truenumbpnts=array([1]), \
                       trueminmdefs=5e-2/3600./180.*pi, \
                       truemaxmdefs=1e-1/3600./180.*pi, \
                       condcatl=False, \
                       elemtype='lens', \
                       exprtype='hubb', \
                      )
    

def pcat_lens_mock():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       numbswep=10000, \
                       factthin=100, \
                       inittype='pert', \
                       #checprio=True, \
                       makeplotintr=True, \
                       #shrtfram=True, \
                       condcatl=False, \
                       elemtype='lens', \
                       exprtype='hubb', \
                      )
    

globals().get(sys.argv[1])()
