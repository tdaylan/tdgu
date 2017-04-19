from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(verbtype=0, elemtype='lens', exprtype='hubb', defa=True)

    numbcnfg = 3
    numbiter = 1
   
    varbinpt = tdpy.util.varb(numbcnfg)
    #varbinpt.defn_para('expo', 1.1e3 / gdat.hubbexpofact, 4.4e3 / gdat.hubbexpofact, scal='logt')
    varbinpt.defn_para('truebacpbac0ene0', 3e-8, 3e-7, scal='logt')
    #varbinpt.defn_para('truespecsour', 1e-1 * gdat.hubbexpofact, 1e1 * gdat.hubbexpofact, scal='logt')

    listnameoutpvarb = ['defsdistsloppop0', 'meanpntspop0', 'dotsassc', 'spechost', 'beinhost']
    numboutpvarb = len(listnameoutpvarb)
    liststrgoutpvarb = []
    listscaloutpvarb = []

    liststrgvarbinpt = [r'$\theta_{E,min}$ [arcsec]', '$A$ [1/cm$^2$/s]', r'$\epsilon$ [cm$^2$ s]', '$f_s$ [erg/cm$^2$/s]', '$f_h$ [erg/cm$^2$/s]']
    grid = empty((4, numboutpvarb, varbinpt.size, numbcnfg))
    
    dictvarb = dict()
    dictvarb['numbswep'] = 200000
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
    

def pcat_lens_mock_syst():
   
    liststrgvarboutp = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']
    
    numbswepnomi = 200000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['inittype'] = 'pert'
    dictargs['checprio'] = False
    dictargs['condcatl'] = False
                
    namecnfg = 'pcat_lens_mock_syst'

    listlablinpt = ['Nominal', '$N=1$', r'$\alpha_{s,min}$', 'No Penalty', 'Long']
    dictargsvari = {}
    dictargsvari['numbswep'] =         [numbswepnomi,      numbswepnomi,      numbswepnomi,       numbswepnomi,      30*numbswepnomi  ]
    dictargsvari['fittminmdefs'] =     [None,              None,              1e-4/3600./180.*pi, None,              None             ]
    dictargsvari['fittminmnumbpnts'] = [None,              array([1]),        None,               None,              None             ]
    dictargsvari['fittmaxmnumbpnts'] = [None,              array([1]),        None,               None,              None             ]
    dictargsvari['priofactdoff'] =     [1.,                1.,                1.,                 0.,                1.               ]
    dictargsvari['checprio'] =         [True,              False,             False,              False,             False            ]
    dictargsvari['strgcnfg'] =         [namecnfg + '0000', namecnfg + '0001', namecnfg + '0002',  namecnfg + '0003', namecnfg + '0004']

    dictglob = pcat.main.initarry( \
                                  liststrgvarboutp, \
                                  dictargsvari, \
                                  dictargs, \
                                  sameseed=True, \
                                  listlablinpt=listlablinpt, \
                                  makeplotarry=True, \
                                 )
    
    print 'dictglob'
    print dictglob


def pcat_lens_mock_intr():
   
    pcat.main.init( \
                   intreval=True, \
                   elemtype='lens', \
                   inittype='pert', \
                   exprtype='hubb', \
                  )
    

def pcat_lens_mock():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       numbswep=1000, \
                       factthin=100, \
                       inittype='pert', \
                       checprio=True, \
                       makeplotintr=True, \
                       #shrtfram=True, \
                       #makeplotfram=False, \
                       condcatl=False, \
                       elemtype='lens', \
                       exprtype='hubb', \
                      )
    

globals().get(sys.argv[1])()
