from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(verbtype=0, elemtype='lens', exprtype='hubb', defa=True)

    numbcnfg = 3
    numbiter = 1
   
    varbinpt = tdpy.util.varb(numbcnfg)
    varbinpt.defn_para('expo', 1.1e3 / gdat.hubbexpofact, 4.4e3 / gdat.hubbexpofact, scal='logt')
    varbinpt.defn_para('bacp', 1e-1 * gdat.hubbexpofact / (pi * gdat.truesizesour**2), 1e1 * gdat.hubbexpofact / (pi * gdat.truesizesour**2), scal='logt')
    varbinpt.defn_para('truespecsour', 1e-1 * gdat.hubbexpofact, 1e1 * gdat.hubbexpofact, scal='logt')

    listnameoutpvarb = ['defsdistsloppop0', 'meanpntspop0', 'dotsassc', 'spechost', 'beinhost']
    numboutpvarb = len(listnameoutpvarb)
    liststrgoutpvarb = []
    listscaloutpvarb = []

    liststrgvarbinpt = [r'$\theta_{E,min}$ [arcsec]', '$A$ [1/cm$^2$/s]', r'$\epsilon$ [cm$^2$ s]', '$f_s$ [erg/cm$^2$/s]', '$f_h$ [erg/cm$^2$/s]']
    grid = empty((4, numboutpvarb, varbinpt.size, numbcnfg))
    
    for k in range(numbiter):
        seedstat = get_state()
        for m in range(varbinpt.size):
            for l in range(varbinpt.numb):
                
                if m > 0 and l == varbinpt.numb / 2:
                    grid[:, :, m, l] = grid[:, :, 0, numbcnfg / 2]
                    continue
    
                for p in range(varbinpt.size):
                    if p == m:
                        varb = varbinpt.para[p][l]
                    else:
                        varb = varbinpt.para[p][numbcnfg / 2]
                     
                    if varbinpt.name[p] == 'expo':
                        expo = varb
                    if varbinpt.name[p] == 'bacp':
                        truebacp = array([varb])
                    if varbinpt.name[p] == 'truespecsour':
                        truespecsour = varb
                
                gdat = pcat.main.init( \
                                      seedstat=seedstat, \
                                      numbswep=200000, \
                                      condcatl=False, \
                                      elemtype='lens', \
                                      exprtype='hubb', \
                                      inittype='refr', \
                                      expo=expo, \
                                      truebacp=truebacp, \
                                      truespecsour=truespecsour, \
                                     )
                
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
    dictargs['inittype'] = 'refr'
    dictargs['checprio'] = False
    dictargs['condcatl'] = False

    listlablinpt = ['Nominal', '$N=1$', r'$\alpha_{s,min}$', 'Long', 'Penalty', 'Prior Check']
    dictargsvari = {}
    dictargsvari['numbswep'] = [numbswepnomi, numbswepnomi, numbswepnomi, 3*numbswepnomi, numbswepnomi, numbswepnomi]
    dictargsvari['fittminmdefs'] = [None, None, 6e-3/3600./180.*pi, None, None, None]
    dictargsvari['inittype'] = ['refr', 'refr', 'pert', 'refr', 'refr', 'refr']
    dictargsvari['fittminmnumbpnts'] = [None, array([1]), None, None, None, None]
    dictargsvari['fittmaxmnumbpnts'] = [None, array([1]), None, None, None, None]
    dictargsvari['priofactdoff'] = [0., 0., 0., 0., 1., 0.]
    dictargsvari['checprio'] = [False, False, False, False, False, True]

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
                   inittype='refr', \
                   exprtype='hubb', \
                  )
    

def pcat_lens_mock():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       numbswep=1000, \
                       factthin=100, \
                       makeplotintr=True, \
                       makeplotfram=False, \
                       condcatl=False, \
                       elemtype='lens', \
                       inittype='refr', \
                       exprtype='hubb', \
                      )
    

globals().get(sys.argv[1])()
