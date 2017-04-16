from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(verbtype=0, elemtype='lens', exprtype='hubb', defa=True)

    numbcnfg = 3
    numbiter = 1
   
    varbinpt = tdpy.util.varb(numbcnfg)
    varbinpt.defn_para('minmdefs', 1e-3, 1e-2, scal='logt')
    varbinpt.defn_para('bacp', 1e-1 * gdat.hubbexpofact / (pi * gdat.truesizesour**2), 1e1 * gdat.hubbexpofact / (pi * gdat.truesizesour**2), scal='logt')
    varbinpt.defn_para('truespecsour', 1e-1 * gdat.hubbexpofact, 1e1 * gdat.hubbexpofact, scal='logt')
    varbinpt.defn_para('truespechost', 1e-1 * gdat.hubbexpofact, 1e1 * gdat.hubbexpofact, scal='logt')

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
                     
                    if varbinpt.name[p] == 'minmdefs':
                        trueminmdefs = varb / gdat.anglfact
                    if varbinpt.name[p] == 'bacp':
                        truebacp = array([varb])
                    if varbinpt.name[p] == 'truespecsour':
                        truespecsour = varb
                    if varbinpt.name[p] == 'truespechost':
                        truespechost = varb
                
                gdat = pcat.main.init( \
                                      seedstat=seedstat, \
                                      numbswep=10000, \
                                      factthin=100, \
                                      makeplot=False, \
                                      elemtype='lens', \
                                      exprtype='hubb', \
                                      truebacp=truebacp, \
                                      trueminmdefs=trueminmdefs, \
                                      truespecsour=truespecsour, \
                                      truespechost=truespechost, \
                                     )
                
                for n in range(numboutpvarb):
                    if listnameoutpvarb[n] == 'dotsassc':
                        grid[0, n, m, l] = gdat.defsfactplot * gdat.medidotsassc[0][0]
                        grid[1:3, n, m, l] = gdat.defsfactplot * gdat.errrdotsassc[0][:, 0]
                        grid[3, n, m, l] = gdat.defsfactplot * gdat.truedots[0][0, 0]

                        if k == 0 and m == 0 and l == 0:
                            liststrgoutpvarb.append(r'$\theta_{E,0}$')
                            listscaloutpvarb.append('logt')
                    else:
                        indx = where(gdat.namefixp == listnameoutpvarb[n])[0]

                        grid[0, n, m, l] = gdat.factfixpplot[indx] * getattr(gdat, 'medifixp')[indx]
                        grid[1:3, n, m, l] = gdat.factfixpplot[indx] * getattr(gdat, 'errrfixp')[:, indx].flatten()
                        grid[3, n, m, l] = gdat.factfixpplot[indx] * getattr(gdat, 'truefixp')[indx]
                        if k == 0 and m == 0 and l == 0:
                            liststrgoutpvarb.append(getattr(gdat, 'strgfixp')[indx][0])
                            listscaloutpvarb.append(getattr(gdat, 'scalfixp')[indx])
                            
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
    
    numbswepnomi = 1000000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['inittype'] = 'refr'
    dictargs['checprio'] = False
    dictargs['condcatl'] = False

    listlablinpt = ['Nominal', '$N=1$', r'$\alpha_{s,min}$', 'Long', 'Penalty']
    dictargsvari = {}
    dictargsvari['numbswep'] = [numbswepnomi, numbswepnomi, numbswepnomi, 3*numbswepnomi, None]
    dictargsvari['fittminmdefs'] = [None, 6e-3/3600./180.*pi, None, None, None]
    dictargsvari['fittminmnumbpnts'] = [None, array([1]), None, None, None]
    dictargsvari['fittmaxmnumbpnts'] = [None, array([1]), None, None, None]
    dictargsvari['priofactdoff'] = [None, None, None, None, 1.]

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
                       numbswep=1000000, \
                       makeplotintr=True, \
                       condcatl=False, \
                       #checprio=True, \
                       elemtype='lens', \
                       inittype='refr', \
                       exprtype='hubb', \
                      )
    

globals().get(sys.argv[1])()
