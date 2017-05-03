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
                       numbswep=500000, \
                       truenumbpnts=array([1]), \
                       trueminmdefs=1e-2/3600./180.*pi, \
                       truemaxmdefs=1e-1/3600./180.*pi, \
                      )
    

def pcat_lens_mock_spmr():
   
    pcat.main.init( \
                   elemtype='lens', \
                   numbswep=500000, \
                   truelgalimps=array([0.]), \
                   truebgalimps=array([0.]), \
                   truedefsimps=array([1e-2 / anglfact]), \
                   truenumbpnts=array([1]), \
                   probtran=1., \
                   probbrde=0., \
                  )
    

def pcat_lens_mock_syst():
   
    numbswepnomi = 500000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['fittampldisttype'] = 'igam'
                
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['numbswep'] =         [numbswepnomi, numbswepnomi,       numbswepnomi, 3*numbswepnomi]
    dictargsvari['truenumbpnts'] =     [None,         320,                None,         None          ]
    dictargsvari['trueminmdefs'] =     [None,         5e-4/anglfact,      None,         None          ]
    dictargsvari['fittminmdefs'] =     [None,         2e-3/anglfact,      None,         None          ]
    dictargsvari['fittminmnumbpnts'] = [None,         None,               array([1]),   None          ]
    dictargsvari['fittmaxmnumbpnts'] = [None,         None,               array([1]),   None          ]
    dictargsvari['checprio'] =         [True,         False,              False,        False         ]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )
    

def pcat_lens_mock_dotn():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 500000
    dictargs['fittampldisttype'] = 'igam'
    dictargsvari = {}
    dictargsvari['dotnpowr'] = [2., 1., 0.]
    dictargsvari['spatdisttype'] = [['grad'], ['grad'], ['unif']]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_doff():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 500000
    dictargsvari = {}
    dictargsvari['priofactdoff'] = [1., 2., 0.5]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_perf():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 500000
    dictargs['fittampldisttype'] = 'igam'
    
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['trueminmdefs'] = [1e-3 / anglfact, 2e-3 / anglfact, 4e-3 / anglfact, None,          None,          None]
    dictargsvari['truemaxmdefs'] = [1e-2 / anglfact, 2e-2 / anglfact, 4e-2 / anglfact, None,          None,          None]
    dictargsvari['truebacp']     = [None,            None,            None,            array([3e-8]), array([1e-7]), array([3e-7])]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )
    

def pcat_lens_mock_sign():
    
    numbiter = 3

    dictargs = {}
    dictargs['mockonly'] = True
    dictargs['trueminmdefs'] = 5e-4 / 3600. / 180. * pi
    dictargs['truenumbpnts'] = array([1000])
    dictargs['makeplot'] = False
    dictargs['truemaxmnumbpnts'] = array([1000])
    
    liststrgvarboutp = ['truehistdefssign', 'truehistdefs', 'truedefssign', 'truedefs', 'trueindxelemsign']
    
    dictargsvari = {}
    dictargsvari['elemtype'] = ['lens' for n in range(5)]
    
    numbinpt = len(dictargsvari['elemtype'])

    for k in range(numbiter):
        listgdat, dictglob = pcat.main.initarry( \
                                      dictargsvari, \
                                      dictargs, \
                                      randseedelem=True, \
                                      liststrgvarboutp=liststrgvarboutp, \
                                     )
        
        gdat = listgdat[0]

        hist = zeros((numbinpt, gdat.numbbinsplot))
        histsign = zeros((numbinpt, gdat.numbbinsplot))
        histfitt = zeros((numbinpt, gdat.numbbinsplot))
        factsign = zeros((numbinpt, gdat.numbbinsplot))
        for m in range(numbinpt):
            hist[m, :] = dictglob['truehistdefs'][m][0, :]
            histsign[m, :] = dictglob['truehistdefssign'][m][0, :]
            factsign[m, :] = histsign[m, :] / hist[m, :]
            alph, loca, beta = sp.stats.invgamma.fit(dictglob['truedefssign'][m][0], floc=0.)
            histfitt[m, :] = len(dictglob['trueindxelemsign'][m][0]) * sp.stats.invgamma.pdf(gdat.meandefs, alph, loc=loca, scale=beta) * gdat.deltdefs
        
        meanfactsign = mean(factsign, 0)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        for m in range(numbinpt):
            axis.loglog(gdat.meandefs * gdat.factdefsplot, hist[m, :], color='g', alpha=0.1, ls='-')
            axis.loglog(gdat.meandefs * gdat.factdefsplot, histsign[m, :], color='g', alpha=0.1, ls='--')
            axis.loglog(gdat.meandefs * gdat.factdefsplot, histfitt[m, :], color='g', alpha=0.1, ls='-.')
        axis.loglog(gdat.meandefs * gdat.factdefsplot, mean(hist, 0), color='g', ls='-')
        axis.loglog(gdat.meandefs * gdat.factdefsplot, mean(histsign, 0), color='g', ls='--')
        axis.loglog(gdat.meandefs * gdat.factdefsplot, mean(histfitt, 0), color='g', ls='-.')

        axis.set_xlabel(gdat.labldefstotl)
        axis.set_ylabel('$N$')
        axis.set_ylim([0.5, None])
        plt.tight_layout()
        figr.savefig(gdat.pathimag + 'histdefs%04d.pdf' % k)
        plt.close(figr)

        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.loglog(gdat.meandefs * gdat.factdefsplot, meanfactsign, color='black')
        for m in range(numbinpt):
            axis.loglog(gdat.meandefs * gdat.factdefsplot, factsign[m, :], alpha=0.1, color='g')
        axis.set_xlabel(gdat.labldefstotl)
        axis.set_ylabel('$f$')
        plt.tight_layout()
        figr.savefig(gdat.pathimag + 'factsign%04d.pdf' % k)
        plt.close(figr)


def pcat_lens_mock_tmpr():
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['burntmpr'] = True
    dictargs['truemaxmnumbpnts'] = array([0])
    dictargs['truenumbpnts'] = array([0])
    dictargsvari = {}
    dictargsvari['inittype'] = ['refr', 'rand']

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_macr():
    
    pcat.main.init( \
                   elemtype='lens', \
                   truemaxmnumbpnts=array([0]), \
                   truenumbpnts=array([0]), \
                  )


def pcat_lens_mock_test():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargsvari = {}
    dictargsvari['inittype'] = ['refr', 'pert', 'pert']
    dictargsvari['variasca'] = [False,  False,  True, ]
    dictargsvari['variacut'] = [False,  False,  True, ]

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
                       fittampldisttype='igam', \
                       variasca=False, \
                       variacut=False, \
                       inittype='refr', \
                      )
   

globals().get(sys.argv[1])()
