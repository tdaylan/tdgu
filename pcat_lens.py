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
    dictvarb['numbswep'] = 100
    dictvarb['makeplot'] = False
    
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
    

def pcat_lens_intrevalmodlcnts():
   
    pcat.main.init( \
                   elemtype='lens', \
                   makeplotinit=False, \
                   intrevalmodlcnts=True, \
                   truenumbpnts=array([0]), \
                   truemaxmnumbpnts=array([0]), \
                  )
    

def pcat_lens_mock_sing():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       truenumbpnts=array([1]), \
                       trueminmdefs=1e-2/3600./180.*pi, \
                       truemaxmdefs=1e-1/3600./180.*pi, \
                      )
    

def pcat_lens_mock_spmr():
   
    pcat.main.init( \
                   elemtype='lens', \
                   truelgalimps=array([0.]), \
                   truebgalimps=array([0.]), \
                   truedefsimps=array([1e-2 / anglfact]), \
                   truenumbpnts=array([1]), \
                   probtran=1., \
                   probbrde=0., \
                  )
    

def pcat_lens_mock_syst():
   
    seed(4)
    
    numbswepnomi = 2000000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['factthin'] = 200
    dictargs['numbburn'] = 1000000
    #dictargs['makeplot'] = False
    #dictargs['factthin'] = 100
    #dictargs['emptsamp'] = True
    #dictargs['numbproc'] = 2
    #dictargs['initlgalsour'] = -2.41108523236e-07
    #dictargs['initbgalsour'] = 1.2690938483e-07
    #dictargs['initlgalhost'] = 1.1090807602e-07
    #dictargs['initbgalhost'] = 2.26345624113e-08
 
    numbelem = array([25. * 10.**0.9], dtype=int)
    print 'numbelem'
    print numbelem
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['numbswep'] =             [numbswepnomi, numbswepnomi, numbswepnomi,   numbswepnomi,  numbswepnomi, numbswepnomi]
    dictargsvari['truenumbpnts'] =         [None,         None,         None,           numbelem,      None,         None,       ]
    dictargsvari['trueminmdefs'] =         [None,         None,         None,           3e-4/anglfact, None,         None,       ]
    dictargsvari['fittminmdefs'] =         [None,         None,         3e-3/anglfact,  None,          None,         None,       ]
    dictargsvari['fittminmnumbpnts'] =     [None,         array([1]),   None,           None,          None,         None,       ]
    dictargsvari['fittmaxmnumbpnts'] =     [None,         array([1]),   None,           None,          None,         None,       ]
    dictargsvari['truestdvdefsdistslop'] = [0.5,          0.5,          0.5,            0.5,           0.5,          'none',     ]
    dictargsvari['truescalmeanpnts']     = ['self',       'logt',       'self',         'self',        'logt',       'self',     ]
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  indxruns=0, \
                                 )
    

def pcat_lens_mock_perf():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 2000000
    
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['trueminmdefs'] = [2e-3 / anglfact, 4e-3 / anglfact, 8e-3 / anglfact, None,          None,          None]
    dictargsvari['truemaxmdefs'] = [1e-2 / anglfact, 2e-2 / anglfact, 4e-2 / anglfact, None,          None,          None]
    dictargsvari['truebacp']     = [None,            None,            None,            array([3e-8]), array([1e-7]), array([3e-7])]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )
    

def pcat_lens_mock_reln():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['relnindx'] = 1.
    dictargs['liststrgfeatmodu'] = ['lgalbgal']
    dictargs['liststrgpdfnmodu'] = ['tmplreln']
    dictargsvari = {}
    dictargsvari['relnpowr'] = [2., 1., 0.]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_doff():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargsvari = {}
    dictargsvari['priofactdoff'] = [1., 2., 0.5]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_sele():
    
    numbitermacr = 30
    numbiterelem = 10
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['mockonly'] = True
    dictargs['trueminmdefs'] = 5e-4 / anglfact
    dictargs['truenumbpnts'] = array([1000])
    dictargs['variasca'] = False
    dictargs['variacut'] = False
    dictargs['allwfixdtrue'] = False
    dictargs['verbtype'] = 0
    dictargs['makeplot'] = False
    dictargs['truemaxmnumbpnts'] = array([1000])
    
    listnamesele = ['pars', 'nrel']
    numbsele = len(listnamesele)
    listnamefeatsele = ['defs', 'mcut', 'rele']
    numbfeatsele = len(listnamefeatsele)

    dictargsvari = {}
    dictargsvari['elemtype'] = ['lens' for n in range(numbiterelem)]
    
    matrcutf = empty((numbitermacr, numbiterelem, numbfeatsele))
    
    liststrgvarboutp = []
    for strgvarbelem in listnamefeatsele:
        liststrgvarboutp += ['truehist' + strgvarbelem]
        for namesele in listnamesele:
            liststrgvarboutp += ['true' + strgvarbelem + namesele]
            liststrgvarboutp += ['truehist' + strgvarbelem + namesele]
    
    for k in range(numbitermacr):
        listgdat, dictglob = pcat.main.initarry( \
                                      dictargsvari, \
                                      dictargs, \
                                      randseedelem=True, \
                                      liststrgvarboutp=liststrgvarboutp, \
                                     )
        
        gdat = listgdat[0]
        
        if k == 0:
            hist = zeros((numbiterelem, gdat.numbbinsplot))
            histsele = zeros((numbiterelem, gdat.numbbinsplot))
            histfitt = zeros((numbiterelem, gdat.numbbinsplot))
            factsele = zeros((numbiterelem, gdat.numbbinsplot))
            for namesele in listnamesele:
                pathimagsele = gdat.pathimag + 'sele/' + namesele + '/'
                os.system('mkdir -p %s' % pathimagsele)
            truefixp = empty((numbitermacr, gdat.truenumbfixp))
            corrfixpcutf = empty(gdat.truenumbfixp)
            pvalfixpcutf = empty(gdat.truenumbfixp)
        
        for namefixp in gdat.truenamefixp:
            truefixp[k, :] = gdat.truefixp

        for b, namefeat in enumerate(listnamefeatsele):
            lablvarbtotl = getattr(gdat, 'labl' + namefeat + 'totl')
            meanvarb = getattr(gdat, 'mean' + namefeat)
            deltvarb = getattr(gdat, 'delt' + namefeat)
            factvarbplot = getattr(gdat, 'fact' + namefeat + 'plot')
            limtvarb = [factvarbplot * getattr(gdat, 'minm' + namefeat), factvarbplot * getattr(gdat, 'maxm' + namefeat)]
            for m in range(numbiterelem):
                hist[m, :] = dictglob['truehist' + namefeat][m][0, :]
            for namesele in listnamesele:
                pathimagsele = gdat.pathimag + 'sele/' + namesele + '/'
                for m in range(numbiterelem):
                    histsele[m, :] = dictglob['truehist' + namefeat + namesele][m][0, :]
                    factsele[m, :] = histsele[m, :] / hist[m, :]
                    alph, loca, cutf = sp.stats.invgamma.fit(dictglob['true' + namefeat + namesele][m][0], floc=0., f0=1.9)
                    histfitt[m, :] = sum(histsele[m, :]) * sp.stats.invgamma.pdf(meanvarb, alph, loc=loca, scale=cutf) * deltvarb
                    matrcutf[k, m, b] = cutf
            
                meanfactsele = mean(factsele, 0)
                meanmatrcutf = mean(matrcutf, axis=1)

                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                for m in range(numbiterelem):
                    axis.loglog(meanvarb * factvarbplot, hist[m, :], color='g', alpha=0.1, ls='-.')
                    axis.loglog(meanvarb * factvarbplot, histsele[m, :], color='g', alpha=0.1, ls='--')
                    axis.loglog(meanvarb * factvarbplot, histfitt[m, :], color='g', alpha=0.1, ls='-')
                    axis.axvline(matrcutf[k, m, b] * factvarbplot, color='m', alpha=0.1)
                axis.axvline(np.power(prod(matrcutf[k, m, b]), array([1. / numbiterelem])) * factvarbplot, color='m')
    
                axis.loglog(meanvarb * factvarbplot, mean(hist, 0), color='g', ls='-.')
                axis.loglog(meanvarb * factvarbplot, mean(histsele, 0), color='g', ls='--')
                axis.loglog(meanvarb * factvarbplot, mean(histfitt, 0), color='g', ls='-')
    
                axis.set_xlabel(lablvarbtotl)
                axis.set_ylabel('$N$')
                axis.set_ylim([0.5, None])
                axis.set_xlim(limtvarb)
                plt.tight_layout()
                figr.savefig(pathimagsele + 'hist' + namefeat + '%04d.pdf' % k)
                plt.close(figr)
    
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.loglog(meanvarb * factvarbplot, meanfactsele, color='black')
                for m in range(numbiterelem):
                    axis.loglog(meanvarb * factvarbplot, factsele[m, :], alpha=0.1, color='g')
                    axis.axvline(matrcutf[k, m, b] * factvarbplot, color='m', alpha=0.1)
                axis.axvline(np.power(prod(matrcutf[k, :, b]), array([1. / numbiterelem])) * factvarbplot, color='m')
                axis.set_xlabel(lablvarbtotl)
                axis.set_ylabel('$f$')
                axis.set_xlim(limtvarb)
                plt.tight_layout()
                figr.savefig(pathimagsele + 'fact' + namefeat + '%04d.pdf' % k)
                plt.close(figr)
    
    for namesele in listnamesele:
        pathimagsele = gdat.pathimag + 'sele/' + namesele + '/'
        for b, namefeat in enumerate(listnamefeatsele):
            for a, namefixp in enumerate(gdat.truenamefixp):
                corrfixpcutf[a], pvalfixpcutf[a] = sp.stats.stats.pearsonr(truefixp[:, a], meanmatrcutf[:, b])
                print gdat.truelablfixp[a]
                print corrfixpcutf
                print
            indx = where(isfinite(corrfixpcutf) & (pvalfixpcutf < 0.1))[0]
            numb = indx.size
            figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
            for k in range(numb):
                size = 100. * (0.1 - pvalfixpcutf[k]) + 5.
                axis.plot(k + 0.5, corrfixpcutf[indx][k], ls='', marker='o', markersize=size, color='black')
            axis.set_xticks(arange(numb) + 0.5)
            axis.set_xticklabels(gdat.truelablfixp[indx])
            axis.set_xlim([0., numb])
            axis.set_ylabel(r'$\xi$')
            plt.tight_layout()
            figr.savefig(pathimagsele + 'corr' + namefeat + namesele + '.pdf')
            plt.close(figr)
    

def pcat_lens_mock_init():
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 200000
    dictargs['variasca'] = False
    dictargs['variacut'] = False
    dictargs['makeplotfram'] = False
    dictargsvari = {}
    dictargsvari['inittype'] = ['refr', 'pert']

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


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


def pcat_lens_mock_zero():

    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['fittmaxmnumbpnts'] = array([0])
    dictargsvari = {}
    dictargsvari['truenumbpnts'] = [array([0]), array([20])]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_testvari():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['variasca'] = False
    dictargs['numbswep'] = 100000
    dictargs['makeplotfram'] = False
    dictargs['variacut'] = False
    dictargsvari = {}
    dictargsvari['inittype'] = ['refr', 'pert']

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_test():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargsvari = {}
    dictargsvari['inittype'] =         ['refr', 'pert', 'pert', None,  ]
    dictargsvari['variasca'] =         [False,  False,  True,   None,  ]
    dictargsvari['variacut'] =         [False,  False,  True,   None,  ]
    dictargsvari['fittampldisttype'] = [False,  False,  True,  'igam', ]
    dictargsvari['priofactdoff'] =     [   1.,     1.,    1.,      0., ]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_dofftest():
   
    pcat.main.init( \
                   elemtype='lens', \
                   #verbtype=2, \
                   #truenumbpnts=array([5]), \
                   #truemaxmnumbpnts=array([10]), \
                   #checprio=True, \
                   shrtfram=True, \
                   numbswep=20000, \
                   priofactdoff=1., \
                   numbswepplot=5000, \
                   factthin=1000, \
                   makeplotinit=False, \
                  )
   

def pcat_lens_mockonly():
   
    pcat.main.init( \
                   elemtype='lens', \
                   truenumbpnts=array([20]), \
                   truemaxmnumbpnts=array([400]), \
                   mockonly=True, \
                  )


def pcat_lens_mock():
    
    seed(4)
    
    anglfact = 3600. * 180. / pi

    numbiter = 1
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       #verbtype=2, \
                       #truenumbpnts=array([5]), \
                       #truemaxmnumbpnts=array([10]), \
                       #checprio=True, \
                       #evoltype='maxmllik', \
                       #makeplot=False, \
                       #shrtfram=True, \
                       #optihess=False, \
                       #makeplot=False, \
                       #makeplotinit=False, \
                       #mockonly=True, \
                       makeplot=False, \
                       #truemaxmnumbpnts=array([5]), \
                       #truenumbpnts=array([5]), \
                       #verbtype=2, \
                       numbswep=100000, \
                       #truestdvdefsdistslop=0.01, \
                       factthin=30000, \
                       #rtagredo='20170610_133749_pcat_lens_mock_10000', \
                       #inittype='rand', \
                       #optihess=False, \
                       #savestat=True, \
                       #recostat=True, \
                       #initlgalsour=-1e-1 / anglfact, \
                       #initbgalsour=1e-1 / anglfact, \
                       #burntmpr=True, \
                       #fittmaxmnumbpnts=array([5]), \
                       #shrtfram=True, \
                       numbswepplot=6000, \
                      )
   

def pcat_lens_intrevalresicnts():

    anglfact = 3600. * 180. / pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21

    # half-size of the image
    numbside = 400
    maxmgangdata = numbside * 0.5 * sizepixl
    
    # name of the data file
    strgexprflux = namedatasets + '_%04d.fits' % numbside
    
    if namedatasets == 'lens29075550':
        initbacpbac0ene0 = 1.1e-7
        fittmeanbacpbac0ene0 = 1.1e-7
        fittstdvbacpbac0ene0 = fittmeanbacpbac0ene0 * 1e-3
        fittscalbacpbac0ene0 = 'gaus'
    else:
        initbacpbac0ene0 = None
        fittmeanbacpbac0ene0 = None
        fittstdvbacpbac0ene0 = None
        fittscalbacpbac0ene0 = None

    pcat.main.init( \
                   elemtype='lens', \
                   makeplotinit=False, \
                   intrevalresicnts=True, \
                   strgexpo=strgexpo, \
                   initbacpbac0ene0=initbacpbac0ene0, \
                   fittmeanbacpbac0ene0=fittmeanbacpbac0ene0, \
                   fittstdvbacpbac0ene0=fittstdvbacpbac0ene0, \
                   fittscalbacpbac0ene0=fittscalbacpbac0ene0, \
                   recostat='pcat_lens_inpt', \
                   fittmaxmnumbpnts=array([0]), \
                   maxmgangdata=maxmgangdata, \
                   strgexprflux=strgexprflux, \
                  )
    

def pcat_lens_inpt():
    
    anglfact = 3600. * 180. / pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21
    
    # half-size of the image in pixels
    numbside = 400
    maxmgangdata = numbside * 0.5 * sizepixl
    
    # name of the data file
    strgexprflux = namedatasets + '_%04d.fits' % numbside
    
    if namedatasets == 'lens29075550':
        initbacpbac0ene0 = 1.1e-7
        fittmeanbacpbac0ene0 = 1.1e-7
        fittstdvbacpbac0ene0 = fittmeanbacpbac0ene0 * 1e-3
        fittscalbacpbac0ene0 = 'gaus'
    else:
        initbacpbac0ene0 = None
        fittmeanbacpbac0ene0 = None
        fittstdvbacpbac0ene0 = None
        fittscalbacpbac0ene0 = None

    mask = array([-0.3, 0.1, -0.1, 0.2]) / anglfact
    
    pcat.main.init( \
                   elemtype='lens', \
                   numbswep=20000, \
                   numbswepplot=1000, \
                   #shrtfram=True, \
                   #mockonly=True, \
                   #makeplotintr=True, \
                   #burntmpr=True, \
                   #initsizesour=1.5/anglfact, \
                   #initspecsourene0=1.5e-18, \
                   #verbtype=2, \
                   #mask=mask, \
                   initbacpbac0ene0=initbacpbac0ene0, \
                   fittmeanbacpbac0ene0=fittmeanbacpbac0ene0, \
                   fittstdvbacpbac0ene0=fittstdvbacpbac0ene0, \
                   fittscalbacpbac0ene0=fittscalbacpbac0ene0, \
                   serstype='intp', \
                   optihess=False, \
                   savestat=True, \
                   recostat=True, \
                   #makeplotinit=False, \
                   #makeplotfram=False, \
                   makeplotlpri=False, \
                   strgexpo=strgexpo, \
                   fittmaxmnumbpnts=array([0]), \
                   maxmgangdata=maxmgangdata, \
                   strgexprflux=strgexprflux, \
                  )
   

def pcat_lens_psfn():
    
    strgexpo = 7.37487548893e21
    anglfact = 3600. * 180. / pi
    maxmgangdata = 50. * 0.05 / anglfact
    numbiter = 1
    
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       numbswep=50000, \
                       factthin=500, \
                       numbswepplot=10000, \
                       shrtfram=True, \
                       #mockonly=True, \
                       #makeplotintr=True, \
                       #burntmpr=True, \
                       optihess=False, \
                       #savestat=True, \
                       #recostat=True, \
                       #makeplotinit=False, \
                       #makeplotfram=False, \
                       makeplotlpri=False, \
                       strgexpo=strgexpo, \
                       fittmaxmnumbpnts=array([0]), \
                       maxmgangdata=maxmgangdata, \
                       strgexprflux='lens29065407.fits', \
                      )
   

globals().get(sys.argv[1])()
