from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(elemtype='lens', defa=True, verbtype=0)

    numbcnfg = 3
    numbiter = 1
   
    varbinpt = tdpy.util.varb(numbcnfg)
    #varbinpt.defn_para('expo', 1.1e3 / gdat.hubbexpofact, 4.4e3 / gdat.hubbexpofact, scal='logt')
    varbinpt.defn_para('bacpbac0ene0', 3e-8, 3e-7, scal='logt')
    varbinpt.defn_para('specsour', 3e-20, 3e-19, scal='logt')

    listnameoutpvarb = ['defsdistsloppop0', 'meanpntspop0', 'dotsassc', 'spechost', 'beinhost']
    numboutpvarb = len(listnameoutpvarb)
    liststrgoutpvarb = []
    listscaloutpvarb = []

    grid = empty((4, numboutpvarb, varbinpt.size, numbcnfg))
    
    dictvarb = dict()
    dictvarb['elemtype'] = 'lens'
    dictvarb['exprtype'] = 'hubb'
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
                    
                    if varbinpt.name[p] == 'bacpbac0ene0':
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
                   exprtype='hubb', \
                   makeplotinit=False, \
                   intrevalmodlcnts=True, \
                   numbelemreg0pop0=0, \
                   maxmnumbelemreg0pop0=0, \
                  )
    

def pcat_lens_mock_sing():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       exprtype='hubb', \
                       numbelemreg0pop0=1, \
                       minmdefs=1e-2/3600./180.*pi, \
                       maxmdefs=1e-1/3600./180.*pi, \
                      )
    

def pcat_lens_mock_spmr():
   
    pcat.main.init( \
                   elemtype='lens', \
                   exprtype='hubb', \
                   lgalimps=array([0.]), \
                   bgalimps=array([0.]), \
                   defsimps=array([1e-2 / anglfact]), \
                   numbelemreg0pop0=1, \
                   probtran=1., \
                   probbrde=0., \
                  )
    

def pcat_lens_mock_next():
   
    seed(4)
    
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['diagmode'] = True
    dictargs['numbswep'] = 100000
 
    numbelem = array([25. * 10.**0.9], dtype=int)
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['numbelemreg0pop0']     = [None,        0,  0,  25,   int(25. * 0.1**0.9), int(25. * 10.**0.9)]
    dictargsvari['trueminmdefs']     = [None,        None,        None,        3e-3/anglfact, 3e-2/anglfact,                      3e-4/anglfact]
    dictargsvari['fittminmdefs']     = [None,        None,        None,        3e-4/anglfact, 3e-4/anglfact,                      3e-4/anglfact]
    dictargsvari['priofactdoff']     = [0.,          0.,          1.,          1.,            1.,                                 1.]
    dictargsvari['scalmeanpnts'] = ['logt',      'logt',      'logt',      'logt',        'logt',                            'logt']
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  #indxruns=0, \
                                 )
    

def pcat_lens_mock_syst():
   
    seed(4)
    
    numbswepnomi = 100000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['diagmode'] = True
    #dictargs['makeplot'] = False
    #dictargs['verbtype'] = 2
    #dictargs['inittype'] = 'refr'
    #dictargs['diagmode'] = True

    dictargs['truenumbpntsreg0pop0'] = 25
    dictargs['truemeanpntspop0'] = 25
    dictargs['truedefsdistsloppop0'] = 1.9
    dictargs['truesigcene0evt0'] = 4.21788e-07
    dictargs['truebacp0000reg0ene0'] = 2e-07
    dictargs['truelgalsourreg0'] = -2.41109e-07
    dictargs['truebgalsourreg0'] = 1.26909e-07
    dictargs['truefluxsourreg0'] = 1e-18
    dictargs['truesizesourreg0'] = 1.45444e-06
    dictargs['trueellpsourreg0'] = 0.2
    dictargs['trueanglsourreg0'] = 2.4485
    dictargs['truelgalhostreg0'] = 1.10908e-07
    dictargs['truebgalhostreg0'] = 2.26346e-08
    dictargs['truefluxhostreg0'] = 1e-16
    dictargs['truesizehostreg0'] = 4.84814e-06
    dictargs['truebeinhostreg0'] = 7.27221e-06
    dictargs['trueellphostreg0'] = 0.2
    dictargs['trueanglhostreg0'] = 1.21445
    dictargs['trueserihostreg0'] = 4
    dictargs['truesherextrreg0'] = 0.0956653
    dictargs['truesangextrreg0'] = 1.5708

    dictargs['truelgalreg0pop00000'] = 8.70681e-06
    dictargs['truebgalreg0pop00000'] = 5.5522e-06
    dictargs['truedefsreg0pop00000'] = 4.46996e-07
    dictargs['trueascareg0pop00000'] = 8.3953e-08
    dictargs['trueacutreg0pop00000'] = 7.26722e-07
    dictargs['truelgalreg0pop00001'] = 1.95366e-06
    dictargs['truebgalreg0pop00001'] = -6.43887e-06
    dictargs['truedefsreg0pop00001'] = 2.0933e-07
    dictargs['trueascareg0pop00001'] = 1.98019e-07
    dictargs['trueacutreg0pop00001'] = 5.11875e-06
    dictargs['truelgalreg0pop00002'] = 8.48563e-06
    dictargs['truebgalreg0pop00002'] = 4.20743e-07
    dictargs['truedefsreg0pop00002'] = 5.50444e-08
    dictargs['trueascareg0pop00002'] = 7.67089e-08
    dictargs['trueacutreg0pop00002'] = 5.28643e-06
    dictargs['truelgalreg0pop00003'] = 4.73257e-07
    dictargs['truebgalreg0pop00003'] = 2.66861e-06
    dictargs['truedefsreg0pop00003'] = 8.56312e-08
    dictargs['trueascareg0pop00003'] = 3.15034e-07
    dictargs['trueacutreg0pop00003'] = 3.84845e-06
    dictargs['truelgalreg0pop00004'] = 2.40305e-06
    dictargs['truebgalreg0pop00004'] = 5.18566e-06
    dictargs['truedefsreg0pop00004'] = 6.03287e-08
    dictargs['trueascareg0pop00004'] = 1.82084e-07
    dictargs['trueacutreg0pop00004'] = 4.8727e-06
    dictargs['truelgalreg0pop00005'] = 3.61995e-06
    dictargs['truebgalreg0pop00005'] = -4.77678e-06
    dictargs['truedefsreg0pop00005'] = 1.18797e-07
    dictargs['trueascareg0pop00005'] = 3.02975e-07
    dictargs['trueacutreg0pop00005'] = 8.68302e-06
    dictargs['truelgalreg0pop00006'] = -2.65962e-06
    dictargs['truebgalreg0pop00006'] = 2.66758e-06
    dictargs['truedefsreg0pop00006'] = 6.1361e-08
    dictargs['trueascareg0pop00006'] = 2.41337e-07
    dictargs['trueacutreg0pop00006'] = 1.76904e-06
    dictargs['truelgalreg0pop00007'] = 8.11351e-06
    dictargs['truebgalreg0pop00007'] = -1.32214e-06
    dictargs['truedefsreg0pop00007'] = 3.43939e-07
    dictargs['trueascareg0pop00007'] = 2.02059e-07
    dictargs['trueacutreg0pop00007'] = 8.7719e-06
    dictargs['truelgalreg0pop00008'] = -1.84568e-06
    dictargs['truebgalreg0pop00008'] = -3.27396e-06
    dictargs['truedefsreg0pop00008'] = 1.24152e-07
    dictargs['trueascareg0pop00008'] = 4.09883e-07
    dictargs['trueacutreg0pop00008'] = 8.34863e-06
    dictargs['truelgalreg0pop00009'] = 1.85564e-06
    dictargs['truebgalreg0pop00009'] = -8.05447e-06
    dictargs['truedefsreg0pop00009'] = 1.32745e-07
    dictargs['trueascareg0pop00009'] = 1.18999e-07
    dictargs['trueacutreg0pop00009'] = 7.10343e-06
    dictargs['truelgalreg0pop00010'] = 7.65329e-06
    dictargs['truebgalreg0pop00010'] = 2.85729e-07
    dictargs['truedefsreg0pop00010'] = 1.35078e-07
    dictargs['trueascareg0pop00010'] = 3.15458e-08
    dictargs['trueacutreg0pop00010'] = 5.23671e-06
    dictargs['truelgalreg0pop00011'] = -7.19101e-06
    dictargs['truebgalreg0pop00011'] = 2.22167e-06
    dictargs['truedefsreg0pop00011'] = 8.00093e-08
    dictargs['trueascareg0pop00011'] = 3.7222e-07
    dictargs['trueacutreg0pop00011'] =  4.706e-07
    dictargs['truelgalreg0pop00012'] = -7.56662e-06
    dictargs['truebgalreg0pop00012'] = 3.56868e-06
    dictargs['truedefsreg0pop00012'] = 1.07991e-07
    dictargs['trueascareg0pop00012'] = 2.7714e-07
    dictargs['trueacutreg0pop00012'] = 8.18081e-06
    dictargs['truelgalreg0pop00013'] = -2.37798e-07
    dictargs['truebgalreg0pop00013'] = 6.01449e-06
    dictargs['truedefsreg0pop00013'] = 1.06915e-07
    dictargs['trueascareg0pop00013'] = 4.49287e-07
    dictargs['trueacutreg0pop00013'] = 6.46671e-06
    dictargs['truelgalreg0pop00014'] = -6.81208e-06
    dictargs['truebgalreg0pop00014'] = -2.62666e-06
    dictargs['truedefsreg0pop00014'] = 4.45121e-07
    dictargs['trueascareg0pop00014'] = 1.69823e-07
    dictargs['trueacutreg0pop00014'] = 1.83285e-06
    dictargs['truelgalreg0pop00015'] = -5.30943e-07
    dictargs['truebgalreg0pop00015'] = -2.07925e-06
    dictargs['truedefsreg0pop00015'] = 1.41112e-07
    dictargs['trueascareg0pop00015'] = 2.1175e-07
    dictargs['trueacutreg0pop00015'] = 2.52997e-06
    dictargs['truelgalreg0pop00016'] = -1.69739e-06
    dictargs['truebgalreg0pop00016'] = -1.57014e-06
    dictargs['truedefsreg0pop00016'] = 6.30512e-07
    dictargs['trueascareg0pop00016'] = 4.74931e-07
    dictargs['trueacutreg0pop00016'] = 6.04629e-06
    dictargs['truelgalreg0pop00017'] = -8.08312e-06
    dictargs['truebgalreg0pop00017'] = 4.51844e-06
    dictargs['truedefsreg0pop00017'] = 1.70373e-07
    dictargs['trueascareg0pop00017'] = 4.00467e-07
    dictargs['trueacutreg0pop00017'] = 3.36898e-06
    dictargs['truelgalreg0pop00018'] = -8.55444e-06
    dictargs['truebgalreg0pop00018'] = 2.16851e-06
    dictargs['truedefsreg0pop00018'] = 5.61476e-08
    dictargs['trueascareg0pop00018'] = 3.6823e-07
    dictargs['trueacutreg0pop00018'] = 7.70295e-06
    dictargs['truelgalreg0pop00019'] = -1.77196e-06
    dictargs['truebgalreg0pop00019'] = 8.60636e-06
    dictargs['truedefsreg0pop00019'] = 5.99085e-08
    dictargs['trueascareg0pop00019'] = 4.56979e-07
    dictargs['trueacutreg0pop00019'] = 4.51352e-06
    dictargs['truelgalreg0pop00020'] = 5.10012e-06
    dictargs['truebgalreg0pop00020'] = 4.90283e-06
    dictargs['truedefsreg0pop00020'] = 3.59422e-06
    dictargs['trueascareg0pop00020'] = 4.71858e-07
    dictargs['trueacutreg0pop00020'] = 3.88735e-07
    dictargs['truelgalreg0pop00021'] = -9.38655e-06
    dictargs['truebgalreg0pop00021'] = 8.39509e-07
    dictargs['truedefsreg0pop00021'] =  5.037e-08
    dictargs['trueascareg0pop00021'] = 3.18758e-07
    dictargs['trueacutreg0pop00021'] = 5.18656e-06
    dictargs['truelgalreg0pop00022'] = 3.22999e-06
    dictargs['truebgalreg0pop00022'] = -1.13755e-06
    dictargs['truedefsreg0pop00022'] = 6.89786e-08
    dictargs['trueascareg0pop00022'] = 4.60412e-07
    dictargs['trueacutreg0pop00022'] = 7.0615e-06
    dictargs['truelgalreg0pop00023'] = -9.57364e-06
    dictargs['truebgalreg0pop00023'] = -7.77006e-06
    dictargs['truedefsreg0pop00023'] = 1.46526e-07
    dictargs['trueascareg0pop00023'] = 1.39644e-07
    dictargs['trueacutreg0pop00023'] = 3.95241e-06
    dictargs['truelgalreg0pop00024'] = 5.45961e-06
    dictargs['truebgalreg0pop00024'] = -2.82849e-06
    dictargs['truedefsreg0pop00024'] = 8.48926e-07
    dictargs['trueascareg0pop00024'] = 3.49285e-07
    dictargs['trueacutreg0pop00024'] = 5.35163e-06
    dictargs['verbtype'] = 2

    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['numbswep']                 = [numbswepnomi, numbswepnomi, numbswepnomi,  numbswepnomi,  numbswepnomi, numbswepnomi]
    dictargsvari['truenumbelemreg0pop0']     = [None,         None,         None,          numbelem,      None,         None,       ]
    dictargsvari['trueminmdefs']             = [None,         None,         None,          3e-4/anglfact, None,         None,       ]
    dictargsvari['fittminmdefs']             = [None,         None,         3e-3/anglfact, None,          None,         None,       ]
    dictargsvari['fittminmnumbelemreg0pop0'] = [None,         1,            None,          None,          None,         None,       ]
    dictargsvari['fittmaxmnumbelemreg0pop0'] = [None,         1,            None,          None,          None,         None,       ]
    dictargsvari['stdvdefsdistslop']         = [0.5,          0.5,          0.5,           0.5,           0.5,          'none',     ]
    dictargsvari['scalmeanpnts']             = ['self',       'logt',       'self',        'self',        'logt',       'self',     ]
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  indxruns=1, \
                                 )
    

def pcat_lens_mock_perf():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
    dictargs['numbswep'] = 2000000
    
    anglfact = 3600. * 180. / pi
    dictargsvari = {}
    dictargsvari['minmdefs'] = [2e-3 / anglfact, 4e-3 / anglfact, 8e-3 / anglfact, None,          None,          None]
    dictargsvari['maxmdefs'] = [1e-2 / anglfact, 2e-2 / anglfact, 4e-2 / anglfact, None,          None,          None]
    dictargsvari['bacp']     = [None,            None,            None,            array([3e-8]), array([1e-7]), array([3e-7])]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )
    

def pcat_lens_mock_reln():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
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
    dictargs['exprtype'] = 'hubb'
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
    dictargs['exprtype'] = 'hubb'
    dictargs['minmdefs'] = 5e-4 / anglfact
    dictargs['numbelemreg0pop0'] = 1000
    dictargs['variasca'] = False
    dictargs['variacut'] = False
    dictargs['allwfixdtrue'] = False
    dictargs['verbtype'] = 0
    dictargs['makeplot'] = False
    dictargs['maxmnumbelemreg0pop0'] = 1000
    
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
    dictargs['exprtype'] = 'hubb'
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
    dictargs['exprtype'] = 'hubb'
    dictargs['burntmpr'] = True
    dictargs['maxmnumbelemreg0pop0'] = 0
    dictargs['numbelemreg0pop0'] = 0
    dictargsvari = {}
    dictargsvari['inittype'] = ['refr', 'rand']

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_macr():
    
    pcat.main.init( \
                   elemtype='lens', \
                   exprtype='hubb', \
                   maxmnumbelemreg0pop0=0, \
                   numbelemreg0pop0=0, \
                  )


def pcat_lens_mock_testvari():
   
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['exprtype'] = 'hubb'
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
    dictargs['exprtype'] = 'hubb'
    dictargsvari = {}
    dictargsvari['inittype'] =         ['refr', 'pert', 'pert', None,  ]
    dictargsvari['variasca'] =         [False,  False,  True,   None,  ]
    dictargsvari['variacut'] =         [False,  False,  True,   None,  ]
    dictargsvari['ampldisttype'] = [False,  False,  True,  'igam', ]
    dictargsvari['priofactdoff'] =     [   1.,     1.,    1.,      0., ]

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_dofftest():
   
    pcat.main.init( \
                   elemtype='lens', \
                   exprtype='hubb', \
                   #verbtype=2, \
                   #checprio=True, \
                   shrtfram=True, \
                   numbswep=20000, \
                   priofactdoff=1., \
                   numbswepplot=5000, \
                   factthin=1000, \
                   makeplotinit=False, \
                  )
   

def pcat_lens_mock_many():
   
    pcat.main.init( \
                   elemtype='lens', \
                   exprtype='hubb', \
                   numbswep=100000, \
                   factthin=1000, \
                   numbburn=0, \
                   inittype='refr', \
                   numbelemreg0pop0=25, \
                   maxmnumbelemreg0pop0=100, \
                   numbelemreg1pop0=25, \
                   maxmnumbelemreg1pop0=100, \
                   #makeplot=False, \
                   #explprop=True, \
                   #verbtype=0, \
                   numbregi=2, \
                   #mockonly=True, \
                  )


def pcat_lens_mockonly():
   
    pcat.main.init( \
                   elemtype='lens', \
                   exprtype='hubb', \
                   numbelemreg0pop0=20, \
                   maxmnumbelemreg0pop0=400, \
                   mockonly=True, \
                  )


def pcat_lens_mock():
    
    seed(4)
    
    anglfact = 3600. * 180. / pi

    numbiter = 1
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       exprtype='hubb', \
                       inittype='refr', \
                       numbswep=100000, \
                       #numbswepplot=6000, \
                       verbtype=2, \
                       maxmnumbelemreg0pop0=10, \
                       numbelemreg0pop0=5, \
                       #shrtfram=True, \
                       #optihess=False, \
                       #makeplot=False, \
                       #makeplotinit=False, \
                       #mockonly=True, \
                       #makeplot=False, \
                       #verbtype=2, \
                       #stdvdefsdistslop=0.01, \
                       #rtagredo='20170610_133749_pcat_lens_mock_10000', \
                       #inittype='rand', \
                       #optihess=False, \
                       #savestat=True, \
                       #inittype='reco', \
                       #initlgalsour=-1e-1 / anglfact, \
                       #initbgalsour=1e-1 / anglfact, \
                       #burntmpr=True, \
                       #shrtfram=True, \
                      )
   

globals().get(sys.argv[1])()
