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
    dictargs['numbswep'] = 1000000
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
    

def pcat_lens_mock_syst():
   
    numbswepnomi = 1000000
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['fittampldisttype'] = 'igam'
                
    dictargsvari = {}
    dictargsvari['numbswep'] =         [numbswepnomi, numbswepnomi,       numbswepnomi, 3*numbswepnomi]
    dictargsvari['fittminmdefs'] =     [None,         6e-3/3600./180.*pi, None,         None          ]
    dictargsvari['fittminmnumbpnts'] = [None,         None,               array([1]),   None          ]
    dictargsvari['fittmaxmnumbpnts'] = [None,         None,               array([1]),   None          ]
    dictargsvari['checprio'] =         [True,         False,              False,        False         ]

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
   
    pcat.main.init( \
                   elemtype='lens', \
                   numbswep=1000, \
                   factthin=100, \
                   makeplot=False, \
                   probbrde=0., \
                  )
    

def pcat_lens_mock_dotn():
  
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 1000000
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
            
            namepdfn, listpara = retr_pdfnscipbest(dictglob['truedefssign'][m][0], gdat.binsdefs, gdat.meandefs)
            print 'namepdfn'
            print namepdfn
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


def retr_pdfnscipbest(data, binsxdat, meanxdat):
    
    print 'binsxdat'
    print binsxdat
    print 'meanxdat'
    print meanxdat

    pdfndata = histogram(data, bins=binsxdat, normed=True)[0]
    st = sp.stats 
    listpdfnscip = [        
        st.alpha,st.anglit,st.arcsine,st.beta,st.betaprime,st.bradford,st.burr,st.cauchy,st.chi,st.chi2,st.cosine,
        st.dgamma,st.dweibull,st.erlang,st.expon,st.exponnorm,st.exponweib,st.exponpow,st.f,st.fatiguelife,st.fisk,
        st.foldcauchy,st.foldnorm,st.frechet_r,st.frechet_l,st.genlogistic,st.genpareto,st.gennorm,st.genexpon,
        st.genextreme,st.gausshyper,st.gamma,st.gengamma,st.genhalflogistic,st.gilbrat,st.gompertz,st.gumbel_r,
        st.gumbel_l,st.halfcauchy,st.halflogistic,st.halfnorm,st.halfgennorm,st.hypsecant,st.invgamma,st.invgauss,
        st.invweibull,st.johnsonsb,st.johnsonsu,st.ksone,st.kstwobign,st.laplace,st.levy,st.levy_l,st.levy_stable,
        st.logistic,st.loggamma,st.loglaplace,st.lognorm,st.lomax,st.maxwell,st.mielke,st.nakagami,st.ncx2,st.ncf,
        st.nct,st.norm,st.pareto,st.pearson3,st.powerlaw,st.powerlognorm,st.powernorm,st.rdist,st.reciprocal,
        st.rayleigh,st.rice,st.recipinvgauss,st.semicircular,st.t,st.triang,st.truncexpon,st.truncnorm,st.tukeylambda,
        st.uniform,st.vonmises,st.vonmises_line,st.wald,st.weibull_min,st.weibull_max,st.wrapcauchy
    ]

    pdfnscipbest = st.norm
    listparabest = (0.0, 1.0)
    chsqbest = inf

    for pdfnscip in listpdfnscip:

        try:
            listpara = pdfnscip.fit(data, floc=0.)
            arg = listpara[:-2]
            loca = listpara[-2]
            scal = listpara[-1]

            pdfnfitt = pdfnscip.pdf(meanxdat, loc=loca, scale=scal, *arg)
            chsq = sum((pdfndata - pdfnfitt)**2)

            if chsqbest > chsq > 0:
                print 'Updating best as ' + pdfnscip.name
                pdfnscipbest = pdfnscip
                listparabest = listpara
                chsqbest = chsq

        except Exception:
            pass

    return pdfnscipbest.name, listparabest


def pcat_lens_mock_comp():
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 250000
    dictargs['numbswepplot'] = 50000
    dictargs['factthin'] = 500
    dictargs['numbburn'] = 200000
    dictargs['verbtype'] = 1
    dictargs['burntmpr'] = True
    dictargs['shrtfram'] = True
    dictargs['trueminmdefs'] = 5e-2 / anglfact
    dictargs['truemaxmdefs'] = 1e-1 / anglfact
    dictargsvari = {}
    dictargsvari['inittype'] = ['pert', 'refr']

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_macr():
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['elemtype'] = 'lens'
    dictargs['numbswep'] = 250000
    dictargs['numbswepplot'] = 50000
    dictargs['factthin'] = 500
    dictargs['numbburn'] = 200000
    dictargs['verbtype'] = 1
    dictargs['burntmpr'] = True
    dictargs['shrtfram'] = True
    dictargs['truemaxmnumbpnts'] = array([0])
    dictargs['truenumbpnts'] = array([0])
    dictargsvari = {}
    dictargsvari['inittype'] = ['pert', 'refr']

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
                       numbswep=1000, \
                       fittampldisttype='igam', \
                       priofactdoff=0., \
                       #verbtype=2, \
                       numbswepplot=10000, \
                       #variasca=False, \
                       #variacut=False, \
                       factthin=100, \
                       #trueminmdefs=1e-2/anglfact, \
                       #truemaxmdefs=5e-2/anglfact, \
                       #makeplot=False, \
                       #fittspatdisttype=['unif'], \
                       #shrtfram=True, \
                       #makeplot=False, \
                       #checprio=True, \
                       #makeplotintr=True, \
                      )
   

globals().get(sys.argv[1])()
