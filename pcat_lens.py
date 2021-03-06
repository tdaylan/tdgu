from __init__ import *

def intr_lens_evalresicntp():

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
    strgexprsbrt = namedatasets + '_%04d.fits' % numbside
    
    if namedatasets == 'lens29075550':
        initbacpbac0en00 = 1.1e-7
        fittmeanbacpbac0en00 = 1.1e-7
        fittstdvbacpbac0en00 = fittmeanbacpbac0en00 * 1e-3
        fittscalbacpbac0en00 = 'gaus'
    else:
        initbacpbac0en00 = None
        fittmeanbacpbac0en00 = None
        fittstdvbacpbac0en00 = None
        fittscalbacpbac0en00 = None

    pcat.main.init( \
                   exprtype='hubb', \
                   makeplotinit=False, \
                   intrevalresicntp=True, \
                   strgexpo=strgexpo, \
                   initbacpbac0en00=initbacpbac0en00, \
                   fittmeanbacpbac0en00=fittmeanbacpbac0en00, \
                   fittstdvbacpbac0en00=fittstdvbacpbac0en00, \
                   fittscalbacpbac0en00=fittscalbacpbac0en00, \
                   inittype='reco', \
                   namerecostat='pcat_lens_inpt', \
                   fittmaxmnumbelem=array([0]), \
                   maxmgangdata=maxmgangdata, \
                   strgexprsbrt=strgexprsbrt, \
                  )
    

def intr_lens_evalcntpresi():
   
    pcat.main.init( \
                   exprtype='hubb', \
                   intrevalcntpresi=True, \
                   numbelempop0=1, \
                   maxmnumbelempop0=1, \
                   numbelempop1=1, \
                   maxmnumbelempop1=1, \
                   numbelempop2=1, \
                   maxmnumbelempop2=1, \
                  )
    

def intr_lens_evalcntpmodl():
   
    pcat.main.init( \
                   exprtype='hubb', \
                   intrevalcntpmodl=True, \
                   numbelempop0=1, \
                   maxmnumbelempop0=1, \
                   numbelempop1=1, \
                   maxmnumbelempop1=1, \
                   numbelempop2=1, \
                   maxmnumbelempop2=1, \
                  )
    

def pcat_lens_mock_truesgnl(strgcnfgextnexec=None):
   
    dictargs = {}
    
    dictargs['exprtype'] = 'hubb'
    dictargs['truenumbelempop0'] = 25
    dictargs['fittmaxmnumbelempop0'] = 25
    dictargs['elemtype'] = ['lens']
   
    listnamecnfgextn = ['nomi', 'truenone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['truenone']['truenumbelempop0'] = 0
    
    lablxaxi = ''
    scalxaxi = 'self'
    listtickxaxi = ['Nominal signal', 'No signal']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='truesgnl', \
                                  lablxaxivari=lablxaxi, \
                                  scalxaxivari=scalxaxi, \
                                  tickxaxivari=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_trueminmdefs(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    
    #dictargs['verbtype'] = 3
    #dictargs['probtran'] = 0.
    #dictargs['probspmr'] = 0.
    #dictargs['probspmr'] = 0.
    #dictargs['makeplot'] = False
    dictargs['exprtype'] = 'hubb'
    
    #dictargs['sqzeexpl'] = True
    #dictargs['optitype'] = 'none'
    
    #dictargs['inittype'] = 'refr'
    dictargs['numbswep'] = 20000
    dictargs['numbswepplot'] = 5000
    dictargs['numbburn'] = 1000
    dictargs['numbsamp'] = 1000

    dictargs['truenumbelempop0'] = 4
    dictargs['fittmaxmnumbelempop0'] = 8
    dictargs['elemtype'] = ['lens']
    #dictargs['elemtype'] = []
   
    listnamecnfgextn = ['truevlow', 'trueloww', 'nomi', 'truehigh', 'truevhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['truevlow']['trueminmdefs'] = 3e-4 / anglfact
    dictargsvari['trueloww']['trueminmdefs'] = 1e-3 / anglfact
    
    #dictargsvari['nomi']['fittminmdefs'] = 3e-20 / anglfact
    #dictargsvari['nomi']['trueminmdefs'] = 3e-20 / anglfact
    
    dictargsvari['nomi']['trueminmdefs'] = 3e-2 / anglfact
    dictargsvari['truehigh']['trueminmdefs'] = 1e-2 / anglfact
    dictargsvari['truevhig']['trueminmdefs'] = 3e-2 / anglfact

    lablxaxi = r'$\alpha_{min}$ [^{\prime\prime}]'
    scalxaxi = 'logt'
    listtickxaxi = [tdpy.util.mexp(anglfact * dictargsvari[namecnfgextn]['trueminmdefs']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  boolexecpara=False, \

                                  namexaxivari='trueminmdefs', \
                                  lablxaxivari=lablxaxi, \
                                  scalxaxivari=scalxaxi, \
                                  tickxaxivari=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_pars(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    
    dictargs = {}
    
    dictargs['exprtype'] = 'hubb'
    dictargs['truenumbelempop0'] = 25
    dictargs['elemtype'] = ['lens']
   
    listnamecnfgextn = ['none', 'loww', 'nomi', 'high', 'vhig']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['none']['priofactdoff'] = 0.
    dictargsvari['loww']['priofactdoff'] = 0.25
    dictargsvari['nomi']['boolsampprio'] = True
    dictargsvari['nomi']['priofactdoff'] = 0.5
    dictargsvari['high']['priofactdoff'] = 0.75
    dictargsvari['vhig']['priofactdoff'] = 1.
    
    #dictargs['numbswep'] = 1000
    #dictargs['numbsamp'] = 10
    #dictargs['verbtype'] = 2
    
    lablxaxi = r'$\alpha_{p}$'
    scalxaxi = 'self'
    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['priofactdoff']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                 
                                  namexaxivari='priofactdoff', \
                                  lablxaxivari=lablxaxi, \
                                  scalxaxivari=scalxaxi, \
                                  tickxaxivari=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_fittnumb(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    
    dictargs['exprtype'] = 'hubb'
    dictargs['truenumbelempop0'] = 25
    dictargs['elemtype'] = ['lens']
   
    listnamecnfgextn = ['fittnone', 'fittsing', 'fittdoub', 'fittmany', 'nomi']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittnone']['fittminmnumbelempop0'] = 0
    dictargsvari['fittnone']['fittmaxmnumbelempop0'] = 0

    dictargsvari['fittsing']['fittminmnumbelempop0'] = 1
    dictargsvari['fittsing']['fittmaxmnumbelempop0'] = 1

    dictargsvari['fittdoub']['fittminmnumbelempop0'] = 2
    dictargsvari['fittdoub']['fittmaxmnumbelempop0'] = 2

    dictargsvari['fittmany']['fittminmnumbelempop0'] = 10
    dictargsvari['fittmany']['fittmaxmnumbelempop0'] = 10

    listtickxaxi = ['0', '1', '2', '10', 'Tran.']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='fittnumbelempop0', \
                                  listtickxaxi=listtickxaxi, \
                        
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_trueback(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['truenumbelempop0'] = 25
    dictargs['elemtype'] = ['lens']
    
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 100
    
    anglfact = 3600. * 180. / pi
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['loww']['truebacpbac0en00'] = 1e-7
    dictargsvari['nomi']['truebacpbac0en00'] = 2e-7
    dictargsvari['high']['truebacpbac0en00'] = 4e-7

    listtickxaxi = [tdpy.util.mexp(dictargsvari[namecnfgextn]['truebacpbac0en00']) for namecnfgextn in listnamecnfgextn] 
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='truebacp', \
                                  tickxaxivari=listtickxaxi, \

                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_lens_mock_sour(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    
    dictargs['elemtype'] = ['lens', 'lghtgausbgrd']
    dictargs['numbelempop0'] = 20
    dictargs['maxmnumbelempop0'] = 100
    dictargs['numbelempop1'] = 10
    dictargs['maxmnumbelempop1'] = 100
    dictargs['spatdisttype'] = ['unif', 'dsrcexpo']
    
    numbelem = int(25. * 10.**0.9)

    listnamecnfgextn = ['nomi', 'datanone', 'subhsing', 'truelowr', 'parsnone', 'truevlow', 's2nrhigh', 's2nrvhig', 'amplhigh', 'bgrdunif']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['datanone']['killexpo'] = True
    
    dictargsvari['bgrdunif']['spatdisttype'] = ['unif', 'unif']
    
    dictargsvari['subhsing']['fittminmnumbelempop0'] = 1
    dictargsvari['subhsing']['fittmaxmnumbelempop0'] = 1

    dictargsvari['truevlow']['truenumbelempop0'] = 30
    
    dictargsvari['truelowr']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['truevlow']['truenumbelempop0'] = int(25. * 10.**0.9)
    dictargsvari['truevlow']['trueminmdefs'] = 3e-4 / anglfact
    dictargsvari['truevlow']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['parsnone']['priofactdoff'] = 0.
    
    dictargsvari['s2nrhigh']['strgexpo'] = 1e4 / 1.63050e-19
    
    dictargsvari['s2nrvhig']['strgexpo'] = 1e5 / 1.63050e-19

    dictargsvari['amplhigh']['minmdefs'] = 1e-1 / anglfact
        
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def test_lens_mock_sele():
    
    numbitermacr = 30
    numbiterelem = 10
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['mockonly'] = True
    dictargs['exprtype'] = 'hubb'
    dictargs['minmdefs'] = 5e-4 / anglfact
    dictargs['numbelempop0'] = 1000
    dictargs['variasca'] = False
    dictargs['variacut'] = False
    dictargs['allwfixdtrue'] = False
    dictargs['maxmnumbelempop0'] = 1000
    
    listnamesele = ['pars', 'nrel']
    numbsele = len(listnamesele)
    listnamefeatsele = ['defs', 'mcut', 'rele']
    numbfeatsele = len(listnamefeatsele)

    dictargsvari = {}
    
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
                                      listnamecnfgextn, \
                                      seedelemtype='rand', \
                                      liststrgvarboutp=liststrgvarboutp, \
                                    
                                      strgcnfgextnexec=strgcnfgextnexec, \
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
    

def test_lens_mock_tmpr():
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['burntmpr'] = True
    dictargs['maxmnumbelempop0'] = 0
    dictargs['numbelempop0'] = 0

    listnamecnfgextn = ['nomi', 'refr', 'pert']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_many(strgcnfgextnexec=None):
    
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['burntmpr'] = True
    for k in range(5):
        dictargs['maxmnumbelempop0reg%d' % k] = 0
        dictargs['numbelempop0reg%d' % k] = 0
    
    dictargs['inittype'] = 'pert'
    dictargs['elemtype'] = ['lens']
    dictargs['numbregi'] = 3
    dictargs['backtype'] = [1., 1., 1.]
    
    dictargs['numbswep'] = 10000
    dictargs['numbsamp'] = 100
    
    listnamecnfgextn = ['nomi', 'regising', 'regimany']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['regising']['numbregi'] = 1 
    dictargsvari['regising']['backtype'] = [1.]
    
    dictargsvari['regimany']['numbregi'] = 5
    dictargsvari['regimany']['backtype'] = [1.] * 5
    
    for k in indxregi:
        dictglob = pcat.main.initarry( \
                                      dictargsvari, \
                                      dictargs, \
                                      listnamecnfgextn, \
                                      
                                      strgcnfgextnexec=strgcnfgextnexec, \
                                     )
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_spmr(strgcnfgextnexec=None):
  
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['inittype'] = 'refr'
    dictargs['numbelempop0'] = 1
    dictargs['truelgalpop00000'] = 1. / anglfact
    dictargs['truebgalpop00000'] = 0.5 / anglfact
    dictargs['truedefspop00000'] = 1e-2 / anglfact
    dictargs['probtran'] = 1.
    dictargs['elemtype'] = ['lens']
    dictargs['probspmr'] = 1.
    dictargs['indxenerincl'] = array([0])
    
    listnamecnfgextn = ['nomi', 'tranboth', 'parshigh', 'masshigh', 'massloww', 'trannone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['tranboth']['inittype'] = 'pert'
    dictargsvari['tranboth']['probtran'] = 0.4
    dictargsvari['tranboth']['probspmr'] = 0.3
    
    dictargsvari['parshigh']['priofactdoff'] = 1.
    
    dictargsvari['masshigh']['truedefspop00000'] = 3e-2 / anglfact
    
    dictargsvari['massloww']['truedefspop00000'] = 3e-3 / anglfact
    
    dictargsvari['trannone']['probtran'] = 0.
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def writ_data():

    # RA/DEC lists
    liststrgrade = []
    listrade = [[], []]
    
    pathbase = os.environ["TDGU_DATA_PATH"] + '/pcat_lens_inpt/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    
    # read SLACS tables
    print('Reading SLACS tables...')
    pathslacpara = pathbase + 'data/slacpara.fits'
    pathslacfull = pathbase + 'data/slacfull.fits'
    hdun = pf.open(pathslacfull)
    numbhead = len(hdun)
    print('%s extensions found.' % numbhead)
    for k in range(numbhead):
        print('Extension %d' % k)
        head = hdun[k].header
        data = hdun[k].data
        
        if data == None:
            print('Data is None, skipping...')
            continue
        else:
            pass
    
        arry = array(stack((head.keys(), head.values()), 1))
        listtype = []
    
        for n in range(arry.shape[0]):
            if arry[n, 0].startswith('TTYPE'):
                listtype.append(arry[n, 1])
        
        if len(listtype) != len(data[0]):
            raise Exception('Number of types does not match the number of fields.')
        
    # find the RA/DEC of relevant SLACS systems 
    indxgold = where((data['Mph'] == 'E') & (data['Mul'] == 'S') & (data['Lens'] == 'A'))[0]
    numbslac = indxgold.size
    
    path = pathdata + 'slacdownlist.txt'
    fileobjt = open(path, 'w')
    for k in indxgold:
        
        # construct the delimited RA/DEC string
        strgrade = '%s %s %s %s %s %s' % (data['SDSS'][k][:2], data['SDSS'][k][2:4], data['SDSS'][k][4:9], data['SDSS'][k][9:12], data['SDSS'][k][12:14], data['SDSS'][k][14:])
        
        ## fill the RA/DEC lists
        liststrgrade.append(strgrade)
        listrade[0].append(data['_RA'][k])
        listrade[1].append(data['_DE'][k])
        
        ## write the RA/DEC list of relevant SLACS systems to disc
        strgline = strgrade + ' \n'
        fileobjt.write(strgline)
    
    for k in range(len(indxgold)):
        print('%20s %20s %20g %20g' % (data['SDSS'][indxgold[k]], data['Name'][indxgold][k], data['_RA'][indxgold][k], data['_DE'][indxgold][k]))
    
    # cutout properties
    numbside = 400
    numbsidehalf = numbside / 2
    
    # data path
    pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    
    ## RA/DEC string of the reference star
    #strgradestar = '00 29 12.65 -00 53 59.7'
    strgradestar = '00 29 06.79 -00 54 07.5'
    liststrgrade.append(strgradestar)
    coorstar = ap.coordinates.SkyCoord(strgradestar, unit=(ap.units.hourangle, ap.units.deg))
    listrade[0].append(coorstar.ra.degree)
    listrade[1].append(coorstar.dec.degree)
    
    numbrade = len(listrade[0])
    print('%d coordinates found.' % numbrade)
    
    # list of files to be read
    listnamefile = ['hst_10886_02_acs_wfc_f814w_drz.fits']
    numbfile = len(listnamefile)
    print('%d files found.' % numbfile)
    
    for k, namefile in enumerate(listnamefile):
        
        print('File number %d' % k)
            
        # read the data fields
        pathfile = pathdata + namefile
        listdata = tdpy.util.read_fits(pathfile, verbtype=0)
        
        # read the WCS header
        listhdun = ap.io.fits.open(pathfile)
        wcso = ap.wcs.WCS(listhdun[2].header)
        
        # RA/DEC string
        strgrade = liststrgrade[k]
    
        # iterate over the RA/DEC list    
        for n in range(numbrade):
            
            # RA/DEC
            strgrade = liststrgrade[n]
            indxyaxi, indxxaxi = wcso.wcs_world2pix(listrade[0][n], listrade[1][n], 0)
            # check if the coordinate is inside the image
            if not isfinite(indxyaxi) or not isfinite(indxxaxi) or indxxaxi - numbsidehalf < 0 or indxyaxi - numbsidehalf < 0 or \
                                                                        indxxaxi + numbsidehalf > listdata[1].shape[1] or indxyaxi + numbsidehalf > listdata[1].shape[0]:
                continue
                #raise Exception('')
    
            path = pathdatapcat + 'lens%s%s%s%s_%04d.fits' % (liststrgrade[n][3:5], liststrgrade[n][6:8], liststrgrade[n][16:18], liststrgrade[n][19:21], numbside)
           
            indxxaxi = int(indxxaxi)
            indxyaxi = int(indxyaxi)
           
            # cut out the image
            rate = listdata[1][indxxaxi-numbsidehalf:indxxaxi+numbsidehalf, indxyaxi-numbsidehalf:indxyaxi+numbsidehalf] # s^-1
    
            # gather different bands
            rate = rate[None, :, :, None]
            
            # find the number of photons per area per time per A per solid angle
            effa = 1. / listdata[4]['PHOTFLAM'][0] # erg^-1 cm^2 A
            timeobsv = listdata[4]['EXPTIME'][0] # s
            apix = (0.05 * pi / 3600. / 180.)**2 # sr^2
            expo = effa * timeobsv # erg^-1 cm^2 s A 
            sbrt = rate / effa / apix
            cntp = sbrt * expo * apix
            
            print('Writing to %s...' % path)
            pf.writeto(path, sbrt, clobber=True)
            
        

def pcat_lens_inpt(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21
    
    # half-size of the image in pixels
    maxmgangdata = 100 * 0.5 * sizepixl
    maxmgangdatalarg = 400 * 0.5 * sizepixl

    # name of the data file
    strgexprsbrt = namedatasets + '_0100.fits'
    strgexprsbrtlarg = namedatasets + '_0400.fits'
    
    if namedatasets == 'lens29075550':
        initlgalsour = -0.1 / anglfact
        initbgalsour = 0.1 / anglfact
        initbacpbac0en00 = 1e-7
        fittmeanbacpbac0en00 = 1.115e-7
        fittstdvbacpbac0en00 = fittmeanbacpbac0en00 * 1e-3
        fittscalbacpbac0en00 = 'gaus'
    else:
        initlgalsour = None
        initbgalsour = None
        initbacpbac0en00 = None
        fittmeanbacpbac0en00 = None
        fittstdvbacpbac0en00 = None
        fittscalbacpbac0en00 = None

    listmask = [
                ['sqre', -0.3, 0.1, -0.1, 0.2] , \
                ['circ', -9, 8, 1] , \
               ]
    for k, mask in enumerate(listmask):
        for n, valu in enumerate(mask):
            if not isinstance(valu, str):
                listmask[k][n] = valu / anglfact
    
    dictargs = {}
    dictargs['elemtype'] = ['lens']
    dictargs['exprtype'] = 'hubb'
    dictargs['strgexpo'] = strgexpo
    dictargs['indxenerincl'] = array([0])
    dictargs['savestat'] = True
    
    # temp
    dictargs['sqzeprop'] = True
    dictargs['numbswep'] = 10000
    dictargs['numbswepplot'] = 1000
    
    dictargs['serstype'] = 'intp'
    #dictargs['verbtype'] = 2
    dictargs['numbsamp'] = 100
    dictargs['inittype'] = 'reco'
    
    #dictargs['inittype'] = 'rand'
    
    listnamecnfgextn = ['largrofi', 'largrofimask', 'nomi', 'mask', 'sour', 'dsrcexpo', 'sourmask', 'hostmult']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargs['initlgalhostisf0'] = -0.1 / anglfact
    dictargs['initbgalhostisf0'] = 0.
    dictargs['initlgalhostisf1'] = -0.1 / anglfact
    dictargs['initbgalhostisf1'] = 0.
    dictargs['initlgalsour'] = -0.2 / anglfact
    dictargs['initbgalsour'] = 0.2 / anglfact
    
    dictargsvari['largrofi']['maxmnumbelempop0'] = 0
    dictargsvari['largrofi']['maxmnumbelempop1'] = 0
    dictargsvari['largrofi']['maxmnumbelempop2'] = 0
    dictargsvari['largrofi']['strgexprsbrt'] = strgexprsbrtlarg
    dictargsvari['largrofi']['maxmgangdata'] = maxmgangdatalarg
    #dictargsvari['largrofi']['thindata'] = True
    
    dictargsvari['largrofimask']['maxmnumbelempop0'] = 0
    dictargsvari['largrofimask']['maxmnumbelempop1'] = 0
    dictargsvari['largrofimask']['maxmnumbelempop2'] = 0
    dictargsvari['largrofimask']['listmask'] = listmask
    dictargsvari['largrofimask']['strgexprsbrt'] = strgexprsbrtlarg
    dictargsvari['largrofimask']['maxmgangdata'] = maxmgangdatalarg
    
    dictargsvari['dsrcexpo']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['dsrcexpo']['maxmnumbelempop0'] = 2
    dictargsvari['dsrcexpo']['elemtype'] = ['lghtgausbgrd']
    #dictargsvari['dsrcexpo']['spatdisttype'] = ['unif', 'dsrcexpo']
    dictargsvari['dsrcexpo']['spatdisttype'] = ['dsrcexpo']
    dictargsvari['dsrcexpo']['dsrcdisttype'] = ['expo']
    dictargsvari['dsrcexpo']['dsrcdistsexppop0'] = 1e-2 / anglfact
    #dictargsvari['dsrcexpo']['numbburn'] = 0
    #dictargsvari['dsrcexpo']['forcsavestat'] = True
    
    dictargsvari['sour']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['sour']['maxmnumbelempop0'] = 1
    dictargsvari['sour']['elemtype'] = ['lghtgausbgrd']
    #dictargsvari['sour']['forcsavestat'] = True
    
    dictargsvari['hostmult']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['hostmult']['numbsersfgrd'] = array([2])
    #dictargsvari['hostmult']['proppsfp'] = False
    dictargsvari['hostmult']['maxmnumbelempop0'] = 0
    dictargsvari['hostmult']['maxmnumbelempop1'] = 0
    dictargsvari['hostmult']['elemtype'] = ['lens', 'lghtgausbgrd']
    #dictargsvari['hostmult']['forcsavestat'] = True
    #dictargsvari['hostmult']['listmask'] = listmask
    
    dictargsvari['sourmask']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['sourmask']['maxmnumbelempop0'] = 0
    dictargsvari['sourmask']['maxmnumbelempop1'] = 0
    dictargsvari['sourmask']['elemtype'] = ['lens', 'lghtgausbgrd']
    dictargsvari['sourmask']['listmask'] = listmask
    
    dictargsvari['nomi']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['nomi']['maxmgangdata'] = maxmgangdata
    dictargsvari['nomi']['maxmnumbelempop0'] = 0
    
    dictargsvari['mask']['listmask'] = listmask
    dictargsvari['mask']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['mask']['maxmgangdata'] = maxmgangdata
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def test_lens_mock_psfn():
   
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['truenumbelempop0'] = 25
    dictargs['elemtype'] = ['lens']
    dictargs['truesigcen00evt0'] = 4e-7
    
    anglfact = 3600. * 180. / pi
    
    listnamecnfgextn = ['loww', 'nomi', 'high']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['psfnfrwr']['truesigcen00evt0'] = 1e-7
    dictargsvari['psfnfrwr']['fittsigcen00evt0'] = 2e-7
    
    dictargsvari['psfnfxwr']['proppsfp'] = False
    dictargsvari['psfnfxwr']['initsigcen00evt0'] = 2e-7
    
    dictargsvari['psfnfrtr']['truesigcen00evt0'] = 4e-7

    dictargsvari['psfnfxtr']['proppsfp'] = False
    dictargsvari['psfnfxtr']['initsigcen00evt0'] = 4e-7

    listtickxaxi = ['Free, Wrong PSF', 'Fixed, Wrong PSF', 'Free, True PSF', 'Fixed, True PSF']
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  namexaxivari='psfnvari', \
                                  tickxaxivari=listtickxaxi, \

                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_papr(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['elemtype'] = ['lens']
    dictargs['seedtype'] = 4
    dictargs['priofactdoff'] = 0.
   
    # temp
    dictargs['limtydathistfeat'] = [0.5, 10.]
    
    dictargs['truemaxmnumbelempop0'] = 25
    dictargs['truenumbelempop0'] = 25
    #dictargs['maxmnumbelempop0'] = 0
    #dictargs['numbelempop0'] = 0
    
    dictargs['truemeanpntspop0'] = 25
    dictargs['initnumbelempop0'] = 25
    dictargs['truedefsdistsloppop0'] = 1.9
    dictargs['truesigcen00evt0'] = 4.21788e-07
    dictargs['truebacp0000en00'] = 2e-07
    dictargs['truelgalsour'] = -2.41109e-07
    dictargs['truebgalsour'] = 1.26909e-07
    dictargs['truefluxsour'] = 1e-18
    dictargs['truesizesour'] = 1.45444e-06
    dictargs['trueellpsour'] = 0.2
    dictargs['trueanglsour'] = 2.4485
    dictargs['truelgalhostisf0'] = 1.10908e-07
    dictargs['truebgalhostisf0'] = 2.26346e-08
    dictargs['truefluxhostisf0'] = 1e-16
    dictargs['truesizehostisf0'] = 4.84814e-06
    dictargs['truebeinhostisf0'] = 7.27221e-06
    dictargs['trueellphostisf0'] = 0.2
    dictargs['trueanglhostisf0'] = 1.21445
    dictargs['trueserihostisf0'] = 4
    dictargs['truesherextrisf0'] = 0.0956653
    dictargs['truesangextrisf0'] = 1.5708

    #dictargs['truelgalpop00000'] = 8.70681e-06
    #dictargs['truebgalpop00000'] = 5.5522e-06
    #dictargs['truedefspop00000'] = 4.46996e-07
    #dictargs['trueascapop00000'] = 8.3953e-08
    #dictargs['trueacutpop00000'] = 7.26722e-07
    #dictargs['truelgalpop00001'] = 1.95366e-06
    #dictargs['truebgalpop00001'] = -6.43887e-06
    #dictargs['truedefspop00001'] = 2.0933e-07
    #dictargs['trueascapop00001'] = 1.98019e-07
    #dictargs['trueacutpop00001'] = 5.11875e-06
    #dictargs['truelgalpop00002'] = 8.48563e-06
    #dictargs['truebgalpop00002'] = 4.20743e-07
    #dictargs['truedefspop00002'] = 5.50444e-08
    #dictargs['trueascapop00002'] = 7.67089e-08
    #dictargs['trueacutpop00002'] = 5.28643e-06
    #dictargs['truelgalpop00003'] = 4.73257e-07
    #dictargs['truebgalpop00003'] = 2.66861e-06
    #dictargs['truedefspop00003'] = 8.56312e-08
    #dictargs['trueascapop00003'] = 3.15034e-07
    #dictargs['trueacutpop00003'] = 3.84845e-06
    #dictargs['truelgalpop00004'] = 2.40305e-06
    #dictargs['truebgalpop00004'] = 5.18566e-06
    #dictargs['truedefspop00004'] = 6.03287e-08
    #dictargs['trueascapop00004'] = 1.82084e-07
    #dictargs['trueacutpop00004'] = 4.8727e-06
    #dictargs['truelgalpop00005'] = 3.61995e-06
    #dictargs['truebgalpop00005'] = -4.77678e-06
    #dictargs['truedefspop00005'] = 1.18797e-07
    #dictargs['trueascapop00005'] = 3.02975e-07
    #dictargs['trueacutpop00005'] = 8.68302e-06
    #dictargs['truelgalpop00006'] = -2.65962e-06
    #dictargs['truebgalpop00006'] = 2.66758e-06
    #dictargs['truedefspop00006'] = 6.1361e-08
    #dictargs['trueascapop00006'] = 2.41337e-07
    #dictargs['trueacutpop00006'] = 1.76904e-06
    #dictargs['truelgalpop00007'] = 8.11351e-06
    #dictargs['truebgalpop00007'] = -1.32214e-06
    #dictargs['truedefspop00007'] = 3.43939e-07
    #dictargs['trueascapop00007'] = 2.02059e-07
    #dictargs['trueacutpop00007'] = 8.7719e-06
    #dictargs['truelgalpop00008'] = -1.84568e-06
    #dictargs['truebgalpop00008'] = -3.27396e-06
    #dictargs['truedefspop00008'] = 1.24152e-07
    #dictargs['trueascapop00008'] = 4.09883e-07
    #dictargs['trueacutpop00008'] = 8.34863e-06
    #dictargs['truelgalpop00009'] = 1.85564e-06
    #dictargs['truebgalpop00009'] = -8.05447e-06
    #dictargs['truedefspop00009'] = 1.32745e-07
    #dictargs['trueascapop00009'] = 1.18999e-07
    #dictargs['trueacutpop00009'] = 7.10343e-06
    #dictargs['truelgalpop00010'] = 7.65329e-06
    #dictargs['truebgalpop00010'] = 2.85729e-07
    #dictargs['truedefspop00010'] = 1.35078e-07
    #dictargs['trueascapop00010'] = 3.15458e-08
    #dictargs['trueacutpop00010'] = 5.23671e-06
    #dictargs['truelgalpop00011'] = -7.19101e-06
    #dictargs['truebgalpop00011'] = 2.22167e-06
    #dictargs['truedefspop00011'] = 8.00093e-08
    #dictargs['trueascapop00011'] = 3.7222e-07
    #dictargs['trueacutpop00011'] =  4.706e-07
    #dictargs['truelgalpop00012'] = -7.56662e-06
    #dictargs['truebgalpop00012'] = 3.56868e-06
    #dictargs['truedefspop00012'] = 1.07991e-07
    #dictargs['trueascapop00012'] = 2.7714e-07
    #dictargs['trueacutpop00012'] = 8.18081e-06
    #dictargs['truelgalpop00013'] = -2.37798e-07
    #dictargs['truebgalpop00013'] = 6.01449e-06
    #dictargs['truedefspop00013'] = 1.06915e-07
    #dictargs['trueascapop00013'] = 4.49287e-07
    #dictargs['trueacutpop00013'] = 6.46671e-06
    #dictargs['truelgalpop00014'] = -6.81208e-06
    #dictargs['truebgalpop00014'] = -2.62666e-06
    #dictargs['truedefspop00014'] = 4.45121e-07
    #dictargs['trueascapop00014'] = 1.69823e-07
    #dictargs['trueacutpop00014'] = 1.83285e-06
    #dictargs['truelgalpop00015'] = -5.30943e-07
    #dictargs['truebgalpop00015'] = -2.07925e-06
    #dictargs['truedefspop00015'] = 1.41112e-07
    #dictargs['trueascapop00015'] = 2.1175e-07
    #dictargs['trueacutpop00015'] = 2.52997e-06
    #dictargs['truelgalpop00016'] = -1.69739e-06
    #dictargs['truebgalpop00016'] = -1.57014e-06
    #dictargs['truedefspop00016'] = 6.30512e-07
    #dictargs['trueascapop00016'] = 4.74931e-07
    #dictargs['trueacutpop00016'] = 6.04629e-06
    #dictargs['truelgalpop00017'] = -8.08312e-06
    #dictargs['truebgalpop00017'] = 4.51844e-06
    #dictargs['truedefspop00017'] = 1.70373e-07
    #dictargs['trueascapop00017'] = 4.00467e-07
    #dictargs['trueacutpop00017'] = 3.36898e-06
    #dictargs['truelgalpop00018'] = -8.55444e-06
    #dictargs['truebgalpop00018'] = 2.16851e-06
    #dictargs['truedefspop00018'] = 5.61476e-08
    #dictargs['trueascapop00018'] = 3.6823e-07
    #dictargs['trueacutpop00018'] = 7.70295e-06
    #dictargs['truelgalpop00019'] = -1.77196e-06
    #dictargs['truebgalpop00019'] = 8.60636e-06
    #dictargs['truedefspop00019'] = 5.99085e-08
    #dictargs['trueascapop00019'] = 4.56979e-07
    #dictargs['trueacutpop00019'] = 4.51352e-06
    #dictargs['truelgalpop00020'] = 5.10012e-06
    #dictargs['truebgalpop00020'] = 4.90283e-06
    #dictargs['truedefspop00020'] = 3.59422e-06
    #dictargs['trueascapop00020'] = 4.71858e-07
    #dictargs['trueacutpop00020'] = 3.88735e-07
    #dictargs['truelgalpop00021'] = -9.38655e-06
    #dictargs['truebgalpop00021'] = 8.39509e-07
    #dictargs['truedefspop00021'] =  5.037e-08
    #dictargs['trueascapop00021'] = 3.18758e-07
    #dictargs['trueacutpop00021'] = 5.18656e-06
    #dictargs['truelgalpop00022'] = 3.22999e-06
    #dictargs['truebgalpop00022'] = -1.13755e-06
    #dictargs['truedefspop00022'] = 6.89786e-08
    #dictargs['trueascapop00022'] = 4.60412e-07
    #dictargs['trueacutpop00022'] = 7.0615e-06
    #dictargs['truelgalpop00023'] = -9.57364e-06
    #dictargs['truebgalpop00023'] = -7.77006e-06
    #dictargs['truedefspop00023'] = 1.46526e-07
    #dictargs['trueascapop00023'] = 1.39644e-07
    #dictargs['trueacutpop00023'] = 3.95241e-06
    #dictargs['truelgalpop00024'] = 5.45961e-06
    #dictargs['truebgalpop00024'] = -2.82849e-06
    #dictargs['truedefspop00024'] = 8.48926e-07
    #dictargs['trueascapop00024'] = 3.49285e-07
    #dictargs['trueacutpop00024'] = 5.35163e-06
    
    anglfact = 3600. * 180. / pi
    
    numbelem = int(25. * 10.**0.9)

    listnamecnfgextn = ['fittlhig', 'fitthigh', 'fittvhig', 'truenone']
    #listnamecnfgextn = ['nomi', 'fittlhig', 'fitthigh', 'fittvhig', 'truenone', 'truenoneparsnomi', 'parsnomi', \
    #                                                                            'subhsing', 'trueloww', 's2nrhigh', 's2nrvhig', 'boolsampprio']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['fittlhig']['fittminmdefs'] = 0.005 / anglfact
    
    dictargsvari['fitthigh']['fittminmdefs'] = 0.01 / anglfact

    dictargsvari['fittvhig']['fittminmdefs'] = 0.02 / anglfact
    
    dictargsvari['truenone']['fittminmdefs'] = 0.01 / anglfact
    dictargsvari['truenone']['truenumbelempop0'] = 0
    
    #dictargsvari['parsnomi']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['parsnomi']['priofactdoff'] = 1.
    
    #dictargsvari['subhsing']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['subhsing']['probtran'] = 0.
    #dictargsvari['subhsing']['initnumbelempop0'] = 1
    #dictargsvari['subhsing']['fittminmnumbelempop0'] = 0
    #dictargsvari['subhsing']['fittmaxmnumbelempop0'] = 1
    #
    #dictargsvari['truenoneparsnomi']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['truenoneparsnomi']['truenumbelempop0'] = 0
    #dictargsvari['truenoneparsnomi']['priofactdoff'] = 1.
    #
    #dictargsvari['trueloww']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['trueloww']['truenumbelempop0'] = int(25. * 10.**0.9)
    #dictargsvari['trueloww']['trueminmdefs'] = 3e-4 / anglfact
    #
    #dictargsvari['s2nrhigh']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['s2nrhigh']['strgexpo'] = 1e4 / 1.63050e-19
    #
    #dictargsvari['s2nrvhig']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['s2nrvhig']['strgexpo'] = 1e5 / 1.63050e-19

    #dictargsvari['prio']['fittminmdefs'] = 0.01 / anglfact
    #dictargsvari['prio']['boolsampprio'] = True
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )

    

globals().get(sys.argv[1])(*sys.argv[2:])
