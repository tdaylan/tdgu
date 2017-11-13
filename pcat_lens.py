from __init__ import *

def pcat_lens_mock_grid():

    listnameoutpvarb = ['maxmllik', 'medilliktotl', 'stdvlliktotl', 'levi', 'info']

    gdat = pcat.main.init(exprtype='hubb', defa=True, verbtype=0)

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
    dictvarb['exprtype'] = 'hubb'
    dictvarb['numbswep'] = 100
    dictvarb['makeplot'] = False
    
    cntrcnfg = 0
    for k in range(numbiter):
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
    

def pcat_lens_intrevalcntpresi():
   
    pcat.main.init( \
                   exprtype='hubb', \
                   makeplotinit=False, \
                   intrevalcntpresi=True, \
                   numbelempop0reg0=1, \
                   maxmnumbelempop0reg0=1, \
                   numbelempop1reg0=1, \
                   maxmnumbelempop1reg0=1, \
                   numbelempop2reg0=1, \
                   maxmnumbelempop2reg0=1, \
                  )
    

def pcat_lens_intrevalcntpmodl():
   
    pcat.main.init( \
                   exprtype='hubb', \
                   makeplotinit=False, \
                   intrevalcntpmodl=True, \
                   numbelempop0reg0=1, \
                   maxmnumbelempop0reg0=1, \
                   numbelempop1reg0=1, \
                   maxmnumbelempop1reg0=1, \
                   numbelempop2reg0=1, \
                   maxmnumbelempop2reg0=1, \
                  )
    

def pcat_lens_mock_sing():
   
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       exprtype='hubb', \
                       numbelempop0reg0=1, \
                       minmdefs=1e-2/3600./180.*pi, \
                       maxmdefs=1e-1/3600./180.*pi, \
                      )
    

def pcat_lens_mock_next(strgcnfgextnexec=None):
   
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    #dictargs['numbswep'] = 1000
    
    dictargs['exprtype'] = 'hubb'
    dictargs['truenumbelempop0reg0'] = 25
    dictargs['elemtype'] = ['lens']
   
    # temp
    #dictargs['makeplotinit'] = True
    #dictargs['shrtfram'] = False
    dictargs['numbswep'] = 50000
    dictargs['numbburn'] = 20000
    dictargs['factthin'] = 300
    dictargs['probspmr'] = 0.
    dictargs['priofactdoff'] = 0.
 
    #dictargsvari['numbelempop0reg0']     = [None,        0,  0,  25,   int(25. * 0.1**0.9), int(25. * 10.**0.9)]
    #dictargsvari['trueminmdefs']     = [None,        None,        None,        3e-3/anglfact, 3e-2/anglfact,                      3e-4/anglfact]
    #dictargsvari['fittminmdefs']     = [None,        None,        None,        3e-4/anglfact, 3e-4/anglfact,                      3e-4/anglfact]
    #dictargsvari['priofactdoff']     = [0.,          0.,          1.,          1.,            1.,                                 1.]
    #dictargsvari['scalmeanpnts'] = ['logt',      'logt',      'logt',      'logt',        'logt',                            'logt']
   
    listnamecnfgextn = ['nomi', 'parsnone', 'parsloww', 'parsmore', 'parshigh', 'zerosgnl']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['zerosgnl']['truenumbelempop0reg0'] = 0
    dictargsvari['zerosgnl']['truemaxmnumbelempop0reg0'] = 0
    dictargsvari['zerosgnl']['fittmaxmnumbelempop0reg0'] = 100
    
    dictargsvari['parsnone']['priofactdoff'] = 0.
    dictargsvari['parsloww']['priofactdoff'] = 0.5
    dictargsvari['parsmore']['priofactdoff'] = 1.5
    dictargsvari['parshigh']['priofactdoff'] = 2.
    
    #dictargsvari['subhsing']['fittminmnumbelempop0reg0'] = 1
    #dictargsvari['subhsing']['fittmaxmnumbelempop0reg0'] = 1

    #dictargsvari['truelowr']['truenumbelempop0reg0'] = int(25. * 10.**0.9)
    #dictargsvari['truelowr']['trueminmdefs'] = 3e-3 / anglfact
    #dictargsvari['truelowr']['fittminmdefs'] = 3e-3 / anglfact
    #

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock_syst(strgcnfgextnexec=None):
   
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['elemtype'] = ['lens']
    dictargs['seedtype'] = 4
    
    dictargs['initnumbelempop0reg0'] = 25
    
    dictargs['truenumbelempop0reg0'] = 25
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

    dictargs['truelgalpop0reg00000'] = 8.70681e-06
    dictargs['truebgalpop0reg00000'] = 5.5522e-06
    dictargs['truedefspop0reg00000'] = 4.46996e-07
    dictargs['trueascapop0reg00000'] = 8.3953e-08
    dictargs['trueacutpop0reg00000'] = 7.26722e-07
    dictargs['truelgalpop0reg00001'] = 1.95366e-06
    dictargs['truebgalpop0reg00001'] = -6.43887e-06
    dictargs['truedefspop0reg00001'] = 2.0933e-07
    dictargs['trueascapop0reg00001'] = 1.98019e-07
    dictargs['trueacutpop0reg00001'] = 5.11875e-06
    dictargs['truelgalpop0reg00002'] = 8.48563e-06
    dictargs['truebgalpop0reg00002'] = 4.20743e-07
    dictargs['truedefspop0reg00002'] = 5.50444e-08
    dictargs['trueascapop0reg00002'] = 7.67089e-08
    dictargs['trueacutpop0reg00002'] = 5.28643e-06
    dictargs['truelgalpop0reg00003'] = 4.73257e-07
    dictargs['truebgalpop0reg00003'] = 2.66861e-06
    dictargs['truedefspop0reg00003'] = 8.56312e-08
    dictargs['trueascapop0reg00003'] = 3.15034e-07
    dictargs['trueacutpop0reg00003'] = 3.84845e-06
    dictargs['truelgalpop0reg00004'] = 2.40305e-06
    dictargs['truebgalpop0reg00004'] = 5.18566e-06
    dictargs['truedefspop0reg00004'] = 6.03287e-08
    dictargs['trueascapop0reg00004'] = 1.82084e-07
    dictargs['trueacutpop0reg00004'] = 4.8727e-06
    dictargs['truelgalpop0reg00005'] = 3.61995e-06
    dictargs['truebgalpop0reg00005'] = -4.77678e-06
    dictargs['truedefspop0reg00005'] = 1.18797e-07
    dictargs['trueascapop0reg00005'] = 3.02975e-07
    dictargs['trueacutpop0reg00005'] = 8.68302e-06
    dictargs['truelgalpop0reg00006'] = -2.65962e-06
    dictargs['truebgalpop0reg00006'] = 2.66758e-06
    dictargs['truedefspop0reg00006'] = 6.1361e-08
    dictargs['trueascapop0reg00006'] = 2.41337e-07
    dictargs['trueacutpop0reg00006'] = 1.76904e-06
    dictargs['truelgalpop0reg00007'] = 8.11351e-06
    dictargs['truebgalpop0reg00007'] = -1.32214e-06
    dictargs['truedefspop0reg00007'] = 3.43939e-07
    dictargs['trueascapop0reg00007'] = 2.02059e-07
    dictargs['trueacutpop0reg00007'] = 8.7719e-06
    dictargs['truelgalpop0reg00008'] = -1.84568e-06
    dictargs['truebgalpop0reg00008'] = -3.27396e-06
    dictargs['truedefspop0reg00008'] = 1.24152e-07
    dictargs['trueascapop0reg00008'] = 4.09883e-07
    dictargs['trueacutpop0reg00008'] = 8.34863e-06
    dictargs['truelgalpop0reg00009'] = 1.85564e-06
    dictargs['truebgalpop0reg00009'] = -8.05447e-06
    dictargs['truedefspop0reg00009'] = 1.32745e-07
    dictargs['trueascapop0reg00009'] = 1.18999e-07
    dictargs['trueacutpop0reg00009'] = 7.10343e-06
    dictargs['truelgalpop0reg00010'] = 7.65329e-06
    dictargs['truebgalpop0reg00010'] = 2.85729e-07
    dictargs['truedefspop0reg00010'] = 1.35078e-07
    dictargs['trueascapop0reg00010'] = 3.15458e-08
    dictargs['trueacutpop0reg00010'] = 5.23671e-06
    dictargs['truelgalpop0reg00011'] = -7.19101e-06
    dictargs['truebgalpop0reg00011'] = 2.22167e-06
    dictargs['truedefspop0reg00011'] = 8.00093e-08
    dictargs['trueascapop0reg00011'] = 3.7222e-07
    dictargs['trueacutpop0reg00011'] =  4.706e-07
    dictargs['truelgalpop0reg00012'] = -7.56662e-06
    dictargs['truebgalpop0reg00012'] = 3.56868e-06
    dictargs['truedefspop0reg00012'] = 1.07991e-07
    dictargs['trueascapop0reg00012'] = 2.7714e-07
    dictargs['trueacutpop0reg00012'] = 8.18081e-06
    dictargs['truelgalpop0reg00013'] = -2.37798e-07
    dictargs['truebgalpop0reg00013'] = 6.01449e-06
    dictargs['truedefspop0reg00013'] = 1.06915e-07
    dictargs['trueascapop0reg00013'] = 4.49287e-07
    dictargs['trueacutpop0reg00013'] = 6.46671e-06
    dictargs['truelgalpop0reg00014'] = -6.81208e-06
    dictargs['truebgalpop0reg00014'] = -2.62666e-06
    dictargs['truedefspop0reg00014'] = 4.45121e-07
    dictargs['trueascapop0reg00014'] = 1.69823e-07
    dictargs['trueacutpop0reg00014'] = 1.83285e-06
    dictargs['truelgalpop0reg00015'] = -5.30943e-07
    dictargs['truebgalpop0reg00015'] = -2.07925e-06
    dictargs['truedefspop0reg00015'] = 1.41112e-07
    dictargs['trueascapop0reg00015'] = 2.1175e-07
    dictargs['trueacutpop0reg00015'] = 2.52997e-06
    dictargs['truelgalpop0reg00016'] = -1.69739e-06
    dictargs['truebgalpop0reg00016'] = -1.57014e-06
    dictargs['truedefspop0reg00016'] = 6.30512e-07
    dictargs['trueascapop0reg00016'] = 4.74931e-07
    dictargs['trueacutpop0reg00016'] = 6.04629e-06
    dictargs['truelgalpop0reg00017'] = -8.08312e-06
    dictargs['truebgalpop0reg00017'] = 4.51844e-06
    dictargs['truedefspop0reg00017'] = 1.70373e-07
    dictargs['trueascapop0reg00017'] = 4.00467e-07
    dictargs['trueacutpop0reg00017'] = 3.36898e-06
    dictargs['truelgalpop0reg00018'] = -8.55444e-06
    dictargs['truebgalpop0reg00018'] = 2.16851e-06
    dictargs['truedefspop0reg00018'] = 5.61476e-08
    dictargs['trueascapop0reg00018'] = 3.6823e-07
    dictargs['trueacutpop0reg00018'] = 7.70295e-06
    dictargs['truelgalpop0reg00019'] = -1.77196e-06
    dictargs['truebgalpop0reg00019'] = 8.60636e-06
    dictargs['truedefspop0reg00019'] = 5.99085e-08
    dictargs['trueascapop0reg00019'] = 4.56979e-07
    dictargs['trueacutpop0reg00019'] = 4.51352e-06
    dictargs['truelgalpop0reg00020'] = 5.10012e-06
    dictargs['truebgalpop0reg00020'] = 4.90283e-06
    dictargs['truedefspop0reg00020'] = 3.59422e-06
    dictargs['trueascapop0reg00020'] = 4.71858e-07
    dictargs['trueacutpop0reg00020'] = 3.88735e-07
    dictargs['truelgalpop0reg00021'] = -9.38655e-06
    dictargs['truebgalpop0reg00021'] = 8.39509e-07
    dictargs['truedefspop0reg00021'] =  5.037e-08
    dictargs['trueascapop0reg00021'] = 3.18758e-07
    dictargs['trueacutpop0reg00021'] = 5.18656e-06
    dictargs['truelgalpop0reg00022'] = 3.22999e-06
    dictargs['truebgalpop0reg00022'] = -1.13755e-06
    dictargs['truedefspop0reg00022'] = 6.89786e-08
    dictargs['trueascapop0reg00022'] = 4.60412e-07
    dictargs['trueacutpop0reg00022'] = 7.0615e-06
    dictargs['truelgalpop0reg00023'] = -9.57364e-06
    dictargs['truebgalpop0reg00023'] = -7.77006e-06
    dictargs['truedefspop0reg00023'] = 1.46526e-07
    dictargs['trueascapop0reg00023'] = 1.39644e-07
    dictargs['trueacutpop0reg00023'] = 3.95241e-06
    dictargs['truelgalpop0reg00024'] = 5.45961e-06
    dictargs['truebgalpop0reg00024'] = -2.82849e-06
    dictargs['truedefspop0reg00024'] = 8.48926e-07
    dictargs['trueascapop0reg00024'] = 3.49285e-07
    dictargs['trueacutpop0reg00024'] = 5.35163e-06
   
    # temp
    #dictargs['makeplotinit'] = True
    #dictargs['shrtfram'] = False
    dictargs['numbswep'] = 100000
    #dictargs['verbtype'] = 2
    
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / pi

    listnamecnfgextn = ['nomi', 'truelowr', 'truelowrparsnone', 'truenone', 'truenoneparsnone', 'subhsing', 'truevlow', 's2nrhigh', 's2nrvhig', 'datanone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['truelowr']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['truelowrparsnone']['fittminmdefs'] = 0.01 / anglfact
    dictargsvari['truelowrparsnone']['priofactdoff'] = 0.

    dictargsvari['truenone']['truenumbelempop0reg0'] = 0
    
    dictargsvari['truenoneparsnone']['truenumbelempop0reg0'] = 0
    dictargsvari['truenoneparsnone']['priofactdoff'] = 0.
    
    dictargsvari['subhsing']['probtran'] = 0.
    dictargsvari['subhsing']['initnumbelempop0reg0'] = 1
    dictargsvari['subhsing']['fittminmnumbelempop0reg0'] = 0
    dictargsvari['subhsing']['fittmaxmnumbelempop0reg0'] = 1
    
    dictargsvari['truevlow']['truenumbelempop0reg0'] = int(25. * 10.**0.9)
    dictargsvari['truevlow']['trueminmdefs'] = 3e-4 / anglfact
    dictargsvari['truevlow']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['s2nrhigh']['strgexpo'] = 1e4 / 1.63050e-19
    
    dictargsvari['s2nrvhig']['strgexpo'] = 1e5 / 1.63050e-19

    dictargsvari['datanone']['killexpo'] = True
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_lens_mock_sour(strgcnfgextnexec=None):
   
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    
    dictargs['numbelempop0reg0'] = 2
    dictargs['maxmnumbelempop0reg0'] = 20
    dictargs['numbelempop1reg0'] = 50
    dictargs['maxmnumbelempop1reg0'] = 100
    
    dictargs['numbelempop2reg0'] = 0
    dictargs['maxmnumbelempop2reg0'] = 0
    
    dictargs['truelgalsourreg0'] = 0.
    dictargs['truebgalsourreg0'] = 0.
    dictargs['truefluxsourreg0'] = 3e-18
    dictargs['truefluxpop0reg00000'] = 1e-18
    dictargs['truefluxpop0reg00001'] = 3e-18
    dictargs['truelgalpop2reg00000'] = -0.2 / anglfact
    dictargs['truebgalpop2reg00000'] = 0.
    dictargs['truefluxpop2reg00000'] = 1e-19
    dictargs['truelgalpop2reg00001'] = -0.1 / anglfact
    dictargs['truebgalpop2reg00001'] = 0.
    dictargs['truefluxpop2reg00001'] = 3e-19
    dictargs['truelgalpop2reg00002'] = 0.1 / anglfact
    dictargs['truebgalpop2reg00002'] = 0.
    dictargs['truefluxpop2reg00002'] = 3e-19
    dictargs['truelgalpop2reg00003'] = 0.2 / anglfact
    dictargs['truebgalpop2reg00003'] = 0.
    dictargs['truefluxpop2reg00003'] = 1e-19
    #dictargs['refrlegdpopl'] = ['PS', 'Subhalo', 'Blob']
   
    # temp
    #dictargs['makeplotinit'] = True
    #dictargs['shrtfram'] = False
    dictargs['numbswep'] = 40000
    dictargs['numbswepplot'] = 4000
    dictargs['numbsamp'] = 100
    dictargs['probspmr'] = 0.
    #dictargs['spatdisttype'] = ['unif', 'unif', 'gangprop']
    #dictargs['makeplot'] = False
    
    #dictargs['inittype'] = 'refr'
    #dictargs['optitype'] = 'none'
    #dictargs['verbtype'] = 2
    
    numbelem = int(25. * 10.**0.9)
    anglfact = 3600. * 180. / pi

    listnamecnfgextn = ['nomi', 'datanone', 'subhsing', 'truelowr', 'pars', 'truevlow', 's2nrhigh', 's2nrvhig', 'amplhigh']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['datanone']['killexpo'] = True
    
    dictargsvari['subhsing']['fittminmnumbelempop0reg0'] = 1
    dictargsvari['subhsing']['fittmaxmnumbelempop0reg0'] = 1

    dictargsvari['truevlow']['truenumbelempop0reg0'] = 30
    
    dictargsvari['truelowr']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['truevlow']['truenumbelempop0reg0'] = int(25. * 10.**0.9)
    dictargsvari['truevlow']['trueminmdefs'] = 3e-4 / anglfact
    dictargsvari['truevlow']['fittminmdefs'] = 0.01 / anglfact
    
    dictargsvari['pars']['priofactdoff'] = 0.
    
    dictargsvari['s2nrhigh']['strgexpo'] = 1e4 / 1.63050e-19
    
    dictargsvari['s2nrvhig']['strgexpo'] = 1e5 / 1.63050e-19

    dictargsvari['amplhigh']['minmdefs'] = 1e-1 / anglfact
        
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

def pcat_lens_mock_perf():
   
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    
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


def pcat_lens_mock_sele():
    
    numbitermacr = 30
    numbiterelem = 10
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['mockonly'] = True
    dictargs['exprtype'] = 'hubb'
    dictargs['minmdefs'] = 5e-4 / anglfact
    dictargs['numbelempop0reg0'] = 1000
    dictargs['variasca'] = False
    dictargs['variacut'] = False
    dictargs['allwfixdtrue'] = False
    dictargs['verbtype'] = 0
    dictargs['makeplot'] = False
    dictargs['maxmnumbelempop0reg0'] = 1000
    
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
                                      seedelemtype='rand', \
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
    

def pcat_lens_mock_tmpr():
    
    anglfact = 3600. * 180. / pi
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['burntmpr'] = True
    dictargs['maxmnumbelempop0reg0'] = 0
    dictargs['numbelempop0reg0'] = 0
    dictargsvari = {}
    dictargsvari['inittype'] = ['refr', 'rand']

    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_macr():
    
    pcat.main.init( \
                   exprtype='hubb', \
                   maxmnumbelempop0reg0=0, \
                   numbelempop0reg0=0, \
                  )


def pcat_lens_mock_test():
   
    dictargs = {}
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


def pcat_lens_mock_many():
    
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['burntmpr'] = True
    dictargs['maxmnumbelempop0reg0'] = 0
    dictargs['numbelempop0reg0'] = 0
    
    #dictargs['makeplotinit'] = False
    #dictargs['shrtfram'] = False
    #dictargs['verbtype'] = 2
    dictargs['inittype'] = 'pert'
    dictargs['numbswep'] = 100000
    dictargs['numbswepplot'] = 10000
    dictargs['numbsamp'] = 100
    dictargs['elemtype'] = ['lens']
    dictargs['numbregi'] = 3
    
    listnamecnfgextn = ['nomi', 'regising', 'regimany']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['regising']['numbregi'] = 1 
    dictargsvari['regimany']['numbregi'] = 5
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                 )


def pcat_lens_mock_spmr(strgcnfgextnexec=None):
  
    anglfact = 3600. * 180. / pi
    
    dictargs = {}
    dictargs['exprtype'] = 'hubb'
    dictargs['inittype'] = 'refr'
    dictargs['numbelempop0reg0'] = 1
    dictargs['truelgalpop0reg00000'] = 1. / anglfact
    dictargs['truebgalpop0reg00000'] = 0.5 / anglfact
    dictargs['truedefspop0reg00000'] = 1e-2 / anglfact
    dictargs['priofactdoff'] = 0.
    dictargs['probtran'] = 1.
    dictargs['elemtype'] = ['lens']
    dictargs['probspmr'] = 1.
    dictargs['indxenerincl'] = array([0])
    
    # temp
    #dictargs['makeplotinit'] = False
    #dictargs['shrtfram'] = True
    dictargs['numbswep'] = 10000
    dictargs['numbburn'] = 0
    dictargs['numbsamp'] = 1000
    dictargs['numbswepplot'] = 2000
    #dictargs['verbtype'] = 2
    
    listnamecnfgextn = ['nomi', 'tranboth', 'parshigh', 'masshigh', 'massloww', 'trannone']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    dictargsvari['tranboth']['inittype'] = 'pert'
    dictargsvari['tranboth']['probtran'] = 0.4
    dictargsvari['tranboth']['probspmr'] = 0.3
    
    dictargsvari['parshigh']['priofactdoff'] = 1.
    
    dictargsvari['masshigh']['truedefspop0reg00000'] = 3e-2 / anglfact
    
    dictargsvari['massloww']['truedefspop0reg00000'] = 3e-3 / anglfact
    
    dictargsvari['trannone']['probtran'] = 0.
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_mock():
   
    pcat.main.init( \
                   exprtype='hubb', \
                   numbswep=100, \
                   #numbsamp=10, \
                   numbburn=0, \
                   factthin=100, \
                   inittype='refr', \
                   elemtype=['lens'], \
                   numbelempop0reg0=1, \
                   maxmnumbelempop0reg0=1, \
                   #numbelempop0reg0=0, \
                   #numbelempop1reg0=1, \
                   #numbelempop2reg0=0, \
                   #maxmnumbelempop0reg0=0, \
                   #maxmnumbelempop1reg0=3, \
                   #maxmnumbelempop2reg0=0, \
                   makeplot=False, \
                   
                   #makeplotintr=True, \
                   #makeplotinit=True, \

                   # temp
                   #probtran=0., \
                   #sqzeprop=True, \
                   #maxmnumbelempop0reg0=1, \
                   #maxmnumbelempop1reg0=1, \
                   #maxmnumbelempop2reg0=1, \
                   #makeplot=False, \
                   #makeplotinit=False, \
                   #shrtfram=True, \
                   #explprop=True, \
                   #verbtype=2, \
                   #mockonly=True, \
                  )


def writ_data():

    # RA/DEC lists
    liststrgrade = []
    listrade = [[], []]
    
    pathbase = os.environ["TDGU_DATA_PATH"] + '/pcat_lens_inpt/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    
    # read SLACS tables
    print 'Reading SLACS tables...'
    pathslacpara = pathbase + 'data/slacpara.fits'
    pathslacfull = pathbase + 'data/slacfull.fits'
    hdun = pf.open(pathslacfull)
    numbhead = len(hdun)
    print '%s extensions found.' % numbhead
    for k in range(numbhead):
        print 'Extension %d' % k
        head = hdun[k].header
        data = hdun[k].data
        
        if data == None:
            print 'Data is None, skipping...'
            continue
        else:
            pass
            #print 'data object has keys'
            #print data.names
    
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
        print '%20s %20s %20g %20g' % (data['SDSS'][indxgold[k]], data['Name'][indxgold][k], data['_RA'][indxgold][k], data['_DE'][indxgold][k])
    
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
    print '%d coordinates found.' % numbrade
    
    # list of files to be read
    listnamefile = ['hst_10886_02_acs_wfc_f814w_drz.fits']
    numbfile = len(listnamefile)
    print '%d files found.' % numbfile
    
    for k, namefile in enumerate(listnamefile):
        
        print 'File number %d' % k
            
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
           
            if False:
                print 'listdata[4]'
                print listdata[4].names
                print 'strgrade'
                print strgrade
                print 'listrade[0][n]'
                print listrade[0][n]
                print 'listrade[1][n]'
                print listrade[1][n]
                print 'indxxaxi'
                print indxxaxi
                print 'indxyaxi'
                print indxyaxi
            
            indxxaxi = int(indxxaxi)
            indxyaxi = int(indxyaxi)
           
    
            if False:
                print 'MDRIZSKY'
                print listdata[4]['MDRIZSKY']
                print 'SKYSUB'
                print listdata[4]['SKYSUB']
                print 'indxxaxi'
                print indxxaxi
                print 'indxyaxi'
                print indxyaxi
                print 'listdata[1]'
                summgene(listdata[1])
                print 'EXPTIME'
                print listdata[4]['EXPTIME'][0]
                print 'PHOTFLAM'
                print listdata[4]['PHOTFLAM'][0]
                print 'CCDGAIN'
                print listdata[4]['CCDGAIN'][0]
                print
            
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
            
            if False:
                print 'expo'
                print expo
                print 'rate'
                summgene(rate)
                print 'mean(cntp[:10, :10])'
                print mean(cntp[:10, :10])
                print 'cntp'
                summgene(cntp)
                print 'sbrt'
                summgene(sbrt)
            if True:
                print 'sbrt'
                summgene(sbrt)
                print 'apix'
                print apix
            
            print 'Writing to %s...' % path
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
        initlgalsourreg0 = -0.1 / anglfact
        initbgalsourreg0 = 0.1 / anglfact
        initbacpbac0ene0 = 1e-7
        fittmeanbacpbac0ene0 = 1.115e-7
        fittstdvbacpbac0ene0 = fittmeanbacpbac0ene0 * 1e-3
        fittscalbacpbac0ene0 = 'gaus'
    else:
        initlgalsourreg0 = None
        initbgalsourreg0 = None
        initbacpbac0ene0 = None
        fittmeanbacpbac0ene0 = None
        fittstdvbacpbac0ene0 = None
        fittscalbacpbac0ene0 = None

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
    dictargs['serstype'] = 'intp'
    dictargs['inittype'] = 'reco'
    #dictargs['namerecostat'] = 'pcat_lens_inpt'
    dictargs['maxmnumbelempop0reg0'] = 5
    dictargs['maxmnumbelempop1reg0'] = 5
    dictargs['maxmnumbelempop2reg0'] = 5
    # temp
    dictargs['sqzeprop'] = True
    dictargs['numbswep'] = 1000
    dictargs['numbsamp'] = 100
    dictargs['makeplotinit'] = False
    dictargs['numbswepplot'] = 1000
    dictargs['shrtfram'] = True
    dictargs['verbtype'] = 2
 
    listnamecnfgextn = ['largrofi', 'largrofimask', 'nomi', 'mask']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    
    #dictargsvari['largrofi']['numbswep'] = 10000
    dictargsvari['largrofi']['numbswepplot'] = 1000
    dictargsvari['largrofi']['maxmnumbelempop0reg0'] = 0
    dictargsvari['largrofi']['maxmnumbelempop1reg0'] = 0
    dictargsvari['largrofi']['maxmnumbelempop2reg0'] = 0
    dictargsvari['largrofi']['strgexprsbrt'] = strgexprsbrtlarg
    dictargsvari['largrofi']['maxmgangdata'] = maxmgangdatalarg
    #dictargsvari['largrofi']['thindata'] = True
    #dictargsvari['largrofi']['initlgalhostreg0'] = 0.
    #dictargsvari['largrofi']['initbgalhostreg0'] = 0.
    #dictargsvari['largrofi']['initlgalsourreg0'] = 0.
    #dictargsvari['largrofi']['initbgalsourreg0'] = 0.
    
    #dictargsvari['largrofimask']['numbswep'] = 10000
    dictargsvari['largrofimask']['numbswepplot'] = 1000
    dictargsvari['largrofimask']['maxmnumbelempop0reg0'] = 0
    dictargsvari['largrofimask']['maxmnumbelempop1reg0'] = 0
    dictargsvari['largrofimask']['maxmnumbelempop2reg0'] = 0
    dictargsvari['largrofimask']['listmask'] = listmask
    dictargsvari['largrofimask']['strgexprsbrt'] = strgexprsbrtlarg
    dictargsvari['largrofimask']['maxmgangdata'] = maxmgangdatalarg
    
    dictargsvari['nomi']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['nomi']['maxmgangdata'] = maxmgangdata
    
    dictargsvari['mask']['listmask'] = listmask
    dictargsvari['mask']['strgexprsbrt'] = strgexprsbrt
    dictargsvari['mask']['maxmgangdata'] = maxmgangdata
    
    #fittmeanbacpbac0ene0=fittmeanbacpbac0ene0, \
    #fittstdvbacpbac0ene0=fittstdvbacpbac0ene0, \
    #fittscalbacpbac0ene0=fittscalbacpbac0ene0, \
   
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def pcat_lens_psfn():
    
    strgexpo = 7.37487548893e21
    anglfact = 3600. * 180. / pi
    maxmgangdata = 50. * 0.05 / anglfact
    numbiter = 1
    
    for k in range(numbiter):
        pcat.main.init( \
                       exprtype='hubb', \
                       numbswep=50000, \
                       factthin=500, \
                       numbswepplot=10000, \
                       shrtfram=True, \
                       #mockonly=True, \
                       #makeplotintr=True, \
                       #burntmpr=True, \
                       #savestat=True, \
                       #inittype='reco', \
                       #makeplotinit=False, \
                       makeplotlpri=False, \
                       strgexpo=strgexpo, \
                       fittmaxmnumbelem=array([0]), \
                       maxmgangdata=maxmgangdata, \
                       strgexprsbrt='lens29065407.fits', \
                      )
   

def pcat_lens_intrevalresicntp():

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
                   exprtype='hubb', \
                   makeplotinit=False, \
                   intrevalresicntp=True, \
                   strgexpo=strgexpo, \
                   initbacpbac0ene0=initbacpbac0ene0, \
                   fittmeanbacpbac0ene0=fittmeanbacpbac0ene0, \
                   fittstdvbacpbac0ene0=fittstdvbacpbac0ene0, \
                   fittscalbacpbac0ene0=fittscalbacpbac0ene0, \
                   inittype='reco', \
                   namerecostat='pcat_lens_inpt', \
                   fittmaxmnumbelem=array([0]), \
                   maxmgangdata=maxmgangdata, \
                   strgexprsbrt=strgexprsbrt, \
                  )
    

def pcat_lens_mockonly():
   
    pcat.main.init( \
                   exprtype='hubb', \
                   maxmnumbelempop0reg0=0, \
                   numbelempop0reg0=0, \
                   maxmnumbelempop1reg0=1000, \
                   numbelempop1reg0=1000, \
                   maxmnumbelempop2reg0=0, \
                   numbelempop2reg0=0, \
                   makeplotinit=True, \
                   makeplotintr=True, \
                   mockonly=True, \
                  )


globals().get(sys.argv[1])(*sys.argv[2:])
