from __init__ import *

def pcat_lens_mock_grid():
  
    gdat = pcat.main.init(verbtype=0, pntstype='lens', exprtype='hubb', defa=True)

    truespecsour = gdat.hubbexpofact
    truespechost = htsr * truespecsour
    back = [truespecsour * btsr / (pi * gdat.truefixp[gdat.trueindxfixpsizesour]**2) / 2.]
    
    numbcnfg = 5
    numbiter = 2
   
    varbinpt = tdpy.util.varb(numbcnfg)
    varbinpt.defn_para('minmflux', 3e-4, 3e-2, scal='logt')
    varbinpt.defn_para('back', 1e-1 * fact / (pi * truesizesour**2), 1e1 * fact / (pi * truesizesour**2), scal='logt')
    varbinpt.defn_para('truespecsour', 1e-1 * fact, 1e1 * fact, scal='logt')
    varbinpt.defn_para('truespechost', 1e-1 * fact, 1e1 * fact, scal='logt')

    listnameoutpvarb = ['fluxdistsloppop0', 'meanpntspop0', 'specassc', 'spechost', 'beinhost']
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
                     
                    if varbinpt.name[p] == 'minmflux':
                        minmflux = varb / anglfact
                    if varbinpt.name[p] == 'back':
                        back = varb
                    if varbinpt.name[p] == 'strgexpo':
                        strgexpo = varb
                    if varbinpt.name[p] == 'truespecsour':
                        truespecsour = varb
                    if varbinpt.name[p] == 'truespechost':
                        truespechost = varb
                
                gdat = pcat.main.init( \
                                      seedstat=seedstat, \
                                      numbswep=100000, \
                                      numbswepplot=4000, \
                                      factthin=200, \
                                      pntstype='lens', \
                                      strgexpo=strgexpo, \
                                      exprtype='hubb', \
                                      back=[back], \
                                      minmflux=minmflux, \
                                      truespecsour=truespecsour, \
                                      truespechost=truespechost, \
                                      truenumbpnts=array([10]), \
                                     )
                
                for n in range(numboutpvarb):
                    if listnameoutpvarb[n] == 'specassc':
                        grid[0, n, m, l] = gdat.fluxfactplot * gdat.medispecassc[0][gdat.indxenerfluxdist[0], 0]
                        grid[1:3, n, m, l] = gdat.fluxfactplot * gdat.errrspecassc[0][:, gdat.indxenerfluxdist[0], 0]
                        grid[3, n, m, l] = gdat.fluxfactplot * gdat.truespec[0][0, gdat.indxenerfluxdist[0], 0]

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
    

def pcat_lens_mock_zero():
   
    strgexpo = 0.
    pcat.main.init( \
                   numbswep=10000, \
                   proppsfp=False, \
                   exprinfo=False, \
                   pntstype='lens', \
                   exprtype='hubb', \
                   strgexpo=1e-20, \
                   truenumbpnts=array([40]), \
                  )
    
def pcat_lens_mock():
   
    gdat = pcat.main.init(verbtype=0, pntstype='lens', exprtype='hubb', defa=True, makeplot=False)

    btsr = 0.33
    htsr = 1.

    truespecsour = gdat.hubbexpofact
    truespechost = htsr * truespecsour
    back = [truespecsour * btsr / (pi * gdat.truefixp[gdat.trueindxfixpsizesour]**2) / 2.]
    
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       numbswep=1000, \
                       factthin=100, \
                       minmflux=pi/180./3600.*1e-2, \
                       maxmflux=pi/180./3600.*1e0, \
                       numbswepplot=2000, \
                       proppsfp=False, \
                       exprinfo=False, \
                       pntstype='lens', \
                       exprtype='hubb', \
                       #verbtype=2, \
                       diagmode=False, \
                       back=back, \
                       truenumbpnts=array([40]), \
                       truespecsour=truespecsour, \
                       truespechost=truespechost, \
                      )
    
globals().get(sys.argv[1])()
