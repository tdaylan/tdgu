from __init__ import *

def pcat_lens_mock_grid():
   
    anglfact = 3600. * 180. / pi
    maxmflux = 3e-1 * anglfact
    maxmgang = 2. / anglfact
    mocksizesour = 0.05 / anglfact
    
    numbcnfg = 5
    numbiter = 2
   
    #fact = 1.24510e-19
    fact = 1.63050e-19 # [erg/s/cm^2]
    
    stbr = 10.
    htsr = 10.
    cntssour = 1.

    varbinpt = tdpy.util.varb(numbcnfg)
    varbinpt.defn_para('minmflux', 3e-4, 3e-2, scal='logt')
    varbinpt.defn_para('back', 1e-1 * fact / (pi * mocksizesour**2), 1e1 * fact / (pi * mocksizesour**2), scal='logt')
    varbinpt.defn_para('strgexpo', 1e5 / fact, 1e7 / fact, scal='logt')
    varbinpt.defn_para('mockfluxsour', 1e-1 * fact, 1e1 * fact, scal='logt')
    varbinpt.defn_para('mockfluxhost', 1e-1 * fact, 1e1 * fact, scal='logt')

    listnameoutpvarb = ['fluxdistsloppop0', 'meanpntspop0', 'specassc', 'fluxhost', 'beinhost']
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
                    if varbinpt.name[p] == 'mockfluxsour':
                        mockfluxsour = varb
                    if varbinpt.name[p] == 'mockfluxhost':
                        mockfluxhost = varb
                
                gdat = pcat.main.init( \
                                      seedstat=seedstat, \
                                      numbswep=100000, \
                                      numbswepplot=4000, \
                                      factthin=200, \
                                      proppsfp=False, \
                                      exprinfo=False, \
                                      pntstype='lens', \
                                      indxenerincl=arange(1), \
                                      indxevttincl=arange(1), \
                                      strgexpo=strgexpo, \
                                      exprtype='hubb', \
                                      back=[back], \
                                      minmflux=minmflux, \
                                      maxmflux=maxmflux, \
                                      mockfluxsour=mockfluxsour, \
                                      mockfluxhost=mockfluxhost, \
                                      mocknumbpnts=array([10]), \
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
    

def pcat_lens_mock():
   
    maxmflux = deg2rad(5e-1 / 3600.)
    minmflux = deg2rad(1e-3 / 3600.)
    
    #fact = 1.24510e-19
    fact = 1.63050e-19
    
    anglfact = 3600. * 180. / pi
    maxmgang = 2. / anglfact
    mocksizesour = 0.05 / anglfact

    stbr = 10.
    htsr = 10.
    cntssour = 1.
    strgexpo = 1e4
    
    mockfluxsour = cntssour * fact
    mockfluxhost = htsr * mockfluxsour
    back = [mockfluxsour / stbr / (pi * mocksizesour**2)]
    
    strgexpo /= fact
    
    numbiter = 10
    for k in range(numbiter):
        pcat.main.init( \
                       numbswep=100000, \
                       numbswepplot=10000, \
                       factthin=10, \
                       verbtype=1, \
                       maxmgang=maxmgang, \
                       proppsfp=False, \
                       exprinfo=False, \
                       pntstype='lens', \
                       indxenerincl=arange(1), \
                       indxevttincl=arange(1), \
                       strgexpo=strgexpo, \
                       exprtype='hubb', \
                       back=back, \
                       minmflux=minmflux, \
                       maxmflux=maxmflux, \
                       #maxmnumbpnts=array([3]), \
                       mocknumbpnts=array([20]), \
                       mockfluxsour=mockfluxsour, \
                       mockfluxhost=mockfluxhost, \
                      )
    
globals().get(sys.argv[1])()
