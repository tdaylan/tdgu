from __init__ import *

def pcat_lens_mock_grid():
   
    anglfact = pi / 180. / 3600.
    maxmflux = 3e-1 * anglfact
    
    numbcnfg = 5
    
    varbinpt = tdpy.util.varb(numbcnfg)
    varbinpt.defn_para('minmflux', 3e-4, 3e-2, scal='logt')
    varbinpt.defn_para('strgback', 1e-2, 1e0, scal='logt')
    varbinpt.defn_para('strgexpo', 1e15, 1e17, scal='logt')

    listnameoutpvarb = ['fluxdistslop', 'meanpnts', 'specassc']
    numboutpvarb = len(listnameoutpvarb)
    liststrgoutpvarb = []

    liststrgvarbinpt = ['$R_{min}$ [arcsec]', '$A$ [1/cm^2/s]', r'$\epsilon$ [cm^2 s]']
    grid = empty((4, numboutpvarb, varbinpt.size, numbcnfg))
    
    numbiter = 2
    for k in range(numbiter):
        seedstat = get_state()
        for m in range(varbinpt.size):
            for l in range(varbinpt.numb):
                
                if m > 0 and l == varbinpt.numb / 2:
                    grid[:, :, m, l] = grid[:, :, 0, varbinpt.numb / 2]
                    continue

                for p in range(varbinpt.size):
                    if p == m:
                        varb = varbinpt.para[p][l]
                    else:
                        varb = varbinpt.para[p][numbcnfg / 2]
                     
                    if varbinpt.name[p] == 'minmflux':
                        minmflux = varb * anglfact
                    if varbinpt.name[p] == 'strgback':
                        strgback = varb
                    if varbinpt.name[p] == 'strgexpo':
                        strgexpo = varb
                
                gridchan, gdat = pcat.main.init( \
                                                seedstat=seedstat, \
                                                numbswep=10000, \
                                                numbswepplot=1000, \
                                                factthin=180, \
                                                proplenp=False, \
                                                #makeplot=False, \
                                                diagmode=False, \
                                                proppsfp=False, \
                                                exprinfo=False, \
                                                #probtran=0., \
                                                #proplenp=False, \
                                                pntstype='lens', \
                                                indxenerincl=arange(1), \
                                                indxevttincl=arange(1), \
                                                strgexpo=strgexpo, \
                                                exprtype='hubb', \
                                                strgback=[strgback], \
                                                #maxmnumbpnts=array([20]), \
                                                minmflux=minmflux, \
                                                maxmflux=maxmflux, \
                                                mocknumbpnts=array([10]), \
                                               )
                
                for n in range(numboutpvarb):
                    if listnameoutpvarb[n] == 'specassc':
                        grid[:3, n, m, l] = gdat.fluxfactplot * gdat.postspecassc[:, gdat.indxenerfluxdist[0], 0]
                        grid[3, n, m, l] = gdat.fluxfactplot * gdat.truespec[0][0, gdat.indxenerfluxdist[0], 0]
                        if k == 0 and m == 0 and l == 0:
                            liststrgoutpvarb.append(['$R_0$'])
                    else:
                        grid[:3, n, m, l] = getattr(gdat, 'postfixp')[:, getattr(gdat, 'indxfixp' + listnameoutpvarb[n])].flatten()
                        grid[3, n, m, l] = getattr(gdat, 'truefixp')[getattr(gdat, 'indxfixp' + listnameoutpvarb[n])]
                        if k == 0 and m == 0 and l == 0:
                            liststrgoutpvarb.append(getattr(gdat, 'strgfixp')[getattr(gdat, 'indxfixp' + listnameoutpvarb[n])])
        
        path = os.environ["PCAT_DATA_PATH"] + '/imag/%s_lensgrid/' % gdat.strgtimestmp
        os.system('mkdir -p %s' % path)
        for n in range(numboutpvarb):
            for m in range(varbinpt.size):
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                yerr = tdpy.util.retr_errrvarb(grid[:3, n, m, :])
                axis.errorbar(varbinpt.para[m], grid[0, n, m, :], yerr=yerr, ls='', marker='o')
                axis.plot(varbinpt.para[m], grid[3, n, m, :], marker='x', color='g')
                axis.set_xlabel(liststrgvarbinpt[m])
                axis.set_ylabel(liststrgoutpvarb[n][0])
                
                maxm = amax(varbinpt.para[m])
                minm = amin(varbinpt.para[m])
                if varbinpt.scal[m] == 'logt':
                    axis.set_xscale('log')
                    axis.set_xlim([minm / 2., maxm * 2.])
                else:
                    axis.set_xlim([minm - 1., maxm + 1.])
                plt.tight_layout()
                plt.savefig('%s/%d%d%d.pdf' % (path, n, m, k))
                plt.close(figr)
    

def pcat_lens_mock():
   
    maxmflux = deg2rad(5e-1 / 3600.)
    minmflux = deg2rad(1e-3 / 3600.)
    
    numbiter = 10
    for k in range(numbiter):
        gridchan = pcat.main.init( \
                                  numbswep=100, \
                                  numbswepplot=30000, \
                                  factthin=1, \
                                  numbburn=0, \
                                  verbtype=1, \
                                  #diagmode=False, \
                                  proppsfp=False, \
                                  #proplenp=False, \
                                  #diagmode=False, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  strgback=[1e-1], \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  maxmnumbpnts=array([30]), \
                                  mocknumbpnts=array([10]), \
                                 )
    
globals().get(sys.argv[1])()
