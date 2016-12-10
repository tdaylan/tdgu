from __init__ import *


def pcat_lens_mock_arry_depr():

    #liststrgback = linspace(0.1, 3., 3.)
    liststrgback = logspace(-1., 1., 3.)

    minmflux = deg2rad(5e-2 / 3600.)
    maxmflux = deg2rad(5e-1 / 3600.)

    for k in range(len(liststrgback)):
       
        strgback = [liststrgback[k]]
        gridchan = pcat.main.init( \
                                  numbswep=100000, \
                                  numbswepplot=10000, \
                                  factthin=100, \
                                  diagmode=False, \
                                  optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  strgback=strgback, \
                                  maxmnumbpnts=array([10]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mocknumbpnts=array([5]), \
                                 )

    
def pcat_lens_mock_arry():
   
    anglfact = pi / 180. / 3600.
    maxmflux = 3e-1 * anglfact
    
    numbcnfg = 3
    
    varbinpt = tdpy.util.varb(numbcnfg)
    varbinpt.defn_para('minmflux', 1e-3, 1e-1, scal='logt')
    varbinpt.defn_para('strgback', 1e-2, 1e0, scal='logt')
    varbinpt.defn_para('strgexpo', 1e15, 1e17, scal='logt')

    listnameoutpvarb = ['fluxdistslop', 'meanpnts']
    numboutpvarb = len(listnameoutpvarb)
    liststrgoutpvarb = []

    liststrgvarbinpt = ['$R_{min}$ [arcsec]', '$A$ [1/cm^2/s]', r'$\epsilon$ [cm^2 s]']
    arry = empty((4, numboutpvarb, varbinpt.size, numbcnfg))
    
    numbiter = 1
    for k in range(numbiter):
        seedstat = get_state()
        for m in range(varbinpt.size):
            for l in range(varbinpt.numb):
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
                                                numbswep=1000, \
                                                numbswepplot=3000, \
                                                factthin=10, \
                                                diagmode=False, \
                                                proppsfp=False, \
                                                exprinfo=False, \
                                                probtran=0., \
                                                proplenp=False, \
                                                pntstype='lens', \
                                                indxenerincl=arange(1), \
                                                indxevttincl=arange(1), \
                                                strgexpo=strgexpo, \
                                                exprtype='hubb', \
                                                strgback=[strgback], \
                                                maxmnumbpnts=array([20]), \
                                                minmflux=minmflux, \
                                                maxmflux=maxmflux, \
                                                mocknumbpnts=array([10]), \
                                               )
                
                for n in range(numboutpvarb):
                    arry[:3, n, m, l] = getattr(gdat, 'postfixp')[:, getattr(gdat, 'indxfixp' + listnameoutpvarb[n])].flatten()
                    arry[3, n, m, l] = getattr(gdat, 'truefixp')[getattr(gdat, 'indxfixp' + listnameoutpvarb[n])]
                    print 'arry[3, n, m, l]'
                    print arry[3, n, m, l]
                    print 
                    if k == 0 and m == 0 and l == 0:
                        liststrgoutpvarb.append(getattr(gdat, 'strgfixp')[getattr(gdat, 'indxfixp' + listnameoutpvarb[n])])
        
        for n in range(numboutpvarb):
            for m in range(varbinpt.size):
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                yerr = tdpy.util.retr_errrvarb(arry[:3, n, m, :])
                axis.errorbar(varbinpt.para[m], arry[0, n, m, :], yerr=yerr, ls='', marker='o')
                axis.plot(varbinpt.para[m], arry[3, n, m, :], ls='', marker='x', color='g')
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
                plt.savefig(os.environ["PCAT_DATA_PATH"] + '/imag/lensarry%d%d%d.pdf' % (n, m, k))
                plt.close(figr)
    

def pcat_lens_mock():
    
    maxmflux = deg2rad(5e-1 / 3600.)
    
    numbcnfg = 3
    listminmflux = logspace(-2., -1., numbcnfg)
    liststrgback = logspace(-3., 0., numbcnfg)
    listexpo = logspace(-3., 0., numbcnfg)

    numbiter = 1
    for k in range(numbiter):
        seedstat = get_state()
        for a in range(len(listminmflux)):
            for b in range(len(liststrgback)):
                minmflux = listminmflux[a]
                strgback = [liststrgback[b]]
                strgexpo = listexpo[c]
                gridchan = pcat.main.init( \
                                  #verbtype=2, \ 
                                  seedstat=seedstat, \
                                  #makeplot=False, \
                                  numbswep=10000, \
                                  numbswepplot=3000, \
                                  factthin=10, \
                                  diagmode=False, \
                                  probtran=1., \
                                  #prophypr=False, \
                                  proppsfp=False, \
                                  #propbacp=False, \
                                  #proplenp=False, \
                                  #propcomp=False, \
                                  #stdvprophypr=0.1, \
                                  #stdvpropbacp=5e-3, \
                                  #stdvproplenp=5e-4, \
                                  #probtran=0., \
                                  #optiprop=True, \
                                  exprinfo=False, \
                                  pntstype='lens', \
                                  indxenerincl=arange(1), \
                                  indxevttincl=arange(1), \
                                  strgexpo=1e16, \
                                  exprtype='hubb', \
                                  strgback=strgback, \
                                  maxmnumbpnts=array([20]), \
                                  minmflux=minmflux, \
                                  maxmflux=maxmflux, \
                                  mocknumbpnts=array([10]), \
                                 )

        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.imshow(evid, extent=[minmgain, maxmgain, minmdevi, maxmdevi], cmap='winter', origin='lower')
        cset1 = plt.contourf(gain, devi, evid, cmap='winter')
        axis.set_xlabel('Exposure')
        axis.set_ylabel('Background')
        axis.set_zlabel('SubHalo Mass')
        plt.colorbar(imag, ax=axis, fraction=0.03)
        plt.tight_layout()
        plt.savefig(gdat.pathplot + 'lensgrid%d.pdf' % k)
        plt.close(figr)
    
    
globals().get(sys.argv[1])()
