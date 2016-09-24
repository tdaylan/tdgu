from __init__ import *
 
def cnfg_mocknull():

    sigm = init( \
                mockfracperd=0., \
                numbiter=100000, \
                verbtype=2, \
               )


def cnfg_inpt():

    sigm = init( \
                datatype='inpt', \
                numbiter=100000, \
                verbtype=2, \
               )


def cnfg_mockgrid():

    print 'Initializing grid search...'
    pathimag, pathdata = tdpy.util.retr_path('ferm_time')

    sigmthrs = 3.
    minmfracperd = 0.01
    maxmfracperd = 0.3
    numbfracperd = 20
    listfracperd = logspace(log10(minmfracperd), log10(maxmfracperd), numbfracperd)
    
    minmnumbpuls = 1
    maxmnumbpuls = 6
    numbnumbpuls = 6
    listnumbpuls = linspace(minmnumbpuls, maxmnumbpuls, numbnumbpuls).astype(int)
    
    numbfracperd = listfracperd.size
    numbfluxperd = listnumbpuls.size

    path = pathdata + 'fracdete%04f%04f%04f%04f%04f%04f%04f.fits' % (sigmthrs, -log10(minmfracperd), -log10(maxmfracperd), numbfracperd, minmnumbpuls, maxmnumbpuls, numbnumbpuls)
    if os.path.isfile(path):
        print 'Reading from %s...' % path
        fracdete = pf.getdata(path)
    else:
        fracdete = empty((numbfracperd, numbfluxperd))
        for k in range(listfracperd.size):
            for l in range(listnumbpuls.size):
                tdpy.util.show_memo_simp()
                sigm = init( \
                            mocknumbpuls=listnumbpuls[l], \
                            mockfracperd=listfracperd[k], \
                            numbiter=10, \
                            verbtype=2, \
                           )
                numbiter = sigm.size
                numbdete = where(sigm > sigmthrs)[0].size
                fracdete[k, l] = float(numbdete) / numbiter
        pf.writeto(path, fracdete, clobber=True)

    # plot the grid
    figr, axis = plt.subplots()
    imag = axis.pcolor(listnumbpuls, listfracperd, fracdete, cmap='Greens')
    axis.set_yscale('log')
    axis.set_xlim([minmnumbpuls, maxmnumbpuls])
    axis.set_ylim([minmfracperd, maxmfracperd])
    axis.set_ylabel(r'$\gamma$')
    axis.set_xlabel(r'$N$')
    cbar = plt.colorbar(imag) 
    plt.tight_layout()
    path = pathimag + 'mockgrid_%4f.pdf' % sigmthrs
    plt.savefig(path)
    plt.close(figr)

 
def init( \
         datatype='mock', \
         verbtype=1, \
         numbiter=100, \
         numbfreq=1000, \
         minmfreq=5., \
         maxmfreq=15., \
         maxmtime=1000., \
         stattype='ctft', \
         mockmodltype='sine', \
         mocknumbphot=1000, \
         mocknumbpuls=1, \
         mockfracperd=0.1, \
         makeplot=True, \
         ):
    
    pathimag, pathdata = tdpy.util.retr_path('ferm_time')

    # axes
    ## frequency
    binsfreq, meanfreq, numbfreq, indxfreq = tdpy.util.retr_axis(minmfreq, maxmfreq, numbfreq)
    ## period
    binsperd = 1. / binsfreq
    meanperd = 1. / meanfreq
    numbperd = numbfreq
    indxperd = indxfreq

    # get the time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()

    # get data
    if datatype == 'inpt':
        # read Fermi-LAT data
        path = pathdata + '/fermdata.fits'
        if os.path.isfile(path):
            print 'Reading from %s...' % path
            listphot = pf.getdata(path)
        else:

            rasc = 17. / 24. * 360. + 45. / 60. + 40.04 / 3600.
            decl = -29. + 28.1 / 3600.

            #cmnd = 'ls -d -1 $FERMI_DATA/weekly/spacecraft/* >> spactemp.txt' + spac
            #os.system(cmnd)
            pathinpt = pathdata + 'phottemp.txt'
            strgproc = os.uname()[1]
            if strgproc == 'fink1':
                strg = '$FERMI_DATA/weekly/photon/'
            else:
                strg = pathdata
            cmnd = 'ls -d -1 %s*_w100_* >> %s' % (strg, pathinpt)
            os.system(cmnd)

            cmnd = 'gtselect infile=' + pathinpt + ' outfile=' + path + ' ra=%.3g dec=%.3g rad=0.1 emin=1000 emax=3000 zmax=90 evclass=128 evtype=32' % (rasc, decl)
            os.system(cmnd)

            listphot = pf.getdata(path)
            print 'listphot'
            print listphot
            
        numbphot = listphot.size

    else:
        # generate mock data
        if mockmodltype == 'sine':
            # temp
            mockfreq = choice(meanfreq)

        numbphot = mocknumbphot
    
    # get exposure
    if datatype == 'inpt':
        # read Fermi-LAT exposure
        pass
        numbiter = 1
        indxfreqplot = choice(indxfreq, size=1, replace=False)
    else:
        # generate mock exposure
        expo = 1e11
        indxfreqplot = concatenate((choice(indxfreq, size=1, replace=False), where(meanfreq == mockfreq)[0]))
    
    # construct the run tag
    rtag = '%s_%s' % (strgtimestmp, datatype)
    if datatype == 'mock':
        mockperd = 1. / mockfreq
        mocktimepuls = arange(0., maxmtime - 1., mockperd)
        mocknumbphotperd = int(mockfracperd * mocknumbphot)
        indxphotperd = arange(mocknumbphotperd)
        mocknumbphotrand = mocknumbphot - mocknumbphotperd

        if mockfracperd == 0:
            rtag += '_null'
            if verbtype > 1:
                print 'Null model'
        else:
            rtag += '_%s_%04.f_%04.f' % (mockmodltype, -log10(mockfracperd), log10(mocknumbpuls))
            if verbtype > 1:
                print 'Model %s' % mockmodltype
                print 'mockfracperd'
                print mockfracperd
                print 'mocknumbpuls'
                print mocknumbpuls
            if mockmodltype == 'sine':
                rtag += '_%04.f' % (mockfreq)
            
    minmtime = 0.
    
    numbtime = 100
    binstime = linspace(minmtime, maxmtime, numbtime + 1)
    meantime = (binstime[1:] + binstime[:-1]) / 2.
    
    # initialize the maximum power spectrum
    maxmpsec = zeros(numbiter)
    if makeplot:
        listctftrealtotl = empty(numbfreq)
        listctftimagtotl = empty(numbfreq)
        listpsec = empty(numbfreq)
    fraccomp = 0.
    
    listchrofreq = empty(numbfreq)
    
    for k in range(numbiter):
        
        # generate data
        if datatype == 'mock':
            listphotrand = maxmtime * rand(mocknumbphotrand)
            listindxphotperd = array_split(indxphotperd, mocknumbpuls)
            listphotperd = []
            for m in range(mocknumbpuls):
                listphotperd.append(choice(mocktimepuls, size=listindxphotperd[m].size) + 2. * pi * rand())
            listphotperd = concatenate(listphotperd)
            listphot = concatenate((listphotrand, listphotperd))
    
        # search over the frequency grid
        for n in range(numbfreq):
   
            if k == 0:
                chroinit = time.time()

            if stattype == 'ctft':
                
                # use the power spectrum as the test statistics
                ## take the CTFT
                phas = 2. * pi * meanfreq[n] * listphot
                ctftreal = sin(phas)
                ctftimag = cos(phas)
                ctftrealtotl = sum(ctftreal)
                ctftimagtotl = sum(ctftimag)
                thispsec = (ctftrealtotl**2 + ctftimagtotl**2) / numbphot
                if thispsec > maxmpsec[k]:
                    maxmpsec[k] = thispsec
            else:
                # use the wrapped light curve as the test statistics
                bindphot = histogram(listphot % meanperd[n])

            if makeplot:
                listctftrealtotl[n] = ctftrealtotl
                listctftimagtotl[n] = ctftimagtotl
                listpsec[n] = thispsec
       
                if k == 0 and in1d(n, indxfreqplot):
                
                    strgfreq = '%04g' % (meanfreq[n])
                    figr, axis = plt.subplots()
                    axis.hist(listphot % (1. / meanfreq[n]))
                    axis.set_xlabel(r'$\bar{t}$')
                    plt.tight_layout()
                    path = pathimag + 'histwrap_%s.pdf' % rtag
                    plt.savefig(path)
                    plt.close(figr)
                    
                    figr, axis = plt.subplots()
                    axis.hist(diff(sort(listphot)))
                    axis.set_xlabel(r'$\Delta t$')
                    plt.tight_layout()
                    path = pathimag + 'histdiff_%s.pdf' % rtag
                    plt.savefig(path)
                    plt.close(figr)
                    
                    minm = min(amin(ctftreal), amin(ctftimag))
                    maxm = max(amax(ctftreal), amax(ctftimag))
                    bins = linspace(minm, maxm, 50)
                    figr, axis = plt.subplots()
                    axis.hist(ctftreal, bins=bins, alpha=0.3, label='Re')
                    axis.hist(ctftimag, bins=bins, alpha=0.3, label='Im')
                    axis.set_xlabel('$a_f$')
                    axis.legend(loc=9)
                    plt.tight_layout()
                    path = pathimag + 'histctft_%s.pdf' % rtag
                    plt.savefig(path)
                    plt.close(figr)
    
                    figr, axis = plt.subplots()
                    axis.hist(phas)
                    axis.set_xlabel(r'$\phi$')
                    plt.tight_layout()
                    path = pathimag + 'histphas_%s.pdf' % rtag
                    plt.savefig(path)
                    plt.close(figr)
    
            if k == 0:
                chrofinl = time.time()
                listchrofreq[n] = chrofinl - chroinit

        nextfraccomp = int(100. * k / numbiter)
        if nextfraccomp % 10 == 0 and fraccomp < nextfraccomp and verbtype > 0:
            fraccomp = nextfraccomp
            print '%2d%% completed.' % fraccomp
    
        if makeplot and k == 0:
    
            # bin data
            cntsback = histogram(listphot, bins=binstime)[0] 
            
            figr, axis = plt.subplots()
            axis.plot(meantime, cntsback, marker='o', ls='')
            axis.set_xlabel('$t$ [s]')
            axis.set_ylabel('$N_\gamma$')
            plt.tight_layout()
            path = pathimag + 'cntsback_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)
    
            figr, axis = plt.subplots()
            axis.plot(meanfreq, listpsec)
            axis.set_xlabel('$f$ [Hz]')
            axis.set_ylabel('$S(f)$ [1/Hz]')
            plt.tight_layout()
            path = pathimag + 'psec_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)
    
            figr, axis = plt.subplots()
            bins = logspace(amin(log10(listchrofreq)), amax(log10(listchrofreq)), 20)
            axis.hist(listchrofreq * 1e3, bins=bins)
            axis.set_xscale('log')
            axis.set_yscale('log')
            axis.set_xlabel('$t$ [ms]')
            plt.tight_layout()
            path = pathimag + 'histchrofreq_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)
    
            figr, axis = plt.subplots()
            bins = linspace(amin(listpsec), amax(listpsec), 50)
            delt = bins[1] - bins[0]
            axis.hist(listpsec, bins=bins)
            axis.set_yscale('log')
            para = sp.stats.expon.fit(listpsec)
            axis.plot(bins, numbfreq * delt * sp.stats.expon.pdf(bins, loc=para[0], scale=para[1]))
            plt.tight_layout()
            path = pathimag + 'histpsec_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)
           
            figr, axis = plt.subplots()
            minm = min(amin(listctftrealtotl), amin(listctftimagtotl))
            maxm = max(amax(listctftrealtotl), amax(listctftimagtotl))
            bins = linspace(minm, maxm, 50)
            axis.hist(listctftrealtotl, bins=bins, alpha=0.3, label='Re')
            axis.hist(listctftimagtotl, bins=bins, alpha=0.3, label='Im')
            axis.legend()
            plt.tight_layout()
            path = pathimag + 'histctfttotl_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)
    
    if numbiter > 1:
        figr, axis = plt.subplots()
        bins = linspace(amin(maxmpsec), amax(maxmpsec), 50)
        delt = diff(bins)[0]
        axis.hist(maxmpsec, bins=bins)
        
        axis.plot(bins, numbiter * delt * sp.stats.gumbel_r.pdf(bins, loc=log(numbfreq), scale=1.), label='Gumbel Model')
        para = sp.stats.gumbel_r.fit(maxmpsec)
        axis.plot(bins, numbiter * delt * sp.stats.gumbel_r.pdf(bins, loc=para[0], scale=para[1]), label='Gumbel Fit')
        axis.set_xlabel('max$(S)$')
        axis.set_ylabel('$N$')
        axis.legend()
        plt.tight_layout()
        path = pathimag + 'maxmpsec_%s.pdf' % rtag
        plt.savefig(path)
        plt.close(figr)
        
        pval = sp.stats.gumbel_r.sf(maxmpsec, loc=log(numbfreq), scale=1.)
        indxpvalzero = where(pval == 0.)[0]
        if indxpvalzero.size > 0:
            print 'All p values are zero!'
            plotpval = False
            sigm = zeros_like(pval) + 1000.
        else:
            minm = amin(pval[where(pval > 0.)])
            maxm = amax(pval)
            bins = logspace(log10(minm), log10(maxm), 50)
            plotpval = True

        if plotpval:
            figr, axis = plt.subplots()
            axis.hist(pval, bins=bins)
            axis.set_xlabel('$p$ value')
            axis.set_xscale('log')
            axis.set_ylabel('$N$')
            plt.tight_layout()
            path = pathimag + 'maxmpsecpval_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)
            
            figr, axis = plt.subplots()
            sigm = sp.stats.norm.ppf(1. - pval)
            axis.hist(sigm)
            axis.set_xlabel(r'Significance [$\sigma$]')
            axis.set_ylabel('$N$')
            plt.tight_layout()
            path = pathimag + 'maxmpsecsigm_%s.pdf' % rtag
            plt.savefig(path)
            plt.close(figr)

    if verbtype > 1:
        print '%d iterations perfomed.' % numbiter

    return sigm


if len(sys.argv) > 1:
    name = globals().copy()
    name.get(sys.argv[1])()
else:
    pass

