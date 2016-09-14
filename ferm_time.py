from __init__ import *
 
def cnfg_null():

    sigm = init( \
                mockfracperd=0., \
                numbiter=100000, \
                verbtype=2, \
               )
    

def mock_grid():

    print 'Initializing grid search...'
    pathimag, pathdata = retr_path()

    sigmthrs = 3.
    listfracperd = logspace(-3., -1., 10)
    #listnumbpuls = logspace(0, 1, 2).astype(int)
    listnumbpuls = arange(1, 10)
    numbfracperd = listfracperd.size
    numbfluxperd = listnumbpuls.size
    fracdete = empty((numbfracperd, numbfluxperd))
    
    # null
    sigm = init( \
                mockfracperd=0., \
                verbtype=2, \
               )
    
    for k in range(listfracperd.size):
        for l in range(listnumbpuls.size):
            sigm = init( \
                        mocknumbpuls=listnumbpuls[l], \
                        mockfracperd=listfracperd[k], \
                        verbtype=2, \
                       )
            numbiter = sigm.size
            numbdete = where(sigm > sigmthrs)[0].size
            fracdete[k, l] = float(numbdete) / numbiter

    # plot the grid
    figr, axis = plt.subplots()
    axis.pcolor(listnumbpuls, listfracperd, fracdete, cmap='Greens')
    axis.set_yscale('log')
    axis.set_xscale('log')
    axis.set_ylabel(r'$\gamma$')
    axis.set_xlabel(r'$\phi$')
    plt.tight_layout()
    path = pathimag + 'mockgrid_%4f.pdf' % sigmthrs
    plt.savefig(path)
    plt.close(figr)

 
def retr_path():
    
    pathbase = os.environ["TDGU_DATA_PATH"]

    pathimag = pathbase + '/imag/ferm_time/'
    os.system('mkdir -p %s' % pathimag)

    pathdata = pathbase + '/data/ferm_time'
    os.system('mkdir -p %s' % pathdata)

    return pathimag, pathdata


def init( \
         datatype='mock', \
         verbtype=1, \
         numbiter=1000, \
         numbphot=1000, \
         numbfreq=1000, \
         minmfreq=5., \
         maxmfreq=15., \
         maxmtime=1000., \
         mockmodltype='sine', \
         mocknumbpuls=1, \
         mockfracperd=0.1, \
         makeplot=True, \
         ):
    
    pathimag, pathdata = retr_path()

    # frequency
    binsfreq = linspace(minmfreq, maxmfreq, numbfreq + 1)
    meanfreq = (binsfreq[1:] + binsfreq[:-1]) / 2.
    indxfreq = arange(numbfreq)
    
    # get data
    if datatype == 'inpt':
        # read Fermi-LAT data
        pass
    else:
        # generate mock data
        if mockmodltype == 'sine':
            # temp
            mockfreq = choice(meanfreq)
    
    # get exposure
    if datatype == 'inpt':
        # read Fermi-LAT exposure
        pass
        numbiter = 1
    else:
        # generate mock exposure
        expo = 1e11
    
    # construct the run tag
    rtag = '%s' % (datatype)
    if datatype == 'mock':
        mockperd = 1. / mockfreq
        mocktimepuls = arange(0., maxmtime - 1., mockperd)
        numbphotperd = int(mockfracperd * numbphot)
        indxphotperd = arange(numbphotperd)
        numbphotrand = numbphot - numbphotperd

        if mockfracperd == 0:
            rtag += '_null'
            if verbtype > 1:
                print 'Null model'
        else:
            rtag += '_%s_%04.f_%04.f' % (mockmodltype, -log10(mockfracperd), -log10(mocknumbpuls))
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
    # plotting
    numbfreqplotrand = 1
    indxfreqplot = concatenate((choice(indxfreq, size=numbfreqplotrand, replace=False), where(meanfreq == mockfreq)[0]))
    
    indxnumbphotperd = arange(numbphotperd)
    
    for k in range(numbiter):
        
        # generate data
        if datatype == 'mock':
            listphotrand = maxmtime * rand(numbphotrand)
            listindxphotperd = array_split(indxphotperd, mocknumbpuls)
            listphotperd = []
            for m in range(mocknumbpuls):
                listphotperd.append(choice(mocktimepuls, size=listindxphotperd[m].size) + 2. * pi * rand())
            listphotperd = concatenate(listphotperd)
            listphot = concatenate((listphotrand, listphotperd))
    
        # take the CTFT
        for n in range(numbfreq):
   
            if k == 0:
                chroinit = time.time()

            phas = meanfreq[n] * listphot
            ctftreal = sin(2. * pi * phas)
            ctftimag = cos(2. * pi * phas)
            ctftrealtotl = sum(ctftreal)
            ctftimagtotl = sum(ctftimag)
            thispsec = (ctftrealtotl**2 + ctftimagtotl**2) / numbphot
            if thispsec > maxmpsec[k]:
                maxmpsec[k] = thispsec
    
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
            print '%2d%% completed.' % fraccomp
            fraccomp = nextfraccomp
    
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
        minm = amin(pval[where(pval > 0.)])
        maxm = amax(pval)
        bins = logspace(log10(minm), log10(maxm), 50)
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

mock_grid()
