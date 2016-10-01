from __init__ import *
 
def cnfg_mocknull():

    sigm = init( \
                stattype='maxmkosm', \
                mockfracperd=0., \
                numbiter=10000, \
                verbtype=2, \
               )


def cnfg_inpt():

    sigm = init( \
                datatype='inpt', \
                numbiter=100000, \
                verbtype=2, \
               )


def make_maps():

    # make time series
    ## Vela
    rasc, decl = tdpy.util.conv_rascdecl(8, 35, 20.655, -45, 10, 35.15)
    ## galactic center

    temp, pathdata = tdpy.util.retr_path('tdgu', 'ferm_time/', 'ferm_time/', '')


    path = pathdata + 'listphot'
    if os.path.isfile(path):
        print 'Reading %s...' % path
        listphot = pf.getdata(path)
    else:

        path = 'gtbary '
        os.system(path)
        
        path = 'gtbary '
        os.system(path)
        
        listphot = pf.getdata(path)
        


def cnfg_mockgrid():

    print 'Initializing grid search...'
    pathimag, pathdata = tdpy.util.retr_path('tdgu', 'ferm_time/', 'ferm_time/', '')

    sigmthrs = 3.
    minmfracperd = 0.01
    maxmfracperd = 0.3
    numbfracperd = 2
    listfracperd = logspace(log10(minmfracperd), log10(maxmfracperd), numbfracperd)
    
    minmnumbpuls = 1
    maxmnumbpuls = 6
    numbnumbpuls = 3
    listnumbpuls = linspace(minmnumbpuls, maxmnumbpuls, numbnumbpuls).astype(int)
    
    numbfracperd = listfracperd.size
    numbfluxperd = listnumbpuls.size

    #stattype = 'maxmpsec'
    stattype = 'maxmkosm'

    path = pathdata + 'fracdete_%s_%04f_%04f_%04f_%04f_%04f_%04f_%04f.fits' % (stattype, sigmthrs, -log10(minmfracperd), -log10(maxmfracperd), numbfracperd, \
                                                                                                                            minmnumbpuls, maxmnumbpuls, numbnumbpuls)
    if os.path.isfile(path):
        print 'Reading %s...' % path
        fracdete = pf.getdata(path)
    else:
        fracdete = empty((numbfracperd, numbfluxperd))
        for k in range(listfracperd.size):
            for l in range(listnumbpuls.size):
                tdpy.util.show_memo_simp()
                sigm = init( \
                            stattype=stattype, \
                            mocknumbpuls=listnumbpuls[l], \
                            mockfracperd=listfracperd[k], \
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
    path = pathimag + 'mockgrid_%s%04.f.pdf' % (stattype, sigmthrs)
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
         stattype='maxmpsec', \
         difftimewrap=1e-5, \
         mockmodltype='sine', \
         mocknumbphot=1000, \
         mocknumbpuls=1, \
         mockfracperd=0.1, \
         makeplot=True, \
         ):
   
    # get the time stamp
    strgtimestmp = tdpy.util.retr_strgtimestmp()
    
    # construct the run tag
    rtag = '%s_%s_%s' % (strgtimestmp, stattype, datatype)
    
    # axes
    ## frequency
    binsfreq, meanfreq, difffreq, numbfreq, indxfreq = tdpy.util.retr_axis(minmfreq, maxmfreq, numbfreq)
    ## period
    binsperd = 1. / binsfreq
    meanperd = 1. / meanfreq
    numbperd = numbfreq
    indxperd = indxfreq

    # get data
    if datatype == 'inpt':
        # read Fermi-LAT data
        path = pathdata + '/fermdata.fits'
        if os.path.isfile(path):
            print 'Reading %s...' % path
            listtime = pf.getdata(path)
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

            listtime = pf.getdata(path)
            print 'listtime'
            print listtime
            
        numbphot = listtime.size

    else:
        # generate mock data
        if mockmodltype == 'sine':
            # temp
            mockfreq = choice(meanfreq)

        numbphot = mocknumbphot
   
    indxphot = arange(numbphot, dtype=float) / (numbphot - 1.)

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
            
    # paths
    pathimag, pathdata = tdpy.util.retr_path('tdgu', 'ferm_time/', 'ferm_time/')

    minmtime = 0.
    numbtime = 100
    binstime, meantime, difftime, numbtime, indxtime = tdpy.util.retr_axis(minmtime, maxmtime, numbtime)
    
    fraccomp = 0.
    
    listchrofreq = empty(numbfreq)
    
    if stattype == 'maxmpsec':
        maxmpsec = zeros(numbiter)
        if makeplot:
            listctftrealtotl = empty(numbfreq)
            listctftimagtotl = empty(numbfreq)
            listpsec = empty(numbfreq)
    if stattype == 'maxmkosm':
        maxmkosm = zeros(numbiter)
        if makeplot:
            kosm = empty(numbfreq)
            pvalkosm = empty(numbfreq)
    if stattype == 'maxmcnts':
        numbtimewrap = int(timetotl / difftimewrap)
        bindtimewrapmean = numphot / float(numbtimewrap)

    for k in range(numbiter):
        
        # generate data
        if datatype == 'mock':
            listtimerand = maxmtime * rand(mocknumbphotrand)
            listindxphotperd = array_split(indxphotperd, mocknumbpuls)
            listtimeperd = []
            for m in range(mocknumbpuls):
                listtimeperd.append(choice(mocktimepuls, size=listindxphotperd[m].size) + 2. * pi * rand())
            listtimeperd = concatenate(listtimeperd)
            listtime = concatenate((listtimerand, listtimeperd))
    
        # search over the frequency grid
        for n in range(numbfreq):
   
            if k == 0:
                chroinit = time.time()

            if stattype == 'maxmpsec':
                # use the power spectrum as the test statistics
                ## take the CTFT
                phas = 2. * pi * meanfreq[n] * listtime
                ctftreal = sin(phas)
                ctftimag = cos(phas)
                ctftrealtotl = sum(ctftreal)
                ctftimagtotl = sum(ctftimag)
                ## find the power spectrum
                thispsec = (ctftrealtotl**2 + ctftimagtotl**2) / numbphot
                if thispsec > maxmpsec[k]:
                    maxmpsec[k] = thispsec
            else:
                # use the wrapped light curve as the test statistics
                listtimewrap = (listtime % meanperd[n]) / meanperd[n]
                
                if stattype == 'maxmkosm':
                    listtimewrapsort = sort(listtimewrap)
                    kosm[n] = amax(fabs(indxphot - listtimewrapsort))
                    pvalkosm[n] = 0.

                    if kosm[n] > maxmkosm[k]:
                        maxmkosm[k] = kosm[n]
                if stattype == 'bind':
                    bindtimewrap = histogram(listtimewrap, bins=binstimewrap)[0]
                    pval[n] = sp.stats.poisson.sf(bindtimewrap, lam=bindtimewrapmean[n, :])

            if makeplot:
                if stattype == 'maxmpsec':
                    listctftrealtotl[n] = ctftrealtotl
                    listctftimagtotl[n] = ctftimagtotl
                    listpsec[n] = thispsec
           
                if k == 0 and in1d(n, indxfreqplot):
                    
                    if stattype == 'maxmpsec':
                        minm = min(amin(ctftreal), amin(ctftimag))
                        maxm = max(amax(ctftreal), amax(ctftimag))
                        bins = linspace(minm, maxm, 50)
                        figr, axis = plt.subplots()
                        axis.hist(ctftreal, bins=bins, alpha=0.3, label='Re')
                        axis.hist(ctftimag, bins=bins, alpha=0.3, label='Im')
                        axis.set_xlabel('$a_f$')
                        axis.legend(loc=9)
                        plt.tight_layout()
                        path = pathimag + 'ctfthist_%04d.pdf' % n
                        plt.savefig(path)
                        plt.close(figr)
    
                    if stattype == 'maxmcnts':
                        strgfreq = '%04g' % (meanfreq[n])
                        figr, axis = plt.subplots()
                        axis.hist(listtimewrap)
                        axis.set_xlabel(r'$\bar{t}$')
                        plt.tight_layout()
                        path = pathimag + 'timewraphist_%04d.pdf' % n
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
            cntsback = histogram(listtime, bins=binstime)[0] 
            
            figr, axis = plt.subplots()
            axis.plot(meantime, cntsback, marker='o', ls='')
            axis.set_xlabel('$t$ [s]')
            axis.set_ylabel('$N_\gamma$')
            plt.tight_layout()
            path = pathimag + 'cntsback.pdf'
            plt.savefig(path)
            plt.close(figr)
    
            if stattype == 'maxmkosm':
                figr, axis = plt.subplots()
                axis.hist(kosm, bins=linspace(amin(kosm), amax(kosm), 50))
                axis.set_xlabel(r'TS$_{KS}$')
                plt.tight_layout()
                path = pathimag + 'kosmhist.pdf'
                plt.savefig(path)
                plt.close(figr)
                
                if False:
                    figr, axis = plt.subplots()
                    axis.hist(pvalkosm)
                    axis.set_xlabel(r'$p_{KS}$')
                    plt.tight_layout()
                    path = pathimag + 'pvalkosmhist.pdf'
                    plt.savefig(path)
                    plt.close(figr)

            if stattype == 'maxmpsec':
                figr, axis = plt.subplots()
                axis.plot(meanfreq, listpsec)
                axis.set_xlabel('$f$ [Hz]')
                axis.set_ylabel('$S(f)$ [1/Hz]')
                plt.tight_layout()
                path = pathimag + 'psec.pdf'
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
                path = pathimag + 'psechist.pdf'
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
                path = pathimag + 'ctfttotlhist.pdf'
                plt.savefig(path)
                plt.close(figr)
            
            figr, axis = plt.subplots()
            bins = logspace(amin(log10(listchrofreq * 1e3)), amax(log10(listchrofreq * 1e3)), 20)
            axis.hist(listchrofreq * 1e3, bins=bins)
            axis.set_xscale('log')
            axis.set_yscale('log')
            axis.set_xlabel('$t$ [ms]')
            plt.tight_layout()
            path = pathimag + 'chrofreqhist.pdf'
            plt.savefig(path)
            plt.close(figr)
    
    # calculate p value
    if stattype == 'maxmpsec':
        pval = sp.stats.gumbel_r.sf(maxmpsec, loc=log(numbfreq), scale=1.)
        stat = maxmpsec
    if stattype == 'maxmkosm':
        # retrieve the polynomial fit to the distribution of the maximum KS TS
        path = pathdata + 'poly.npz'
        if os.path.isfile(path):
            print 'Reading %s...' % path
            temp = load(path)
            coefmaxmkosm = temp['arr_0']
            meanmaxmkosm = temp['arr_1']
            deltmaxmkosm = temp['arr_1']
        else:
            if datatype == 'inpt' or mockfracperd != 0.:
                raise Exception('No file found for the null distribution of KS TS maxima...')
            else:
                print 'Fitting a polynomial to the null distribution of KS TS maxima...'
            binsmaxmkosm = linspace(amin(maxmkosm), amax(maxmkosm), 51)
            deltmaxmkosm = diff(binsmaxmkosm)
            meanmaxmkosm = (binsmaxmkosm[1:] + binsmaxmkosm[:-1]) / 2.
            cdfnmaxmkosm = cumsum(histogram(maxmkosm, bins=binsmaxmkosm)[0].astype(float)) / maxmkosm.size
            coefmaxmkosm = polyfit(meanmaxmkosm, cdfnmaxmkosm, 15)
            savez(path, coefmaxmkosm, meanmaxmkosm, deltmaxmkosm)
        cdfnmaxmkosm = poly1d(coefmaxmkosm)(meanmaxmkosm)
        # temp
        cdfnmaxmkosm[where(cdfnmaxmkosm < 0.)] = 0.
        pdfnmaxmkosm = diff(concatenate((array([0.]), cdfnmaxmkosm))) / deltmaxmkosm
        pval = interp(maxmkosm, cdfnmaxmkosm, meanmaxmkosm)
        stat = maxmkosm
    if stattype == 'maxmcnts':
        pval = sp.stats.gumbel_r.sf(stat, loc=log(numbfreq), scale=1.)
        
    # check if any p value is zero
    indxpvalzero = where(pval == 0.)[0]
    
    # calculate the significance
    if indxpvalzero.size > 0:
        plotpval = False
        sigm = zeros_like(pval) + 1000.
    else:
        plotpval = True
        sigm = sp.stats.norm.ppf(1. - pval)
    
    # plot the distribution of the TS and p value if the run is MCMC
    if numbiter > 1:

        if stattype == 'maxmpsec':
            figr, axis = plt.subplots()
            bins = linspace(amin(maxmpsec), amax(maxmpsec), 50)
            delt = diff(bins)[0]
            axis.hist(maxmpsec, bins=bins)
            axis.plot(bins, numbiter * delt * sp.stats.gumbel_r.pdf(bins, loc=log(numbfreq), scale=1.), label='Null (Gumbel)')
            axis.set_xlabel('max$(S)$')
            axis.set_ylabel('$N$')
            axis.legend()
            plt.tight_layout()
            path = pathimag + 'maxmpsec.pdf'
            plt.savefig(path)
            plt.close(figr)
            
        if stattype == 'maxmkosm':
            figr, axis = plt.subplots()
            bins = linspace(amin(maxmkosm), amax(maxmkosm), 50)
            delt = diff(bins)
            axis.hist(maxmkosm, bins=bins)
            axis.plot(meanmaxmkosm, numbiter * deltmaxmkosm * pdfnmaxmkosm, label='Null (Empirical)')
            axis.set_xlabel('max$(D_{KS})$')
            axis.set_ylabel('$N$')
            # temp
            axis.set_ylim([0., None])
            plt.tight_layout()
            path = pathimag + 'maxmkosm.pdf'
            plt.savefig(path)
            plt.close(figr)
    
        if plotpval:

            figr, axis = plt.subplots()
            minmpval = amin(pval)
            maxmpval = amax(pval)
            bins = logspace(log10(minmpval), log10(maxmpval), 50)
            axis.hist(pval, bins=bins)
            axis.set_xlim([minmpval, maxmpval])
            axis.set_xscale('log')
            axis.set_xlabel('$p$')
            axis.set_ylabel('$N$')
            plt.tight_layout()
            path = pathimag + 'pvalhist.pdf'
            plt.savefig(path)
            plt.close(figr)
    
            figr, axis = plt.subplots()
            axis.hist(sigm)
            axis.set_xlabel(r'Significance [$\sigma$]')
            axis.set_ylabel('$N$')
            plt.tight_layout()
            path = pathimag + 'sigmhist.pdf'
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

