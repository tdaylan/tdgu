from __init__ import *
import franlens

anglfact = 180. / pi
maxmgang = deg2rad(20.)
numbsidecart = 200
minmlgal = -maxmgang
maxmlgal = maxmgang
minmbgal = -maxmgang
maxmbgal = maxmgang
pathimag = tdpy.util.retr_path('tdgu', pathextnimag='strg_lens/', onlyimag=True)

# get the pixel grid for the lensing model
grid = franlens.PixelMap(maxmgang - 1.2 * maxmgang / numbsidecart, 2. * maxmgang / numbsidecart)
       
strgpara = ['lgallens', 'bgallens', 'beinlens']
numbiter = 5
maxmgangiter = 2. * maxmgang / 3.
listlgallens = linspace(-maxmgangiter, maxmgangiter, numbiter)
listbgallens = linspace(-maxmgangiter, maxmgangiter, numbiter)
listbeinlens = linspace(-maxmgang / 4., maxmgang / 4., numbiter)
for k in range(numbiter):
    print k
    for l in range(numbiter):
        for m in range(numbiter):

            # set the PSF scale
            psfnscal = 3
            
            # create the source object
            sourtype = 'Gaussian'
            lgalsour = maxmgang / 3.
            bgalsour = maxmgang / 3.
            fluxsour = 1e3
            sizesour = 0.02
            ratisour = 1.
            anglsour = 0.
            sourmodl = franlens.Source(sourtype, lgalsour, bgalsour, fluxsour, sizesour, ratisour, anglsour)
            
            # create the lens object
            lenstype = 'SIE'
            lgallens = listlgallens[k]
            bgallens = listbgallens[l]
            ellplens = 0.
            angllens = 0.
            sherlens = 0.15
            sanglens = -18.435
            beinlens = listbeinlens[m]
            lensmodl = franlens.LensModel(lenstype, lgallens, bgallens, ellplens, angllens, sherlens, sanglens, beinlens)
            
            # solve the lens equation
            mocktotlcnts = franlens.macro_only_image(grid, sourmodl, lensmodl, psfnscal).flatten()[None, :, None]
            
            path = pathimag + 'mocklensimag%d%d%d.pdf' % (k, l, m)
            titl = '%s=%.3g, %s=%.3g, %s=%.3g' % ( \
                                              strgpara[0], anglfact * listlgallens[k], \
                                              strgpara[1], anglfact * listbgallens[l], \
                                              strgpara[2], anglfact * listbeinlens[m], \
                                             )
            scat = [[anglfact * lgalsour, anglfact * bgalsour], [anglfact * lgallens, anglfact * bgallens]]
            tdpy.util.plot_maps(path, mocktotlcnts[0, :, 0], titl=titl, pixltype='cart', scat=scat, minmlgal=anglfact*minmlgal, maxmlgal=anglfact*maxmlgal, \
            #tdpy.util.plot_maps(path, log10(mocktotlcnts[0, :, 0]), titl=titl, pixltype='cart', scat=scat, minmlgal=anglfact*minmlgal, maxmlgal=anglfact*maxmlgal, \
                                                                                                           minmbgal=anglfact*minmbgal, maxmbgal=anglfact*maxmbgal)
        
            

