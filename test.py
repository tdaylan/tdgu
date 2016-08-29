from __init__ import *

numbside = 256    
lgalheal, bgalheal, numbpixl, apix = tdpy.util.retr_healgrid(numbside)

maps = zeros(numbpixl)
maps = exp(-0.5 * (bgalheal / 5.**2) + 3. * exp(-0.5 * (sqrt(lgalheal**2 + bgalheal**2) / 10.)**2)
mapsmask = zeros(numbpixl)
mapsmask[where()] = 
mapsmask = hp.reorder(mapsmask, n2r=True)
indxpixlmask = where(mapsmask == 1)
tdpy.util.plot_maps('./mapsmask_%s.pdf' % rtag, mapsmask, numbsidelgal=1000, numbsidebgal=1000, satu=True)
    
    # plotting settings

