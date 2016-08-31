from __init__ import *

numbside = 2**3
numbpixl = numbside**2 * 12

mapsstor = rand(numbpixl)
maps = copy(mapsstor)
maps = hp.ma(maps)
maps.mask = zeros(numbpixl, dtype=bool)
maps.mask[0:4] = True
print hp.anafast(maps)

maps = copy(mapsstor)
maps[0:4] = 0.
print hp.anafast(maps)

