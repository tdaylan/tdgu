from __init__ import *

stdv = 2. * pi / 180.
norm = tdpy.util.retr_psfngausnorm(stdv)

lgal, bgal, numbpixl, apix = tdpy.util.retr_healgrid(256)
gridheal = array([lgal, bgal])
gridpnts = array([0., 0.])
anglpnts = angdist(gridheal, gridpnts, lonlat=True)
print norm
maps = norm * exp(-0.5 * (anglpnts / stdv)**2)

print sum(maps * apix)

