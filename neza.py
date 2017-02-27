from numpy import *

def intgtdim(xdat, ydat, zdat):

    numb = xdat.shape[0]
    temp = empty(numb)
    for i in range(numb):
        temp[i] = trapz(zdat[i, :], xdat)

    intg = trapz(temp, ydat)

    return intg

xaxi = linspace(0., 1., 100)
yaxi = linspace(0., 1., 100)
xdat, ydat = meshgrid(xaxi, yaxi)
zdat = 0.5 * exp(-((xdat - 0.3)**2 + (ydat - 0.2)**2) / 10.)

intg = intgtdim(xaxi, yaxi, zdat)
print intg

