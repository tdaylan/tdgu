from __init__ import *
from scipy import fftpack

xaxi = linspace(0., 1., 100)
yaxi = linspace(0., 1., 100)
xaxi, yaxi = meshgrid(xaxi, yaxi, indexing='ij')

imag = zeros_like(xaxi)
for k in range(100):
    imag += sin(2. * pi * xaxi / rand()) * sin(2. * pi * yaxi / rand()) + 0.1 * rand(xaxi.size).reshape(sqrt(xaxi.size), sqrt(xaxi.size))

summgene(fft.fft2(imag))

psec = (abs(fft.fft2(imag))**2)[:50, :50]
psec = log10(psec)
# plot the power spectrum
plt.imshow(imag)
plt.savefig('/Users/tansu/Desktop/imag.pdf')
plt.imshow(psec)
plt.savefig('/Users/tansu/Desktop/psec.pdf')

