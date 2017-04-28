from __init__ import *

xdat = logspace(-3., -2., 20)
slop = 1.
pdfn = sp.stats.invgamma.pdf(xdat, slop, scale=0.5e-2)
plt.loglog(xdat, pdfn)
pdfn = sp.stats.invgamma.pdf(xdat, slop, scale=1e-2)
plt.loglog(xdat, pdfn)
pdfn = sp.stats.invgamma.pdf(xdat, slop, scale=2e-3)
plt.loglog(xdat, pdfn)
path = '/Users/tansu/Desktop/figr.pdf'
plt.savefig(path)
os.system('open ' + path)

