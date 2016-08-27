from __init__ import *

def icdf(unit, slop):

    varb = (unit * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))
    
    return varb

def pdfn(varb, slop):
  
    norm = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop))
    pdfn = norm * varb**(-slop)
    
    return pdfn


minm = 1e-9
maxm = 10e-9
numbvarb = 90
numbsamp = 100000
numbbins = 40
alph = 0.5

binssing = linspace(minm, maxm, numbvarb + 1)
meansing = (binssing[:-1] + binssing[1:]) / 2.
deltsing = binssing[1:] - binssing[:-1]

binsdoub = linspace(2. * minm, 2. * maxm, 2 * numbvarb)
meandoub = (binsdoub[:-1] + binsdoub[1:]) / 2.
deltdoub = binsdoub[1:] - binsdoub[:-1]

bins = linspace(minm, 2. * maxm, 2 * numbvarb + 1)

arry = empty((2, numbsamp))

minmslop = 1.5
maxmslop = 3.
numbslop = 4
sloparry = linspace(minmslop, maxmslop, numbslop)
for n in range(numbslop):
    slop = sloparry[n]
    for k in range(2):
        arry[k, :] = icdf(rand(numbsamp), slop)
    
    totl = sum(arry, 0)
    
    powrprob = pdfn(meansing, slop)
    
    convprob = convolve(powrprob, powrprob) * deltdoub[0]
    
    indxdoub = where(meandoub <= maxm)[0]
    convprobpoly = polyval(polyfit(meandoub[indxdoub], convprob[indxdoub], 8), meandoub[indxdoub])
    
    figr, axis = plt.subplots()
    axis.hist(arry[k, :], bins=bins, alpha=alph, label='$f_1$ (Sampled)', color='b')
    axis.hist(totl, bins=bins, alpha=alph, label='$f_0$ (Sampled)', color='g')
    axis.plot(meansing, powrprob * numbsamp * deltsing, label='$f_1$ (Analytic)', color='b')
    axis.plot(meandoub, convprob * numbsamp * deltdoub[0], label='$f_0$ (Numerically convolved)', color='g')
    
    axis.plot(meandoub[indxdoub], convprobpoly * numbsamp * deltdoub[indxdoub], label='$f_0$ (Fit)', color='r')

    axis.set_ylim([0.5, numbsamp])
    axis.set_xlabel('$f$')
    axis.set_xlim([amin(bins), amax(bins)])
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylabel('$N_{samp}$')
    axis.legend()
    plt.tight_layout()
    pathfold = os.environ["TDGU_DATA_PATH"] + '/imag/powrpdfn/'
    os.system('mkdir -p ' + pathfold)
    figr.savefig(pathfold + 'powrpdfn%04d.pdf' % n)
    plt.close(figr)
    
