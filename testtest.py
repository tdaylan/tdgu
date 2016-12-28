from __init__ import *
import numdifftools as nd

xdata = reshape(arange(0,1,0.1),(-1,1))
ydata = 1+2*exp(0.75*xdata)
fun = lambda c: (c[0]+c[1]*exp(c[2]*xdata) - ydata)**2
Jfun = nd.Jacobian(fun)
print allclose(abs(Jfun([1,2,0.75])), 0) # should be numerically zero
print 'Jfun'
print Jfun()
