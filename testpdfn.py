from __init__ import *
class rvarpars(sp.stats.rv_continuous):
    
    def _pdf(self, x, k):
        
        pdfn = x * 0. + 1.
        #pdfn =k * np.exp(-k*x)
        return pdfn

my_rv = rvarpars(name='exp', a=0.)
print my_rv.a, my_rv.b
print my_rv.numargs        # gets figured out automagically
print my_rv.cdf(4, k=3)
print my_rv.rvs(k=3, size=4)
print my_rv.expect(lambda x: 1, args=(2,))    # k=2 here

