from __init__ import *

timeinit = time.time()
flop = 8e10
while True:
    thistime = time.time() - timeinit
    print '    %15d FLOPS' % int(flop * thistime)
    if thistime > 70:
        break

