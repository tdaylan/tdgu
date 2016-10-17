from __init__ import *

def testretr(retrtype):

    numbiter = 1000000

    if retrtype == 'objt':
        gdat = tdpy.util.gdatstrt()
        gdat.strg = 'hey'
    else:
        strg = 'hey'
    listtime = empty(numbiter)
    temp = empty(numbiter, dtype=object)
    
    for k in range(numbiter):
        
        timeinit = time.time()
        
        if retrtype == 'objt':
            temp[k] = gdat.strg
        else:
            temp[k] = strg
    
        timefinl = time.time()
        listtime[k] = timefinl - timeinit

    print 'retrtype'
    print retrtype
    print 'listtime'
    print sum(listtime)
    print

testretr('objt')
testretr('self')
