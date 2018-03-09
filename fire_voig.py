from __init__ import *

def test_fire_mock(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'fire'
    dictargs['spectype'] = ['voig']
    dictargs['strgexpo'] = 1e3
    dictargs['spatdisttype'] = ['line']
    dictargs['elemtype'] = ['lghtlinevoig']
    #dictargs['inittype'] = 'refr'
    # assume a pixel with side 100 arcsec
    anglfact = 3600. * 180. / pi
    dictargs['makeplotinit'] = True
    dictargs['maxmgangdata'] = 100. / anglfact
    dictargs['numbsidecart'] = 1
    dictargs['anlytype'] = 'spec'
    
    # temp
    dictargs['numbelempop0reg0'] = 20
    dictargs['probspmr'] = 0.
    dictargs['numbswep'] = 100000
    dictargs['numbsamp'] = 1000
    
    # true < thrs < modl -- trad 
    # true < modl < thrs -- pcat
    
    # thrs < true < modl -- trad
    # modl < true < thrs -- pcat
    
    # thrs < modl < true -- trad
    # modl < thrs < true -- pcat
    listnamecnfgextn = ['nomi', 's2nrhigh']
    dictargsvari = {}
    for namecnfgextn in listnamecnfgextn:
        dictargsvari[namecnfgextn] = {}
    dictargsvari['s2nrhigh']['strgexpo'] = 1e5
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  listnamecnfgextn, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

globals().get(sys.argv[1])(*sys.argv[2:])

