from __init__ import *
import xspec

def make_spec():

    timeinit = time.time()
    
    modl = xspec.Model("apec")
    #modl = xspec.Model("wabs*pow")
    #fs1 = xspec.FakeitSettings("response1.rsp", exposure = 1500.0)
    #fs1.background = "back1.pha"
    #xspec.AllData.fakeit(nSpectra=1, settings=fs1, applyStats=True, filePrefix="")
    xspec.AllData.fakeit(nSpectra=1, applyStats=True, filePrefix="", noWrite=True)
    xspec.AllData(1).dummyrsp(.3, 30., 100)

    summgene(array(xspec.AllData(1).noticed))
    summgene(array(xspec.AllData(1).values))
    
    timefinl = time.time()
    
    print (timefinl - timeinit) * 1e3



def pcat_chan_spec_mock(strgcnfgextnexec=None):
    
    dictargs = {}
    dictargs['exprtype'] = 'chan'
    dictargs['spatdisttype'] = ['line']
    dictargs['spectype'] = ['edis']
    dictargs['strgexpo'] = 1e4
    dictargs['elemtype'] = ['lghtline']
    dictargs['inittype'] = 'refr'
    # assume a pixel with side 100 arcsec
    anglfact = 3600. * 180. / pi
    dictargs['maxmgangdata'] = 100. / anglfact
    dictargs['numbsidecart'] = 1
    dictargs['anlytype'] = 'spec'
    
    # temp
    dictargs['numbelempop0reg0'] = 10
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
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )
    

globals().get(sys.argv[1])(*sys.argv[2:])

