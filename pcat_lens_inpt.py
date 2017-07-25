from __init__ import *

def writ_data():

    # RA/DEC lists
    liststrgrade = []
    listrade = [[], []]
    
    pathbase = os.environ["TDGU_DATA_PATH"] + '/pcat_lens_inpt/'
    pathdata = pathbase + 'data/'
    pathimag = pathbase + 'imag/'
    
    # read SLACS tables
    print 'Reading SLACS tables...'
    pathslacpara = pathbase + 'data/slacpara.fits'
    pathslacfull = pathbase + 'data/slacfull.fits'
    hdun = pf.open(pathslacfull)
    numbhead = len(hdun)
    print '%s extensions found.' % numbhead
    for k in range(numbhead):
        print 'Extension %d' % k
        head = hdun[k].header
        data = hdun[k].data
        
        if data == None:
            print 'Data is None, skipping...'
            continue
        else:
            pass
            #print 'data object has keys'
            #print data.names
    
        arry = array(stack((head.keys(), head.values()), 1))
        listtype = []
    
        for n in range(arry.shape[0]):
            if arry[n, 0].startswith('TTYPE'):
                listtype.append(arry[n, 1])
        
        if len(listtype) != len(data[0]):
            raise Exception('Number of types does not match the number of fields.')
        
    # find the RA/DEC of relevant SLACS systems 
    indxgold = where((data['Mph'] == 'E') & (data['Mul'] == 'S') & (data['Lens'] == 'A'))[0]
    numbslac = indxgold.size
    
    path = pathdata + 'slacdownlist.txt'
    fileobjt = open(path, 'w')
    for k in indxgold:
        
        # construct the delimited RA/DEC string
        strgrade = '%s %s %s %s %s %s' % (data['SDSS'][k][:2], data['SDSS'][k][2:4], data['SDSS'][k][4:9], data['SDSS'][k][9:12], data['SDSS'][k][12:14], data['SDSS'][k][14:])
        
        ## fill the RA/DEC lists
        liststrgrade.append(strgrade)
        listrade[0].append(data['_RA'][k])
        listrade[1].append(data['_DE'][k])
        
        ## write the RA/DEC list of relevant SLACS systems to disc
        strgline = strgrade + ' \n'
        fileobjt.write(strgline)
    
    for k in range(len(indxgold)):
        print '%20s %20s %20g %20g' % (data['SDSS'][indxgold[k]], data['Name'][indxgold][k], data['_RA'][indxgold][k], data['_DE'][indxgold][k])
    
    # cutout properties
    numbside = 400
    numbsidehalf = numbside / 2
    
    # data path
    pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
    
    ## RA/DEC string of the reference star
    #strgradestar = '00 29 12.65 -00 53 59.7'
    strgradestar = '00 29 06.79 -00 54 07.5'
    liststrgrade.append(strgradestar)
    coorstar = ap.coordinates.SkyCoord(strgradestar, unit=(ap.units.hourangle, ap.units.deg))
    listrade[0].append(coorstar.ra.degree)
    listrade[1].append(coorstar.dec.degree)
    
    numbrade = len(listrade[0])
    print '%d coordinates found.' % numbrade
    
    # list of files to be read
    listnamefile = ['hst_10886_02_acs_wfc_f814w_drz.fits']
    numbfile = len(listnamefile)
    print '%d files found.' % numbfile
    
    for k, namefile in enumerate(listnamefile):
        
        print 'File number %d' % k
            
        # read the data fields
        pathfile = pathdata + namefile
        listdata = tdpy.util.read_fits(pathfile, verbtype=0)
        
        # read the WCS header
        listhdun = ap.io.fits.open(pathfile)
        wcso = ap.wcs.WCS(listhdun[2].header)
        
        # RA/DEC string
        strgrade = liststrgrade[k]
    
        # iterate over the RA/DEC list    
        for n in range(numbrade):
            
            # RA/DEC
            strgrade = liststrgrade[n]
            indxyaxi, indxxaxi = wcso.wcs_world2pix(listrade[0][n], listrade[1][n], 0)
            # check if the coordinate is inside the image
            if not isfinite(indxyaxi) or not isfinite(indxxaxi) or indxxaxi - numbsidehalf < 0 or indxyaxi - numbsidehalf < 0 or \
                                                                        indxxaxi + numbsidehalf > listdata[1].shape[1] or indxyaxi + numbsidehalf > listdata[1].shape[0]:
                continue
                #raise Exception('')
    
            path = pathdatapcat + 'lens%s%s%s%s_%04d.fits' % (liststrgrade[n][3:5], liststrgrade[n][6:8], liststrgrade[n][16:18], liststrgrade[n][19:21], numbside)
           
            if False:
                print 'listdata[4]'
                print listdata[4].names
                print 'strgrade'
                print strgrade
                print 'listrade[0][n]'
                print listrade[0][n]
                print 'listrade[1][n]'
                print listrade[1][n]
                print 'indxxaxi'
                print indxxaxi
                print 'indxyaxi'
                print indxyaxi
            
            indxxaxi = int(indxxaxi)
            indxyaxi = int(indxyaxi)
           
    
            if False:
                print 'MDRIZSKY'
                print listdata[4]['MDRIZSKY']
                print 'SKYSUB'
                print listdata[4]['SKYSUB']
                print 'indxxaxi'
                print indxxaxi
                print 'indxyaxi'
                print indxyaxi
                print 'listdata[1]'
                summgene(listdata[1])
                print 'EXPTIME'
                print listdata[4]['EXPTIME'][0]
                print 'PHOTFLAM'
                print listdata[4]['PHOTFLAM'][0]
                print 'CCDGAIN'
                print listdata[4]['CCDGAIN'][0]
                print
            
            # cut out the image
            rate = listdata[1][indxxaxi-numbsidehalf:indxxaxi+numbsidehalf, indxyaxi-numbsidehalf:indxyaxi+numbsidehalf] # s^-1
    
            # gather different bands
            rate = rate[None, :, :, None]
            
            # find the number of photons per area per time per A per solid angle
            effa = 1. / listdata[4]['PHOTFLAM'][0] # erg^-1 cm^2 A
            timeobsv = listdata[4]['EXPTIME'][0] # s
            apix = (0.05 * pi / 3600. / 180.)**2 # sr^2
            expo = effa * timeobsv # erg^-1 cm^2 s A 
            sbrt = rate / effa / apix
            cntp = sbrt * expo * apix
            
            if False:
                print 'expo'
                print expo
                print 'rate'
                summgene(rate)
                print 'mean(cntp[:10, :10])'
                print mean(cntp[:10, :10])
                print 'cntp'
                summgene(cntp)
                print 'sbrt'
                summgene(sbrt)
            if True:
                print 'sbrt'
                summgene(sbrt)
                print 'apix'
                print apix
            
            print 'Writing to %s...' % path
            pf.writeto(path, sbrt, clobber=True)
            
        

def pcat_lens_inpt():
    
    anglfact = 3600. * 180. / pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21
    
    # half-size of the image in pixels
    numbside = 400
    maxmgangdata = numbside * 0.5 * sizepixl

    # name of the data file
    strgexprsbrt = namedatasets + '_%04d.fits' % numbside
    
    if namedatasets == 'lens29075550':
        initbacpbac0ene0 = 1.1e-7
        fittmeanbacpbac0ene0 = 1.1e-7
        fittstdvbacpbac0ene0 = fittmeanbacpbac0ene0 * 1e-3
        fittscalbacpbac0ene0 = 'gaus'
    else:
        initbacpbac0ene0 = None
        fittmeanbacpbac0ene0 = None
        fittstdvbacpbac0ene0 = None
        fittscalbacpbac0ene0 = None

    mask = array([-0.3, 0.1, -0.1, 0.2]) / anglfact
    
    pcat.main.init( \
                   elemtype='lens', \
                   lensmodltype='none', \
                   numbswep=1000, \
                   numbswepplot=10000, \
                   #burntmpr=True, \
                   #initsizesour=1.5/anglfact, \
                   #initspecsourene0=1.5e-18, \
                   #mask=mask, \
                   indxenerincl=array([0]), \
                   initbacpbac0ene0=initbacpbac0ene0, \
                   fittmeanbacpbac0ene0=fittmeanbacpbac0ene0, \
                   fittstdvbacpbac0ene0=fittstdvbacpbac0ene0, \
                   fittscalbacpbac0ene0=fittscalbacpbac0ene0, \
                   serstype='intp', \
                   optihess=False, \
                   savestat=True, \
                   #inittype='reco', \
                   strgexpo=strgexpo, \
                   fittmaxmnumbpnts=array([0]), \
                   maxmgangdata=maxmgangdata, \
                   strgexprsbrt=strgexprsbrt, \
                  )
   

def pcat_lens_psfn():
    
    strgexpo = 7.37487548893e21
    anglfact = 3600. * 180. / pi
    maxmgangdata = 50. * 0.05 / anglfact
    numbiter = 1
    
    for k in range(numbiter):
        pcat.main.init( \
                       elemtype='lens', \
                       numbswep=50000, \
                       factthin=500, \
                       numbswepplot=10000, \
                       shrtfram=True, \
                       #mockonly=True, \
                       #makeplotintr=True, \
                       #burntmpr=True, \
                       optihess=False, \
                       #savestat=True, \
                       #inittype='reco', \
                       #makeplotinit=False, \
                       #makeplotfram=False, \
                       makeplotlpri=False, \
                       strgexpo=strgexpo, \
                       fittmaxmnumbpnts=array([0]), \
                       maxmgangdata=maxmgangdata, \
                       strgexprsbrt='lens29065407.fits', \
                      )
   

def pcat_lens_intrevalresicntp():

    anglfact = 3600. * 180. / pi
    sizepixl = 0.05 / anglfact
    
    # name of the dataset
    namedatasets = 'lens29075550'
    
    # exposure
    strgexpo = 7.37487548893e21

    # half-size of the image
    numbside = 400
    maxmgangdata = numbside * 0.5 * sizepixl
    
    # name of the data file
    strgexprsbrt = namedatasets + '_%04d.fits' % numbside
    
    if namedatasets == 'lens29075550':
        initbacpbac0ene0 = 1.1e-7
        fittmeanbacpbac0ene0 = 1.1e-7
        fittstdvbacpbac0ene0 = fittmeanbacpbac0ene0 * 1e-3
        fittscalbacpbac0ene0 = 'gaus'
    else:
        initbacpbac0ene0 = None
        fittmeanbacpbac0ene0 = None
        fittstdvbacpbac0ene0 = None
        fittscalbacpbac0ene0 = None

    pcat.main.init( \
                   elemtype='lens', \
                   makeplotinit=False, \
                   intrevalresicntp=True, \
                   strgexpo=strgexpo, \
                   initbacpbac0ene0=initbacpbac0ene0, \
                   fittmeanbacpbac0ene0=fittmeanbacpbac0ene0, \
                   fittstdvbacpbac0ene0=fittstdvbacpbac0ene0, \
                   fittscalbacpbac0ene0=fittscalbacpbac0ene0, \
                   inittype='reco', \
                   namerecostat='pcat_lens_inpt', \
                   fittmaxmnumbpnts=array([0]), \
                   maxmgangdata=maxmgangdata, \
                   strgexprsbrt=strgexprsbrt, \
                  )
    


globals().get(sys.argv[1])()
