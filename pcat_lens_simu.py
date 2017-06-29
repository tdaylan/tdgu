from __init__ import *

# RA/DEC lists
liststrgrade = []
listrade = [[], []]

pathbase = os.environ["TDGU_DATA_PATH"] + '/pcat_lens_simu/'
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
numbside = 100
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
        
        print 'Coordinate number %d' % n
        
        # RA/DEC
        strgrade = liststrgrade[n]
        indxyaxi, indxxaxi = wcso.wcs_world2pix(listrade[0][n], listrade[1][n], 0)
        # check if the coordinate is inside the image
        if not isfinite(indxyaxi) or not isfinite(indxxaxi) or indxxaxi - numbsidehalf < 0 or indxyaxi - numbsidehalf < 0 or \
                                                                    indxxaxi + numbsidehalf > listdata[1].shape[1] or indxyaxi + numbsidehalf > listdata[1].shape[0]:
            continue
            #raise Exception('')

        path = pathdatapcat + 'lens%s%s%s%s_%04d.fits' % (liststrgrade[n][3:5], liststrgrade[n][6:8], liststrgrade[n][16:18], liststrgrade[n][19:21], numbside)
       
        print 'listdata[4]'
        print listdata[4].names
        if False:
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
       
        print 'MDRIZSKY'
        print listdata[4]['MDRIZSKY']
        print 'SKYSUB'
        print listdata[4]['SKYSUB']
        print

        if False:
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
        
        # cut out the image
        rate = listdata[1][indxxaxi-numbsidehalf:indxxaxi+numbsidehalf, indxyaxi-numbsidehalf:indxyaxi+numbsidehalf] # s^-1

        # gather different bands
        rate = rate[None, :, :, None]
        
        # find the number of photons per area per time per A per solid angle
        effa = 1. / listdata[4]['PHOTFLAM'][0] # erg^-1 cm^2 A
        timeobsv = listdata[4]['EXPTIME'][0] # s
        apix = (0.05 * pi / 3600. / 180.)**2 # sr^2
        expo = effa * timeobsv # erg^-1 cm^2 s A 
        flux = rate / effa / apix
        cnts = flux * expo * apix
        
        if False:
            print 'expo'
            print expo
            print 'rate'
            summgene(rate)
            print 'mean(cnts[:10, :10])'
            print mean(cnts[:10, :10])
            print 'cnts'
            summgene(cnts)
            print 'flux'
            summgene(flux)
        
        print 'Writing to %s...' % path
        pf.writeto(path, flux, clobber=True)
        
    
    #globals().get(sys.argv[1])()
        
