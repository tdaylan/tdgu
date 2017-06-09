from __init__ import *

pathbase = os.environ["TDGU_DATA_PATH"] + '/pcat_lens_simu/'
pathdata = pathbase + 'data/'
pathimag = pathbase + 'imag/'

pathslacpara = pathbase + 'data/slacpara.fits'
#tdpy.util.read_fits(pathslacpara, pathimag)

pathslacfull = pathbase + 'data/slacfull.fits'
#tdpy.util.read_fits(pathslacfull, pathimag)

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
        print 'data object has keys'
        print data.names

    arry = array(stack((head.keys(), head.values()), 1))
    listtype = []

    for n in range(arry.shape[0]):
        if arry[n, 0].startswith('TTYPE'):
            listtype.append(arry[n, 1])
    
    if len(listtype) != len(data[0]):
        raise Exception('Number of types does not match the number of fields.')
    

indx = where((data['Mph'] == 'E') & (data['Mul'] == 'S') & (data['Lens'] == 'A'))[0]

path = pathdata + 'list.txt'
fileobjt = open(path, 'w')
for item in data['SDSS'][indx]:
    fileobjt.write("%s %s %s %s %s %s \n" % (item[:2], item[2:4], item[4:9], item[9:12], item[12:14], item[14:]))

for k in range(len(indx)):
    print '%20s %20s' % (data['SDSS'][indx[k]], data['Name'][indx][k])

pathfile = pathdata + 'hst_10886_02_acs_wfc_f814w_drz.fits'
listdata = tdpy.util.read_fits(pathfile, full=True)

listhdun = astropy.io.fits.open(pathfile)
wcso = astropy.wcs.WCS(listhdun[2].header)


indxthis = where(data['SDSS'] == '002907.77-005550.5')[0]

print 'Working with ' + data['SDSS'][indxthis][0]
print data['_RA'][indxthis]
print data['_DE'][indxthis]
print

# temp 0 or 1 makes a difference!
indxyaxi, indxxaxi = wcso.wcs_world2pix(data['_RA'][indxthis], data['_DE'][indxthis], 0)

indxxaxi = int(indxxaxi[0])
indxyaxi = int(indxyaxi[0])
numbside = 100
numbsidehalf = numbside / 2

pathdatapcat = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
path = pathdatapcat + 'lens0029.fits'
print 'Writing to %s...' % path

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
print 'expo'
print expo
print 'rate'
summgene(rate)
flux = rate / effa / apix
cnts = flux * expo * apix
print 'mean(cnts[:10, :10])'
print mean(cnts[:10, :10])
print 'cnts'
summgene(cnts)
print 'flux'
summgene(flux)

pf.writeto(path, flux, clobber=True)

#globals().get(sys.argv[1])()

