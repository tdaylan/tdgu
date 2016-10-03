from __init__ import *

        
datapath = os.environ["PCAT_DATA_PATH"] + '/data/inpt/'
path = datapath + 'fermflux_cmp0_igal.fits'
maps = pf.getdata(path)[2, :, 3]
almcorig = hp.map2alm(maps)

fgl3lgal = deg2rad(5.)
fgl3bgal = deg2rad(10.)

tdpy.util.plot_maps('/Users/tansu/Desktop/maps.pdf', maps, satu=True, igal=True)

numbside = int(sqrt(maps.size / 12))
lgalcntr = deg2rad(array([0., 5., 5.]))
bgalcntr = deg2rad(array([5., 0., 5.]))
for k in range(lgalcntr.size):
    rttr = hp.rotator.Rotator(rot=[rad2deg(lgalcntr[k]), rad2deg(bgalcntr[k]), 0.], deg=True, eulertype='Y')
    fgl3bgal, fgl3lgal = rttr(pi / 2. - fgl3bgal, fgl3lgal)
    fgl3bgal = pi / 2. - fgl3bgal
    print 'lgalcntr'
    print rad2deg(lgalcntr[k])
    print 'bgalcntr'
    print rad2deg(bgalcntr[k])
    print 'fgl3lgal'
    print rad2deg(fgl3lgal)
    print 'fgl3bgal'
    print rad2deg(fgl3bgal)
    print
    
    almc = copy(almcorig)
    hp.rotate_alm(almc, lgalcntr[k], bgalcntr[k], 0.)
    maps = hp.alm2map(almc, numbside)
    tdpy.util.plot_maps('/Users/tansu/Desktop/maps%04d%04d.pdf' % (rad2deg(lgalcntr[k]), rad2deg(bgalcntr[k])), maps, satu=True, igal=True)

