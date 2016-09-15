from tdgu.__init__ import *
import tdpy.util

def corr_catl(lgalseco, bgalseco, fluxseco, lgalfrst, bgalfrst, fluxfrst, anglassc=deg2rad(1.)):

    numbfrst = lgalfrst.size

    indxsecoassc = zeros(numbfrst) - 1
    fluxassc = zeros(numbfrst)
    numbassc = zeros(numbfrst, dtype=int)
    distassc = zeros(numbfrst) + 1000.
    lgalbgalfrst = array([lgalfrst, bgalfrst])
    for k in range(lgalseco.size):
        lgalbgalseco = array([lgalseco[k], bgalseco[k]])
        dist = angdist(lgalbgalfrst, lgalbgalseco, lonlat=True)
        thisindxfrst = where(dist < anglassc)[0]
        
        if thisindxfrst.size > 0:
            
            # if there are multiple associated true PS, sort them
            indx = argsort(dist[thisindxfrst])
            dist = dist[thisindxfrst][indx]
            thisindxfrst = thisindxfrst[indx]
                
            # store the index of the model PS
            numbassc[thisindxfrst[0]] += 1
            if dist[0] < distassc[thisindxfrst[0]]:
                fluxassc[thisindxfrst[0]] = fluxseco[k]
                distassc[thisindxfrst[0]] = dist[0]
                indxsecoassc[thisindxfrst[0]] = k

    return indxsecoassc


#pathdata = '/n/fink2/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/fits/'
#pathimag = '/n/pan/www/tansu/imag/gaia_init/'

strgproc = os.uname()[1]
if strgproc == 'fink1':
    os.system('')
pathbase = os.environ["TDGU_DATA_PATH"]
pathimag = pathbase + '/imag/gaia_init/'
os.system('mkdir -p %s' % pathimag)

pathdata = pathbase + '/data/gaia_init/'
os.system('mkdir -p %s' % pathdata)

strg = 'GaiaSource_000-000-000.fits'

#tdpy.util.read_fits(pathdata + strg, pathimag=pathimag)

os.system('setenexport TESTAT_DIR=/Users/tansu/Documents/work/git/local/tdgu')
os.system('echo $TEST')
os.system('export PS1CAT_DIR=/n/fink2/dfink/decam-ucal-qz-time/chunks-qz-star-v3')
os.system('echo $PS1CAT_DIR')
numbpstr = 10
lgalpstr = 2. * (rand(numbpstr) - 0.5) * 180.
bgalpstr = 2. * (rand(numbpstr) - 0.5) * 90.
fluxpstr = rand(numbpstr)

lgalpstrneww = lgalpstr + randn(numbpstr) * 0.1
bgalpstrneww = bgalpstr + randn(numbpstr) * 0.1
fluxpstrneww = fluxpstr

lgalgaia = pf.getdata(pathdata + strg, 1)['l']
bgalgaia = pf.getdata(pathdata + strg, 1)['b']
fluxgaia = pf.getdata(pathdata + strg, 1)['phot_g_mean_mag']

indxpstr = corr_catl(lgalpstr, bgalpstr, fluxpstr, lgalgaia, bgalgaia, fluxgaia)

