from __init__ import *

# list of PCAT run plot outputs
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listrtag = fnmatch.filter(os.listdir(pathimag), '2*')

for rtag in listrtag:
   
    # get configuration string
    liststrgtemp =  rtag[16:].split('_')[:-1]
    strgcnfg =  liststrgtemp[0]
    for strgtemp in liststrgtemp[1:]:
        strgcnfg += '_' + strgtemp
    
    listrtagprev = pcat.util.retr_listrtagprev(strgcnfg)
    
    print 'strgcnfg'
    print strgcnfg
    print 'listrtagprev'
    for rtagprev in listrtagprev:
        print rtagprev
    for k in range(len(listrtagprev)):
        cmnd = 'rm -rf %s' % listrtagprev[k]
        print cmnd
    cmndtotl = 'rm -rf'
    for k in range(len(listrtagprev)-1):
        cmndtotl += ' %s' % listrtagprev[k]
    #os.system(cmndtotl)
    print cmndtotl
    print
     
