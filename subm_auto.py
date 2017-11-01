from __init__ import *

path = os.environ["TDGU_PATH"] + '/'
pathfileoutp = path + 'subm_auto.log'
os.system('rm %s' % pathfileoutp)
fileoutp = open(pathfileoutp, 'w')
               
# list of PCAT run plot outputs
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listrtag = fnmatch.filter(os.listdir(pathimag), '2*')

for name in os.listdir(path):
    if name.endswith(".py"):
        print name
        fileobjt = open(path + name, 'r')
        for line in fileobjt:
            if line.startswith('def pcat_'):
                
                namefunc = line[4:-1].split('(')[0]
                
                print 'Function %s' % namefunc
                    
                # check the available run outputs
                booltemp = False
                for rtag in listrtag:
                    strgstat = os.environ["PCAT_DATA_PATH"] + '/data/outp/' + rtag
                    if rtag.startswith(namefunc) and os.path.isfile(strgstat):
                        filestat = open(strgstat, 'r')
                        for line in filestat:
                            if line == 'gdatfinl written.\n':
                                booltemp = True
                if booltemp:
                    print 'Found a previously completed run.'
                    continue

                cmnd = 'python $TDGU_PATH/%s %s' % (name, namefunc)
                print cmnd
                try:
                    os.system(cmnd)
                    fileoutp.write('%s successfull.\n' % namefunc)
                except Exception as excp:
                    strg = str(excp)
                    fileoutp.write('%s failed.\n' % namefunc)
                    fileoutp.write(strg)
                print
                print
                print
                print
                print
                print
                print
                print

fileoutp.close()

