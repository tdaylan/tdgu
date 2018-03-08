from __init__ import *
import subprocess as subp

path = os.environ["TDGU_PATH"] + '/'
pathfileoutp = path + 'subm_auto_%s.log' % tdpy.util.retr_strgtimestmp()
os.system('rm %s' % pathfileoutp)
fileoutp = open(pathfileoutp, 'w')
               
listnamefile = []
listnamefunc = []
for namefile in os.listdir(path):
    if namefile.endswith(".py"):
        fileobjt = open(path + namefile, 'r')
        for line in fileobjt:
            if line.startswith('def pcat_'):
                
                namefunc = line[4:-1].split('(')[0]
                
                listnamefile.append(namefile)
                listnamefunc.append(namefunc)

listnamefile = array(listnamefile)
listnamefunc = array(listnamefunc)
numb = len(listnamefunc)
indx = choice(arange(numb), size=numb, replace=False)
listnamefile = listnamefile[indx]
listnamefunc = listnamefunc[indx]

for namefile, namefunc in zip(listnamefile, listnamefunc):
    print 'Processing confugiration %s, file %s...' % (namefunc, namefile)
        
    cmnd = 'python $TDGU_PATH/%s %s' % (namefile, namefunc)
    print cmnd
    try:
        pass
        subp.check_call(cmnd, shell=True)
        #os.system(cmnd)
    except Exception as excp:
        strg = str(excp)
        fileoutp.write('%s failed.\n' % namefunc)
        fileoutp.write(strg)
    print

fileoutp.close()

