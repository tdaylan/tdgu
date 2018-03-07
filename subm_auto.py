from __init__ import *

path = os.environ["TDGU_PATH"] + '/'
pathfileoutp = path + 'subm_auto_%s.log' % tdpy.util.retr_strgtimestmp()
os.system('rm %s' % pathfileoutp)
fileoutp = open(pathfileoutp, 'w')
               
listnamefunc = []
for name in os.listdir(path):
    if name.endswith(".py"):
        print name
        fileobjt = open(path + name, 'r')
        for line in fileobjt:
            if line.startswith('def pcat_'):
                
                namefunc = line[4:-1].split('(')[0]
                
                listnamefunc.append(namefunc)

for namefunc in array(listnamefunc)[choice(arange(len(listnamefunc)), size=len(listnamefunc), replace=False)]:
    print 'Processing confugiration %s...' % namefunc
        
    cmnd = 'python $TDGU_PATH/%s %s' % (name, namefunc)
    print cmnd
    try:
        os.system(cmnd)
        if namefunc == 'pcat_lens_mock_trueminmdefs':
            raise Exception('')
    except Exception as excp:
        strg = str(excp)
        fileoutp.write('%s failed.\n' % namefunc)
        fileoutp.write(strg)
    print

fileoutp.close()

