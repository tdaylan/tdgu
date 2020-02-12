import numpy as np

# get ANTHEM list
objtfile = open('text.dat')
cntr = 0
liststrganth = []
for line in objtfile:
    if cntr % 2 == 1:
        linesplt = line.split(' ')
        liststrganth.extend(linesplt)
    cntr += 1

for k, strg in enumerate(liststrganth):
    if strg[-1] == '\n':
        liststrganth[k] = strg[:-1]

liststrganthtemp = liststrganth[:]
liststrganth = []
for strg in liststrganthtemp:
    if not strg in liststrganth:
        liststrganth.append(strg)

# get mass list
objtfile = open('mass.dat')
cntr = 0
liststrgmass = []
liststrgnumb = []
liststrgcomm = []
for line in objtfile:
    if cntr > 1:
        linesplt = line.split('\t')
        liststrgmass.append(linesplt[1])
        liststrgnumb.append(linesplt[4])
        liststrgcomm.append(linesplt[5][:-1])
    cntr += 1

for k in range(len(liststrgnumb)):
    if liststrgnumb[k] == '':
        liststrgnumb[k] = 0
    elif liststrgnumb[k] == '3+':
        liststrgnumb[k] = 3
    elif liststrgnumb[k] == '5+':
        liststrgnumb[k] = 5
    elif liststrgnumb[k] == '121 (+52 archive)':
        liststrgnumb[k] = 173
    elif liststrgnumb[k] == '47 PFS, 15 HARPS':
        liststrgnumb[k] = 72
    else:
        print('liststrgnumb[k]')
        print(liststrgnumb[k])
        liststrgnumb[k] = int(liststrgnumb[k])

liststrgmass = np.array(liststrgmass)
liststrgnumb = np.array(liststrgnumb)
liststrgcomm = np.array(liststrgcomm)
indxsort = np.argsort(liststrgnumb)[::-1]
print('liststrgmass')
print(liststrgmass)

liststrgmass = liststrgmass[indxsort]
liststrgnumb = liststrgnumb[indxsort]
liststrgcomm = liststrgcomm[indxsort]

for k, strgmass in enumerate(liststrgmass):
    for strg in liststrganth:
        if strg == strgmass:
            
            print('Match! %s, number of obs: %d, comment: %s' % (strg, liststrgnumb[k], liststrgcomm[k]))
