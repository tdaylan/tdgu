import datetime
from __init__ import *

# inputs
dail = -60
numbdays = 90

datenoww = datetime.datetime.now()
budg = zeros(numbdays)
indxdays = arange(numbdays)

listpanl = [ \
            [15, 2420], \
            [8, 1850], \
            [22, 1850], \
            [2, -2780], \
           ]

for k in indxdays:

    if k != 0:
        budg[k] = budg[k-1]
    
    budg[k] += dail

    dateprim = datenoww - datetime.timedelta(days=k)
    
    for panl in listpanl:
        if dateprim.day == panl[0]:
            budg[k] += panl[1]

    print '%s -- %d' % (dateprim.strftime('%Y-%m-%d'), budg[k])


figr, axis = plt.subplots(figsize=(16, 6))
axis.plot(indxdays, budg)
axis.set_xlabel('Days')
axis.set_ylabel('PnL [\$]')
axis.set_title('The Daylans Budget')
axis.axhline(0, color='r')
path = '/Users/tansu/Desktop/budget.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)
        
