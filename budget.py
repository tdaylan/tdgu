import datetime
from __init__ import *

# inputs
dail = -60
numbdays = 90

datenoww = datetime.datetime.now()
budg = zeros(numbdays)
indxdays = arange(numbdays)

listpanl = [ \
            [7,              2150], \
            [21,             2150], \
            [8,              1850], \
            [22,             1850], \
            [2,             -2780], \
           ]

for k in indxdays:

    if k != 0:
        budg[k] = budg[k-1]
    
    budg[k] += dail

    dateprim = datenoww + datetime.timedelta(days=k)

    for panl in listpanl:
        if isinstance(panl[0], int):
            if dateprim.day == panl[0]:
                budg[k] += panl[1]
        else:
            if dateprim.day == int(panl[0][8:10]) and dateprim.month == int(panl[0][5:7]) and dateprim.year == int(panl[0][:4]):
                budg[k] += panl[1]
    #print '%s -- %d' % (dateprim.strftime('%Y-%m-%d'), budg[k])

figr, axis = plt.subplots(figsize=(16, 6))
axis.plot(indxdays, budg)

datetick = axis.get_xticks()
numbtick = len(datetick)
strgindxdays = zeros(numbtick, dtype=object)
for k, date in enumerate(datetick):
    dateprim = datenoww + datetime.timedelta(days=date)
    strgindxdays[k] = dateprim.strftime("%m-%d")
    #print 'k'
    #print k
    #print 'date'
    #print date
    #print 'dateprim'
    #print dateprim
    #print 'strgindxdays[k]'
    #print strgindxdays[k]
    #print
    
axis.set_xticklabels(strgindxdays)
axis.set_xlabel('Days')
axis.set_ylabel('PnL [\$]')
axis.set_title('The Daylans Budget')
axis.axhline(0, color='r')
path = '/Users/tansu/Desktop/budget.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)
        
