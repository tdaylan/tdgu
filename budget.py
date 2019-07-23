import datetime
from __init__ import *

# inputs
numbdays = 365

if strgmode == 'budg':
    dail = -100
    numbacct = 2
else:
    dail = -1800 / 3500.

datenoww = datetime.datetime.now()
budg = zeros((numbacct, numbdays))
indxdays = arange(numbdays)

if strgmode == 'budg':
    # BoA checking
    budg[0, 0] += 10085
    
    # BoA credit card
    budg[0, 0] += 0
    
    # Capital One checking
    budg[1, 0] += 2465
    
    # Capital One credit card
    budg[1, 0] += 0
else:
    budg[1, 0] += 92

listpanlonet = [ \
                #[0, 2018, 9, 1, -3500], \
               ]

if strgmode == 'budg':
    listpanl = [ \
                [0, 15,              2620], \
                [0, 30,             2620], \
                [1, 12,              1907], \
                [1, 27,             1907], \
                [1, 1,             -2780], \
               ]
else:
    listpanl = [ \
               ]

for k in indxdays:

    if k != 0:
        budg[:, k] = budg[:, k-1]
    
    budg[0, k] += dail

    dateprim = datenoww + datetime.timedelta(days=k)

    for panl in listpanlonet:
        if dateprim.day == panl[3] and dateprim.month == panl[2] and dateprim.year == panl[1]:
            budg[panl[0], k] += panl[4]
    for panl in listpanl:
        if dateprim.day == panl[1]:
            budg[panl[0], k] += panl[2]
    #print '%s -- %d' % (dateprim.strftime('%Y-%m-%d'), budg[k])

figr, axis = plt.subplots(figsize=(16, 6))
axis.plot(indxdays, sum(budg, axis=0), color='b')
axis.plot(indxdays, budg[0, :], color='b', alpha=0.5, ls='--')
axis.plot(indxdays, budg[1, :], color='b', alpha=0.5, ls='-.')

datetick = axis.get_xticks()
numbtick = len(datetick)
strgindxdays = zeros(numbtick, dtype=object)
for k, date in enumerate(datetick):
    dateprim = datenoww + datetime.timedelta(days=date)
    strgindxdays[k] = dateprim.strftime("%y-%m-%d")
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
taxx = 7.3e4 * 0.15
axis.axhline(taxx, color='g', ls='--', alpha=0.5)
axis.axhline(taxx + 1e4, color='black', ls='--', alpha=0.5)
axis.axhline(taxx + 2e4, color='black', ls='--', alpha=0.5)
path = '/Users/tansu/Desktop/budget.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)
        
