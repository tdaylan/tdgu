# input PCAT speed per 100x100 pixel region
timeregi = 30. # [min]

# number of time frames in each region
numbtser = 13.7 * 4 * 24 * 60. / 30.

timeregitser = numbtser * timeregi / 60. / 24 # [day]
timeffim = 16.8e6 / 1e4 * timeregi # [day]
timesegm = 4. * timeffim / 7. # [week]
timefsky = 26 * timesegm / 7. # [week]

print 'Full frame, full sky: %d weeks per 1000 cores' % (timefsky / 1000.) 
