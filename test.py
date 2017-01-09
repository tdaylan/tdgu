from __init__ import *

data = rand(100000)

for k in range(5):
    t = time.time()
    for k in range(100000):
        main = rand() * data
    print time.time() - t
    
    t = time.time()
    for k in range(100000):
        main = empty(100000)
        main[:] = rand() * data
    print time.time() - t
    
