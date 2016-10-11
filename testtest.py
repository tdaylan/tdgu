from __init__ import *

class WriteOnceDict(dict):
    def __setattr__(self, key, value):
        if hasattr(self, key):
            raise KeyError('{} has already been set'.format(key))
        super(WriteOnceDict, self).__setattr__(key, value)

y = tdpy.util.gdatstrt()
y.a = 1
y.b = 1
y.a = 1
print y.a, y.b 

