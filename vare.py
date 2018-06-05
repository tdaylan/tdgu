import numpy as np, os#, random
import matplotlib.pyplot as plt

#from tdpy.util import summgene
from __init__ import *

import tensorflow as tf
from keras.models import Sequential, Model
from keras.layers import LSTM, RepeatVector

path = os.environ["TDGU_DATA_PATH"] + '/othr/imag/'

numbsequ = 10
aang = 3

x1 = linspace(-pi, pi)
y1 = sin(x1)
x2 = linspace(-pi + aang, pi + aang)
y2 = sin(x2)

def retr_datalabl(sizebtch):

    datatran = []
    batch_y = []
    for n in range(sizebtch):
        rand = random() * 2. * pi
        arry = linspace(rand, rand + 3. * pi, 2 * numbsequ)
        sig1 = sin(arry)
        sig2 = cos(arry)
        x1 = sig1[:numbsequ]
        y1 = sig1[numbsequ:]
        x2 = sig2[:numbsequ]
        y2 = sig2[numbsequ:]
        x_ = array([x1, x2]).T
        y_ = array([y1, y2]).T
        datatran.append(x_)
        batch_y.append(y_)

    datatran = array(datatran)
    batch_y = array(batch_y)

    return datatran, []

sizebtch = 100
datatran, _ = retr_datalabl(sizebtch)

print 'datatran'
summgene(datatran)

m = Sequential()
m.add(LSTM(2, input_shape=(10, 2)))
m.add(RepeatVector(10))
m.add(LSTM(2, return_sequences=True))

print 'summary'
print

m.summary()

m.compile(loss='mse', optimizer='adam')

history = m.fit(datatran, datatran, epochs=30, batch_size=sizebtch)

