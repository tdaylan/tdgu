import tesstoii.main
import numpy as np

strgtarg = 'TOI1233'
strgmast = '260647166'
labltarg = 'TOI 1233'

numbplan = 4
epocprio = np.array([2458571.335571, 2458586.566895, 2458572.398315, 2458572.111694])
periprio = np.array([14.175671, 19.593409, 6.203183, 3.795304])
duraprio = np.array([3.774857, 4.097541, 3.102675, 2.312993]) / 24. # [days]
tesstoii.main.main( \
                   #strgdata='qlop', \
                   strgtarg=strgtarg, \
                   labltarg=labltarg, \
                   strgmast=strgmast, \
                   numbplan=numbplan, \
                   epocprio=epocprio, \
                   periprio=periprio, \
                   duraprio=duraprio, \
                  )

