import pexo.main
import numpy as np

strgtarg = 'TOI1233'
ticitarg = 260647166
strgmast = 'TIC' + str(ticitarg)
labltarg = 'TOI 1233'

inclprio = np.array([90., 90., 90., 90.])
strgtoii = ['1233.01', '1233.02', '1233.03', '1233.04']

massstar = 0.97 # [M_S]
radistar = 0.888 # [R_S]

pexo.main.main( \
                   strgtarg=strgtarg, \
                   labltarg=labltarg, \
                   
                   strgtoii=strgtoii, \
                   strgmast=strgmast, \
                   ticitarg=ticitarg, \
                   
                   inclprio=inclprio, \

                   massstar=massstar, \
                   radistar=radistar, \
                    
                  )

