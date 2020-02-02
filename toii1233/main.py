import tesstoii.main
import numpy as np

strgtarg = 'TOI1233'
strgmast = '260647166'
labltarg = 'TOI 1233'

strgtoii = '1233'

massstar = 0.97 # [M_S]
radistar = 0.888 # [R_S]

tesstoii.main.main( \
                   #strgdata='qlop', \
                   strgtarg=strgtarg, \
                   labltarg=labltarg, \
                   
                   strgtoii=strgtoii, \
                   strgmast=strgmast, \
                   
                   massstar=massstar, \
                   radistar=radistar, \
                    
                  )

