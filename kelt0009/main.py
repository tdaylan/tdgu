import os
import tesstarg.util
pathbase = os.environ['TDGU_DATA_PATH'] + '/kelt0009/'
pathqlop = pathbase + 'tess/16740101.h5'
pathcsvv = pathbase + 'TESS.csv'
arry = tesstarg.util.read_qlop(pathqlop, pathcsvv=pathcsvv, stdvcons=1e-3)

