from exoworlds.tess import extract_SPOC_data, extract_QLP_data
import os

pathbase = os.environ['TDGU_DATA_PATH'] + '/gj001252/'
pathtess = pathbase + 'tess/'

# extract SPOC data
extract_SPOC_data([pathtess+'tess2019169103026-s0013-0000000370133522-0146-s_lc.fits'], 
                  outdir=pathtess+'TESS_PDC/', PDC=True, auto_correct_dil=True, extract_centd=True, extract_dil=True)

    

