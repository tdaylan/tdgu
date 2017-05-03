#!/bin/bash
dt=`date '+%d/%m/%Y %H:%M:%S'`
echo 'Making animations...'
echo "$dt"
source ~/.bashrc
python $PCAT_PATH/makeanim.py

