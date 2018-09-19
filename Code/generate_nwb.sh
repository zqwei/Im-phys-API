#!/bin/bash

# Created by Ziqiang Wei on 2016-09-03.

# for i in $(ls -a $EIFILE/data_*.mat | cut -d "." -f 1); do
#     if [ ! -e $i.h5 ]; then
#         echo $i
#         python mat2h5.py $i.mat
#     fi
#     # mv $i.h5 GP43_highzoom_h5/$(basename $i).h5
# done

# mv $EIFILE/*.h5 GP517_highzoom_h5

for i in $(ls -a GP517_highzoom_h5/*.h5 | cut -d "." -f 1); do
    if [ ! -e $i.nwb ]; then
        echo $i
        python make_nwb_SingleCellEphysImaging.py $i.h5
    fi
done

mv GP517_highzoom_h5/*.nwb GP517_highzoom_nwb
