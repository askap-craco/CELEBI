#!/bin/bash

final_position=$1
label=$2
startmjd=$3
centre_freq=$4
beamform_dir=$5

# parse final position file for RA and Dec
ra_line=$(head -2 $final_position | tail -1)
dec_line=$(tail -1 $final_position)

ra=$(echo $ra_line | awk -F ' ' '{print $2}' | sed 's/[a-z]//g')
dec=$(echo $dec_line | awk -F ' ' '{print $2}' | sed 's/[a-z]//g')

for npy in $(ls crops/*[IQUV]*); do
    # get tsamp from filename and convert to us
    tsamp=$(echo $npy | awk -F "_" '{print $(NF-1)}')
    tsamp_val=$(echo $tsamp | sed 's/[a-z]//g')
    tsamp_unit=$(echo $tsamp | sed 's/[0-9]//g')
    if [ "$unit" == "ms" ]; then
            tsamp_val=$((tsamp_val*1000))
    fi

    # get Stokes parameter from filename
    par=$(echo $npy | awk -F "_" '{print $NF}')

    outfile=$(echo $npy | sed 's/npy/fil/g' | sed 's/crops\///g')

    args="-s FRB$params.label"
    args="$args --tsamp $tsamp_val"
    args="$args --tstart $startmjd"
    args="$args --f0 $centre_freq"
    args="$args -r $ra"
    args="$args -d=$dec"
    args="$args -o $outfile"
    args="$args $npy"

    python3 $beamform_dir/npy2fil.py $args
done