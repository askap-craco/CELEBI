#!/bin/bash
###############################################################################
# craftpy2 stage 1: get info
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-19
#
# Collates information (e.g. filenames) important for running the rest of the
# pipeline into a single file for convenience.
###############################################################################

prefix=[`basename -- $0`]

if [ $# -lt 8 ]||[ $# -gt 9 ]; then
  echo "$prefix ERROR: incorrect number of arguments"
  echo "$prefix ERROR: usage: $0 label data_frb data_cal fcm ra_frb dec_frb caL_file cal_name [outfile]"
  exit
fi

label=$1
mkdir label

if [ "$9" == "" ]; then
  outfile=${label}.info
else
  outfile=$9
fi

# FRB data
data_frb=$2

# calibrator data
data_cal=$3

# antenna list - all antennas in *both* data_frb and data_cal
# if the list of antennas in the frb and calibrator datasets are different, need
# to manually change the directory names to exclude the inconsistent antennas
antlist_frb=`(cd ${data_frb} && ls -d ak*)`
antlist_frb=`echo $antlist_frb | sed -r 's/ /,/g'` # make list comma separated

antlist_cal=`(cd ${data_cal} && ls -d ak*)`
antlist_cal=`echo $antlist_cal | sed -r 's/ /,/g'` # make list comma seperated

if [ "$antlist_frb" != "$antlist_cal" ]; then
  echo "$prefix ERROR: FRB and calibrator antenna lists are inconsistent!"
  echo "$prefix ERROR: Please rename the ak* directories to make them consistent."
  echo "$prefix ERROR: antlist_frb=$antlist_frb"
  echo "$prefix ERROR: antlist_cal=$antlist_cal"
  exit
fi

antlist=$antlist_frb

# snoopy file
echo "$prefix Please make sure you've put the Bayesian DM in the snoopy file!"
snoopy=${data_frb}/snoopyv2.cand

# primary detection beam + position
# extract the primary beam from the snoopy file
snoopy_entry=`tail -1 $snoopy`
beam_prim=`cut -d' ' -f7 <<< $snoopy_entry`
# in the ak* directories, the X and Y beams are numbered by 2*prim_beam and
# 2*prim_beam + 1 respectively
beam_x=$(( $beam_prim * 2 ))
beam_x0=`printf "%02d" $beam_x`
# get the RA and Dec of the beam from the vcraft header file of an antenna
ant=`echo $antlist | cut -c1-4` # grab first antenna from antlist
hdr=${data_frb}/${ant}/beam${beam_x0}/${ant}_c1_f0.vcraft.hdr
ra_beam=$(cut -d ' ' -f2 <<< `grep BEAM_RA $hdr`)
dec_beam=$(cut -d ' ' -f2 <<< `grep BEAM_DEC $hdr`)

# fcm
# downloaded from aktos11 (may have to just do manually)
fcm=$4

# Bayesian FRB position
ra_frb=$5
dec_frb=$6

# Calibrator position
cal_file=$7
cal_name=$8
cal_line=`grep $cal_name $cal_file`
ra_cal=`echo $cal_line | cut -d " " -f 2`
dec_cal=`echo $cal_line | cut -d " " -f 3`

echo "label=$label" > $outfile
echo "data_frb=$data_frb" >> $outfile
echo "data_cal=$data_cal" >> $outfile
echo "antlist=$antlist" >> $outfile
echo "snoopy=$snoopy" >> $outfile
echo "beam_prim=$beam_prim" >> $outfile
echo "ra_beam=$ra_beam" >> $outfile
echo "dec_beam=$dec_beam" >> $outfile
echo "fcm=$fcm" >> $outfile
echo "ra_frb=$ra_frb" >> $outfile
echo "dec_frb=$dec_frb" >> $outfile
echo "ra_cal=$ra_cal" >> $outfile
echo "dec_cal=$dec_cal" >> $outfile

echo "$prefix Created $outfile"
cat $outfile
