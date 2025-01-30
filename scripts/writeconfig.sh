#!/bin/bash

frbname=$1
fluxcal=${2:-'0407'}
polcal=${3:-'vela'}
echo "FRB name: $frbname";
echo "Fluxcal name: $fluxcal";
echo "Polcal name: $polcal";

ml apptainer
set -o allexport #this is key to passing these to apptainer exec commands!

fcmpath=$(find /fred/oz313/data/frb${frbname}/ -type f -name 'fcm.*txt')
echo $fcmpath
#what if the file is fcm.txt.<numbers>? Then...
if [[ $fcmpath == "" ]] || [[ -z "$fcmpath" ]]; then
	fcmpath=$(find /fred/oz313/data/frb${frbname}/ -type f -name 'fcm.txt*')
fi
snoopypath=$(find /fred/oz313/data/frb${frbname}/frb/ -type f -name 'snoopyv2.cand')
fluxcaldatapath=$(find /fred/oz313/data/frb${frbname}/${fluxcal}/ -type d -name 'D1')
polcaldatapath=$(find /fred/oz313/data/frb${frbname}/${polcal}/ -type d -name 'D1')
frbdatapath=$(find /fred/oz313/data/frb${frbname}/frb/ -type d -name 'C000')
numant=$(find /fred/oz313/data/frb${frbname}/frb/ -type d -name 'ak*' | wc -l)
#args="$args -p $frbdatapath"

#echo $frbdatapath $args
#get central frequency via version of Adam's script
centfreq=$(apptainer exec -B /fred/oz313/:/fred/oz313/ /fred/oz313/sandbox_celebi_02mar24 bash -c 'source /opt/setup_proc_container && python3 /fred/oz313/processing/configs/getcentralfreq.py $frbdatapath')
#echo python getcentralfreq.py
echo "Central frequency is: $centfreq";

#get the first antenna
firstant=$(ls $frbdatapath | head -n 1 | sed -e 's/ak//g' | sed 's/^0*//')

cat >$frbname.config <<EOL
// initialse paramters with default values
includeConfig "default.config"

params.label = "$frbname" // FRB name
params.fcm = "$fcmpath" // FCM file
params.snoopy = "$snoopypath"  // Detection candidate file
params.container = "/fred/oz313/sandbox_celebi_02mar24" // celebi_Nov13.sif path to container/sandbox from which to run CELEBI
includeConfig "racs1low_v1.config" // assumes RACS-low v1

// Flux calibrator
includeConfig "$fluxcal.config"
params.data_fluxcal = "$fluxcaldatapath"


// Polarisation calibrator
includeConfig "$polcal.config"
params.centre_freq_polcal = $centfreq // Central frequency in MHz
params.data_polcal = "$polcaldatapath"


// FRB 
params.data_frb = "$frbdatapath"
params.ra_frb = "" // start with the multibeam localisation, it can be updated later
params.dec_frb = ""
params.dm_frb = 
params.centre_freq_frb = $centfreq // Central frequency in MHz
params.bw = 336

params.refant = $firstant
params.nants = $numant // number of antennas being included - UPDATE IF YOU MOVE ANTENNAS TO EXCLUDE THEM!


// Parameters to add along the way
// flag files
params.fluxflagfile = ""
params.fieldflagfile = ""
params.polflagfile = ""
EOL

cat test.config
