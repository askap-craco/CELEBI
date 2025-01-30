#!/bin/bash

# set up enviroments and variables
ml nextflow/23.04.2

celebi_dir="/fred/oz313/src/users/tyson/CELEBI"
config_dir="/fred/oz313/processing/configs"
work_dir="/fred/oz313/processing/work"


# run nextflow
celebi_command="nextflow ${celebi_dir}/pipelines/main.nf"

# add processes
for proc in ${@:2}; do
    celebi_command+=" ${proc}"
done

# add config file
celebi_command+=" -c ${config_dir}/${1}.config"

# add work directory
celebi_command+=" -w ${work_dir}/${1}"

echo $celebi_command
echo ""


# print infomation of celebi run
date_str=$(date '+%d-%m-%Y %H:%M:%S')

echo "Running [${1}] through CELEBI [${date_str}]"
echo "==========================================="
echo "Processes:"
echo "=========="

for proc in ${@:2}; do 
    echo $proc
done

echo "=========="
echo ""
echo "config file: ${config_dir}/${1}"
echo "work directory: ${work_dir}/${1}"
echo "celebi directory: ${celebi_dir}"
echo ""


# run command
eval $celebi_command > "frb${1}.out"