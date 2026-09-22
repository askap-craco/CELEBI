#!/bin/bash
#-----------------------------------------------------------------------------
#	
#	This shell script is meant for running CELEBI ultimate with single pol support
#	
#	Initially written by Karol Desnos
#		last revised by AB, 8 Sep 26
#
#-----------------------------------------------------------------------------

# module load nextflow/23.10.1
module load java/17.0.4

#	Set the desired directories [Use abosulte paths to avoid confusion]
celebi_dir="/fred/oz313/src/users/abera/CELEBI_ultimate"
config_dir="/fred/oz313/processing/configs/ultconfigs"
log_dir="/fred/oz313/processing/logfiles"
report_dir="/fred/oz313/processing/collecting_reports"

#   Find the current user
current_user=$(id -u -n)
echo "Aloha, ${current_user}"
work_dir="/aphid/scratch-3month/${current_user}/work"

# Default config
rsm=""
unset -v mode
unset -v conf
#nf="/fred/oz313/nextflow_kdesnos/launch.sh"
nf="/fred/oz313/nextflow_kdesnos/build/releases/nextflow-25.03.1-edge-dist"

# FUnction to enable multiple tasks in the same -t field
taskinput() {
  local intask="$1"
  local -a words 
  
  if [[ "$intask" == "magical_leaf" ]]; then
    echo " --fcal --corrpcal --imgpcal --beampcal --calcpcal"
  elif [[ "$intask" == "leaf_storm" ]]; then
    echo " --fcal --corrfld --imgfld --corrpcal --imgpcal"
  elif [[ "$intask" == "charge_beam" ]]; then
    echo " --beamfrb --plotfrb"
  elif [[ "$intask" == "hyper_beam" ]]; then
    echo " --beampcal --beamfrb --calcpcal --plotfrb"  
  elif [[ "$intask" == "confusion" ]]; then
    echo " --mfimage --mfpos"
  elif [[ "$intask" == "psychic" ]]; then
    echo " --fcal --corrpcal --imgpcal --beampcal --calcpcal --corrfld --imgfld --beamfrb --plotfrb"
  else
    tasks=( $intask ) 
    for task in "${tasks[@]}"; do
        echo -n " --$task" 
    done
  fi
}

# Function to print help
print_help() {
  echo "Usage: $0 [options]"
  echo
  echo "Options:"
  echo "  -r                Enable resume mode."
  echo "  -t <task>         Specify the task to perform (required)."
  echo "  -c <config>       Specify the configuration name (required)."
  echo "  -p                Enable preview mode."
  echo "  -s                Enable stub run mode."
  echo "  -d <path>         Specify a work directory"
  echo "  -n <path>         Specify alternative nextflow path. Default is kdesnos' nextflow"
  echo "                    Use '-n nextflow' for native module."
  echo "  -h                Display this help message."
  echo
  echo "Example for resuming task localize on config 250313:"
  echo "  $0 -t localize -c 250313 -r"
}

OPTSTRING=":rt:c:pshn:d:"

while getopts ${OPTSTRING} opt; do
  case ${opt} in
    r)
      echo "Resume mode selected."
      rsm="-resume"
      ;;
    t)
      echo "${OPTARG} task selected."
      modstr="$(taskinput "${OPTARG}")"
      echo "celebi running with ${modstr}"
      mode="${modstr}"
      ;;
    c)
      echo "Config ${OPTARG} selected."
      conf="$OPTARG"
      ;;
    p) 
      echo "Preview mode selected."
      preview="-preview"
      ;;
    s)
      echo "Stub run."
      stub="-stub"
      profile="-profile stub_run"
      ;;
    d) 
	  echo "Used-defined work directory: ${OPTARG}"
	  work_dir="$OPTARG"
	  ;;
	n) 
	  echo "Used-defined nextflow location: ${OPTARG}"
	  nf="$OPTARG"
	  ;;
    h)
      print_help
      exit 0
      ;;
    ?)
      echo "Invalid option: -${OPTARG}."
      print_help
      exit 1
      ;;
  esac
done

if [ -z "$mode" ]; then
  echo "Missing -t option for selecting the task to perform." >&2
  print_help
  exit 1
fi

if [ -z "$conf" ]; then
  echo "Missing -c configuration name." >&2
  print_help
  exit 1
fi

# ======================================================================
# Linter Check 
# ======================================================================
echo "Running Pre-flight Configuration Linter..."
main_config="${config_dir}/$conf.config"

if [ ! -f "$main_config" ]; then
  echo "Error: Main config file $main_config not found!"
  exit 1
fi

# Build array of configurations to check, starting with the main one
config_list=("$main_config")

# Bash-fu: Find uncommented includeConfig lines and extract the path inside the quotes
included_configs=$(grep '^[[:space:]]*includeConfig' "$main_config" 2>/dev/null | sed -E "s/.*includeConfig[[:space:]]*['\"]([^'\"]+)['\"].*/\1/")

if [ -n "$included_configs" ]; then
  for inc in $included_configs; do
    echo " -> Found dynamically included config: $inc"
    config_list+=("$inc")
  done
fi

linter_output=$(python3 "${celebi_dir}/scripts/nf_pipeline_linter.py" \
  "${celebi_dir}/pipelines" \
  "${config_list[@]}" 2>&1)

echo "$linter_output"

# Check if the linter detected any collisions
if echo "$linter_output" | grep -q "\[COLLISION"; then
  echo "CRITICAL ERROR: Configuration collisions detected by linter! Abandoning run to prevent silent overwrites."
  exit 1
fi
echo "Linter passed. Proceeding with run..."
echo "======================================================================"

# Add extra info to logs
trace=nextflow.processor.TaskProcessor,nextflow.config.ConfigBuilder,nextflow.cli.CmdRun

# Check disk quota safely using exact kilobyte columns to avoid GiB/TiB unit string parsing errors
diskspace=$(quota -g oz313)

# 1. '/grp oz313/' waits until the correct group block starts before setting the 'found' flag to 1.
# 2. 'found && $1 == "/fred"' ensures we are in the right block AND exactly on the /fred row.
# 3. We calculate ($4 - $2) / 1073741824 to convert kbytes directly to TiB, then exit.
remaining=$(echo "$diskspace" | awk '/grp oz313/ {found=1} found && $1 == "/fred" { printf "%.2f", ($4 - $2) / 1073741824; exit }')

# Fallback in case awk returns an empty string, preventing `bc` from crashing
if [ -z "$remaining" ]; then
  echo "Warning: Could not parse remaining disk space. Assuming sufficient space."
  remaining=999.0
fi

echo "Remaining disk space = $remaining TB"

if (( $(echo "$remaining < 2.0" | bc -l) )); then
  echo "Disk space lower than 2.0 TB! ABANDON SHIP!!"
  exit 1
fi

# Check how recent the EOPs are
modate=$(stat -c %Y "/fred/oz313/auxfiles/.eops")
current_time=$(date +%s)
diff_seconds=$((current_time - modate))
diff_days=$((diff_seconds / 86400))
echo "EOPS are $diff_days days old"

if (( $(echo "$diff_days > 60" | bc -l) )); then
  echo "EOPs are older than 60 days! Please update EOPs in auxfiles/ using the following commands!!"
  echo "module load apptainer"
  echo "apptainer run /fred/oz313/celebi_container_18sep25.sif /usr/local/difx/bin/update_eop"
  echo "Update /fred/oz313/auxfiles/.eops"
  exit 1
fi

# Main CELEBI command
echo "Command run: $nf -trace $trace -c ${config_dir}/$conf.config run ${celebi_dir}/pipelines/main.nf $preview $rsm $stub -with-dag $mode $profile -w ${work_dir}/$conf >> ${log_dir}/frb$conf.out"

$nf -trace $trace -c ${config_dir}/$conf.config run ${celebi_dir}/pipelines/main.nf $preview $rsm $stub -with-dag $mode $profile -w ${work_dir}/$conf >> ${log_dir}/frb$conf.out

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: Nextflow command failed. Exiting."
  exit 1
fi

# Skip copying reports if the run was a preview, resume, or stub run
if [ -z "$rsm" ] && [ -z "$stub" ]; then
  # Copy the log with the same naming as other report files of the run
  trace_file=$(ls -t ${conf}_????-??-??_??_??_??_trace.txt 2>/dev/null | head -n 1)
  report_base="${trace_file%_trace.txt}"
  cp .nextflow.log ${report_dir}/${report_base}_log.log
  cp ${report_base}_report.html ${report_dir}/${report_base}_report.html
  cp ${report_base}_dag.html ${report_dir}/${report_base}_dag.html
  cp ${report_base}_trace.txt ${report_dir}/${report_base}_trace.txt
  cp ${report_base}_timeline.html ${report_dir}/${report_base}_timeline.html
else
  echo "Skipping report collection due to preview, resume, or stub run mode."
fi










