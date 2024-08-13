#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: April 2022
# Modified: Mar 2023 (new API format)
# Description: run a download script from Bioplatforms on Gadi
# NOTE: this expects there is a `download.sh` script in the current directory
#		along with a `tmp` folder containing a list of URLs and md5sums
#		These can be obtained when selecting bulk download from Bioplatforms
#		Launch this script as qsub -v api_key=_______ download_seqs.sh
######################

# Set PBS directives

#PBS -P nm31
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l walltime=08:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v api_key


# set up a log file
logfile=download_"${PBS_JOBID/\.gadi-pbs/}".log


# check key is present, and export it
if [ -z "$api_key" ]; then
	echo -e "\nPlease provide the API key\n"
	exit 1
else
	export CKAN_API_TOKEN="$api_key"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting sequence download at $(date)\n**********" >> $logfile 2>&1


# launch the script
# the `download.sh` script should automatically download
./download.sh >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished downloading sequences at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
