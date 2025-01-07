#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: February 2023
# Description: download read files using accession numbers and NCBI SRA Toolkit on Gadi
# NOTE:	provide the list of sample_number accession (tab-separated, no header) as an argument:
#		Launch this script as qsub -v accessions=_______ download_ncbi.sh
######################


# Set PBS directives

#PBS -P nm31
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=32GB
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/nm31
#PBS -l wd
#PBS -v accessions


# set the container with the SRA Toolkit
container="/g/data/nm31/ben_anderson/singularity/SRA.sif"


# set up a log file
logfile="$(pwd)"/download_"${PBS_JOBID/\.gadi-pbs/}".log


# check accesions file is present
if [ -z "$accessions" ]; then
	echo -e "\nPlease provide an accessions file\n"
	exit 1
else
	accessions="$(readlink -f $accessions)"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting sequence download at $(date)\n**********" >> $logfile 2>&1


# load the singularity module
module load singularity >> $logfile 2>&1


# read the accessions file and download the sequences
while  IFS=$'\t' read -r -a myArray
do
sample_num="${myArray[0]}"
accession_num="${myArray[1]}"
singularity exec -H "$(pwd)" "$container" fasterq-dump "$accession_num" >> $logfile 2>&1
for file in "$accession_num"*; do mv "$file" "${file/${accession_num}/${sample_num}}"; done
for file in "$sample_num"*.fastq; do gzip "$file"; done
done < $accessions


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished downloading sequences at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
