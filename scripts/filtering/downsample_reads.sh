#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: Feb 2024
# Description: downsample Illumina paired-end reads on Gadi
#	Arguments are two read files and the number of read pairs to sample to  
#	If there are not that many reads in the files, they won't be modified  
######################


# Set PBS directive for if this is launched independently
#PBS -v read1,read2,numreads


# read and check args
if [ -n "$read1" ] && [ -n "$read2" ] && [ -n "$numreads" ]; then
	read1="$(readlink -f $read1)"
	read2="$(readlink -f $read2)"
	numreads="$numreads"
elif [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo -e "\nPlease specify input reads as arguments 1 and 2, and downsampling as argument 3\n"
	exit 1
else
	read1="$(readlink -f $1)"
	read2="$(readlink -f $2)"
	numreads="$3"
fi


# load module
module load seqtk/1.3


# determine how many reads are present in the files (assumes paired read files with the same number)
num1=$(echo $(zcat "$read1" | wc -l)/4 | bc)


# downsample, if needed
if [ $num1 -gt $numreads ]; then
	echo -e "Downsampling $read1 and $read2"
	seqtk sample -s1235 "$read1" "$numreads" | pigz > "${read1/\.fastq\.gz/_downsamp\.fastq\.gz}"
	seqtk sample -s1235 "$read2" "$numreads" | pigz > "${read2/\.fastq\.gz/_downsamp\.fastq\.gz}"
fi

echo "Finished $read1"
