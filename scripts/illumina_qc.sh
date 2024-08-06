#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: April 2022
# Modified: Jun 2023
# Description: run QC on Illumina paired-end reads on Gadi
# NOTE: this script will be called from a separate job using nci-parallel
#	Arguments are two read files
#	Read pair files should be named with the same first number before an underscore
######################


# set a PBS directive in case this is run independently (so args can be specified)
#PBS -v read1,read2


# set max memory to specify (limit to 85% when passing to bbmap scripts)
# it should be 4 GB per cpu
if [ -z $cores_per ]; then
	# this was launched independently from the launch script
	cores_per="$PBS_NCPUS"
fi
max_mem=$(echo "${cores_per}*4*85/100" | bc)g
echo "This job is going to use $cores_per cores with max mem of $max_mem"


# define some variables
bbmap_loc="/g/data/nm31/bin/bbmap"

# bbduk args
min_len="50"		# minimum length of reads to keep after trimming
adapter_r1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
kmer="21"		# kmer value to use when searching for adapter in entire read
mismatch="2"		# how different an adapter sequence kmer can be for the kmer match
minq="20"		# minimum quality to trim to
kshort="8"
adapter_short1="${adapter_r1:0:15}"
adapter_short2="${adapter_r2:0:15}"
forceleft="9"	# trim all reads this many bases at the start of the read; if weirdness detected
dedup="yes"		# specify whether to run deduplification (yes or no)

# set possible reduction (0 for no reduction)
seqtk_num=0

# read and check args
if [ -n "$read1" ] && [ -n "$read2" ]; then
	read1="$(readlink -f $read1)"
	read2="$(readlink -f $read2)"
elif [ -z "$1" ] || [ -z "$2" ]; then
	echo -e "\nPlease specify input reads as arguments 1 and 2\n"
	exit 1
else
	read1="$(readlink -f $1)"
	read2="$(readlink -f $2)"
	if [ -n "$3" ]; then
		dedup="$3"
	fi
	if [ -n "$4" ]; then
		seqtk_num="$4"
	fi
fi


# create a directory for output and start a log file
# it will be based on the read file name
file_name=$(basename ${read1})
prefix=${file_name/_*/}
if [ ! -d "$prefix" ]; then
	mkdir -p "$prefix"
fi
cd "$prefix"
logfile="$prefix".log


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting Illumina QC at $(date)\n**********" >> $logfile 2>&1


# load modules
module load java/jdk-8.40 python3/3.10.0 >> $logfile 2>&1


## start QC
# run readlength stats on the raw files
echo -e "\n\n**********\nAssessing raw readlength distributions\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/readlength.sh in1="$read1" in2="$read2" out=readlength_raw.txt -Xmx"$max_mem" >> $logfile 2>&1


# remove duplicates (optical only; dupedist is for NovSeq)
if [ "$dedup" == "yes" ]; then
	echo -e "\n\n**********\nRemoving optical duplicates with clumpify\n**********\n" >> $logfile 2>&1
	"$bbmap_loc"/clumpify.sh in1="$read1" in2="$read2" out1=dedup_1.fastq.gz out2=dedup_2.fastq.gz \
		dedupe dupedist=12000 optical subs=0 -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
else
	cp "$read1" dedup_1.fastq.gz
	cp "$read2" dedup_2.fastq.gz
fi


# reduce dataset set size if requested
if [ "$seqtk_num" != "0" ]; then
	module load seqtk/1.3 >> $logfile 2>&1
	seqtk sample -s123 dedup_1.fastq.gz "$seqtk_num" | pigz > temp_1_sampled.fastq.gz
	seqtk sample -s123 dedup_2.fastq.gz "$seqtk_num" | pigz > temp_2_sampled.fastq.gz
	mv temp_1_sampled.fastq.gz dedup_1.fastq.gz
	mv temp_2_sampled.fastq.gz dedup_2.fastq.gz
fi


# trim adapters, for quality, and PhiX
echo -e "\n\n**********\nTrimming adapters, for quality, and PhiX\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in1=dedup_1.fastq.gz in2=dedup_2.fastq.gz out1=temp1_1.fastq.gz out2=temp1_2.fastq.gz ktrim=r \
	literal="$adapter_r1","$adapter_r2" k="$kmer" hdist="$mismatch" minlength="$min_len" mink=15 \
	tbo tpe forcetrimleft="$forceleft" -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in1=temp1_1.fastq.gz in2=temp1_2.fastq.gz out1=temp2_1.fastq.gz out2=temp2_2.fastq.gz ktrim=r \
	literal="$adapter_short1","$adapter_short2" k="$kshort" restrictright="$kmer" minlength="$min_len" \
	qtrim=r trimq="$minq" -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in1=temp2_1.fastq.gz in2=temp2_2.fastq.gz out1=trim_dedup_1.fastq.gz out2=trim_dedup_2.fastq.gz \
	ref=phix k=31 hdist="$mismatch" -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
rm temp* dedup*


# correct sequencing errors
echo -e "\n\n**********\nCorrecting sequencing errors\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/reformat.sh in1=trim_dedup_1.fastq.gz in2=trim_dedup_2.fastq.gz out=temp.fastq.gz -Xmx"$max_mem" \
	-t="$cores_per" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbmerge.sh in=temp.fastq.gz out=ecco.fastq.gz ecco mix vstrict ordered \
	prefilter=2 prealloc=t minlength="$min_len" -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
mv ecco.fastq.gz temp.fastq.gz
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/clumpify.sh in=temp.fastq.gz out=eccc.fastq.gz ecc passes=4 -Xmx"$max_mem" \
	-t="$cores_per" >> $logfile 2>&1
mv eccc.fastq.gz temp.fastq.gz
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/tadpole.sh in=temp.fastq.gz out=ecct.fastq.gz ecc k=62 ordered \
	prefilter=2 prealloc=t errormult1=64 -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
mv ecct.fastq.gz temp.fastq.gz
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in=temp.fastq.gz out1=correct_1.fastq.gz out2=correct_2.fastq.gz \
	minlength="$min_len" -Xmx"$max_mem" -t="$cores_per" >> $logfile 2>&1
rm temp* trim_dedup*


# rename the output
mv correct_1.fastq.gz "$prefix"_R1.fastq.gz
mv correct_2.fastq.gz "$prefix"_R2.fastq.gz


# run readlength stats on the final files
echo -e "\n\n**********\nAssessing final readlength distribution\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/readlength.sh in1="$prefix"_R1.fastq.gz in2="$prefix"_R2.fastq.gz out=readlength_qc.txt \
	-Xmx"$max_mem" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished Illumina QC at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
echo "This job finished"
