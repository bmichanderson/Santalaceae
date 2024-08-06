#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: July 2022
# Modified: April 2023 (simplified)
# Description: collate the output locus lengths from HybPiper2 output length summary files (args)
# Note: will create a "combined_lengths.tsv" file in the current directory
##########################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to collate HybPiper2 length stats')


# add arguments to parse
parser.add_argument('len_files', type=str, help='The length files', nargs='*')
parser.add_argument('-o', type = str, dest = 'out_file', help = 'The output file to create [default: combined_lengths.tsv]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
len_files = args.len_files
out_file = args.out_file

if not out_file:
	out_file = 'combined_lengths.tsv'


# read in each file and capture the data (sorted) into a list
master_list = []
for file in len_files:
	with open(file, 'r') as len_file:
		for lineno, line in enumerate(len_file):
			if lineno == 0:		# the first line is the loci
				loci = line.strip().split()[1: ]
			elif lineno == 2:	# the third line is the values
				pieces = line.strip().split()
				sample = pieces[0]
				lengths = pieces[1: ]

		touple_list = [(loci[i], lengths[i]) for i in range(len(loci))]
		touple_list.sort(key = lambda x: x[0])
		master_list.append((sample, touple_list))


# grab the sorted loci from the first sample (should be the same for all)
loci_list = [item[0] for item in master_list[0][1]]


# now write the output
with open(out_file, 'w') as output:
	output.write('Sample\t' + '\t'.join(loci_list) + '\n')
	for entry in master_list:
		output.write(str(entry[0]) + '\t' + '\t'.join([item[1] for item in entry[1]]) + '\n')
