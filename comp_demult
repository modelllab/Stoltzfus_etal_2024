#!/bin/bash

###############################################################
# Demultiplex
###############################################################

# Use Ultraplex to demultiplex fastq files based on the barcode in the first 6 bases of R2 reads
# Ultraplex requires two input files:
  # 1) A fastq
  # 2) A csv containing sample and index information. Create a csv for every fastq file with a matching $base (here, the csv files are located in ../metadata/). Each line should contain a sample/index pair as follows: "index:sample"

# Run from directory containing zipped fastq sequences
if [ $1 == "demult" ]; then
	echo "DEMULTIPLEXING"
	for i in *.fastq.gz;
	do
	echo $i
	base=${i%%_*}
	echo $base
	echo ultraplex -i $i -b ../metadata/$base".csv" -d ../demult_fqs/$base
	ultraplex -i $i -b ../metadata/$base".csv" -d ../demult_fqs/$base;
	done
fi

# Next, moved all demultiplexed fastq.gzs to one directory

# After demultiplexing, extract just the sequence lines from fastq and write to new txt file
if [ $1 == "seq" ]; then
	echo "EXTRACTING SEQS"
	for i in *.fastq.gz;
	do
	echo $i
	base=${i%%.*}
	echo $base
	gunzip $i
	awk '{if(NR%4==2) print $0}' $base".fastq" > $base".fq";
	done
fi

# Remove "ultraplex_demux_" from the beginning of every file name
if [ $1 == "rename" ]; then
	echo "RENAMING"
	for i in *.fq; 
	do 
	mv "$i" "${i#ultraplex_demux_}"; 
	done
fi
