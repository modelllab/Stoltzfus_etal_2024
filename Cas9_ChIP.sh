#!/bin/bash

###############################################################
# Intro
###############################################################

# Used to identify Cas9 ChIP peaks in SF370, SF370 ∆tracr-L, 5448, and NZ131 logarithmic and stationary phase BHI cultures
# Sequencing: MiSeq, PE

# Fasta and GFF3 files from following genbank references:
	# 5448 = CP008776
	# NZ131 = CP000829
	# SF370 = AE004092
	
# Dependencies (version used)
	# cutadapt (3.5)
	# bwa (0.7.17-r1188)
	# samtools (1.16.1)
	# deeptools (3.5.0)
	# macs2 (2.2.9.1)

# Run from command line as follows: $ Cas9_ChIP.sh $1 

###############################################################
# Alignments
###############################################################

# Index genomes
if [ $1 == "index_fa" ]; then
	for i in *.fa;
	do
	bwa index $i
	done
fi

# Trim reads (Q score cutoff = 30)
# $ trim_galore -q 30 --paired --retain_unpaired <mate 1 fastq> <mate 2 fastq>

# Align reads to appropriate reference genome. Run from directory containing trimmed reads.
# Align 5448 samples
if [ $1 == "5448" ]; then
	for i in 5448*R1.fq;
	do
	base=${i%%_R*}
	bwa mem ../genomes/5448_CP008776.fa $base_R1.fq $base_R2.fq > ../aligned/$base.sam;
	done
fi
# Align NZ131 samples
if [ $1 == "NZ131" ]; then
	for i in NZ131*R1.fq;
	do
	base=${i%%_R*}
	bwa mem ../genomes/NZ131_CP000829.fa $base_R1.fq $base_R2.fq > ../aligned/$base.sam;
	done
fi
#Align SF370 and SF370 ∆tracr-L samples
if [ $1 == "SF370" ]; then
	for i in SF370*R1.fq;
	do
	base=${i%%_R*}
	bwa mem ../genomes/SF370_AE004092.fa $base_R1.fq $base_R2.fq > ../aligned/$base.sam;
	done 
fi

# Convert to sams to bams, sort, and index
if [ $1 == "samtools" ]; then
	for i in *.sam;
	do
	base=${i%%.*}
	echo $base
	samtools view -b -S $i > $base.bam && rm $i
	samtools flagstat $base.bam &> $base.flagstat.txt
	samtools sort $base.bam -o $base.sorted.bam && rm $base.bam
	samtools index $base.sorted.bam;
	done
fi


# Make BigWigs from bams for visualization on IGV. Normalized by counts per million in 20 bp bins.
if [ $1 == "bigwig_CPM" ]; then
	for i in *.sorted.bam;
	do
	base=${i%%.*}
	echo $base
	bamCoverage --normalizeUsing CPM -bs 20 -b $i -o ../bigwigs/$base.CPM.bw; 
	done
fi

###############################################################
# Call Peaks
###############################################################

# Call peaks with MACS2 (using 90% genome size for mappable genome)
if [ $1 == "macs2" ]; then
	macs2 callpeak -t 5448_log_ChIP.bam -c 5448_log_input.bam -g 1.667e6 -n 5448_log --outdir ../peaks -f BAMPE
	macs2 callpeak -t 5448_stat_ChIP.bam -c 5448_stat_input.bam -g 1.667e6 -n 5448_stat --outdir ../peaks -f BAMPE
	macs2 callpeak -t NZ131_log_ChIP.bam -c NZ131_log_input.bam -g 1.667e6 -n NZ131_log --outdir ../peaks -f BAMPE
	macs2 callpeak -t NZ131_stat_ChIP.bam -c NZ131_stat_input.bam -g 1.667e6 -n NZ131_stat --outdir ../peaks -f BAMPE
	macs2 callpeak -t SF370_log_ChIP.bam -c SF370_log_input.bam -g 1.667e6 -n SF370_log --outdir ../peaks -f BAMPE
	macs2 callpeak -t SF370_stat_ChIP.bam -c SF370_stat_input.bam -g 1.667e6 -n SF370_stat --outdir ../peaks -f BAMPE
	macs2 callpeak -t SF370dtrL_log_ChIP.bam -c SF370dtrL_log_input.bam -g 1.667e6 -n SF370dtrL_log --outdir ../peaks -f BAMPE
	macs2 callpeak -t SF370dtrL_stat_ChIP.bam -c SF370dtrL_stat_input.bam -g 1.667e6 -n SF370dtrL_stat --outdir ../peaks -f BAMPE
fi