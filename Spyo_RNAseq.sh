#!/bin/bash

###############################################################
# Intro
###############################################################

# Differential expression analysis of WT, ∆trL, and ∆CRISPR-Cas(IIA) SF370 in BHI during logarithmic/stationary phase. Biological triplicate

# Libraries
	# Directional RNA libraries
	# NEBNext rRNA depletion kit and Ultra II Directional RNA Library Prep Kit for illumina

# Sequencing
	# 18 samples multiplexed on NovaSeq6000 SP lane
	# 2X50 reads

# Dependencies (versions used)
	# fastp (0.23.2)
	# bwa (0.7.17-r1188)
	# samtools (1.14)
	# deeptools (3.5.0)
	# salmon (1.8.0)
	# DESeq2 (1.20)

# SF370 fasta/GFF3 from: AE004092

# Run from command line as follows: $ Spyo_RNAseq.sh $1

###############################################################
# QC
###############################################################

# Write fastq quality summaries and trim by quality (quality score cutoff = 20, overrepresented sequence analysis enabled)
if [ $1 == "fastp" ]; then
	for i in *.fq.gz;
	do
	base=${i%%_R*}
	echo $base
	fastp -i $base"_R1.fq.gz" -o ../trimmed_fqs/$base"_trim_R1.fq.gz" -I $base"_R2.fq.gz" -O ../trimmed_fqs/$base"_trim_R2.fq.gz" -q 20 -c -h ../QC/$base.html;
	done
fi

###############################################################
# Alignments
###############################################################

# Index genome
# $ bwa index SF370_AE004092.fa

# Align reads to reference genome
if [ $1 == "align" ]; then
	for i in *R1.fq.gz;
	do
	base=${i%%_t*}
	echo $base
	bwa mem ../genomes/SF370_AE004092.fa $i $base"_trim_R2.fq.gz" > ../alignments/$base.sam;
	done
fi

# Convert sams to bams, sort, and index (sams and unsorted bams are deleted after they've served their purpose)
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

# Make BigWigs to view in igv. Normalized by counts per million, bin size 20bp
if [ $1 == "bigwig" ]; then
	for i in *.sorted.bam;
	do
	base=${i%%.*}
	echo $base
	bamCoverage --normalizeUsing CPM -bs 20 -b $i -o ../bigwigs/$base.bw; 
	done
fi


###############################################################
# Quantifying counts per gene (Salmon)
###############################################################

# Made transcript fasta file with gff3 downloaded from NCBI AE004092.2

# First, had to replace AE004092.2 with AE004092:
# $ sed -i.bakup  's/AE004092.2/AE004092/' SF370_AE004092.gff3

# Then, used gffread to make fasta:
# $ ../../tools/gffread-0.12.7.OSX_x86_64/gffread -w SF370_AE004092_transcripts.fa -g SF370_AE004092.fa SF370_AE004092.gff3
	
# Finally, made index for salmon:
# $ conda activate salmon
# $ salmon index -t SF370_AE004092_transcripts.fa -i SF370_AE004092_transcripts_index
 
# Quantified counts per gene using salmon (using gcBias flag DESeq2 manual recommends for RNASeq data)
if [ $1 == "salmon" ]; then
	for i in *R1.fq.gz;
	do
	base=${i%%_t*}
	echo $base
	salmon quant -i ../genomes/SF370_AE004092_transcripts_index/ -l A -1 $base"_trim_R1.fq.gz" -2 $base"_trim_R2.fq.gz" --validateMappings --gcBias -o ../salmon/$base;
	done
fi

###############################################################
# Importing counts with tximport
###############################################################

# Here, moved into RStudio. Commands entered are indicated with ">". Repeated with appropriate files for all pairwise comparisons; below shows WT vs. ∆tracr-L
# Following this nice tutorial: https://www.hadriengourle.com/tutorials/rna/
 
# Imported data with tximport:
# > setwd('/Users/mstoltz/220321_MS_Spyo_RNAseq')

# Loaded modules:
# > library(tximport)
# > library(GenomicFeatures)
# > library(readr)
 
# Linked transcript names to gene names:
# > txdb <- makeTxDbFromGFF("genomes/SF370_AE004092.gff3")
# > k <- keys(txdb, keytype = "GENEID")
# > tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

# Imported Salmon quantification
# made samples.txt table with two columns: header 'sample' (strain_log#) and 'conditions' (strain_log) 
# made separate files for just WT and ∆trL log and stationary to do pairwise comparisons (samples_log.txt and samples_stat.txt)
# Gene name column in quant.sf files all started with 'gene-' or 'rna-' which didn't match tx2gene names. To fix, removed 'gene-' and 'rna-' from quant.sf gene names:
if [ $1 == "fix_quants" ]; then
	for i in JW*
	do
	cd $i
	echo $i
	sed -i.bakup  's/gene-//' quant.sf
	sed -i.bakup  's/rna-//' quant.sf
	cd ..
	done
fi

# Then, continued with importing (using either samples_log.txt or samples_stat.txt first)
# > samples <- read.table("samples_log.txt", header = TRUE)
# > files <- file.path("salmon", samples$sample, "quant.sf")
# > names(files) <- paste0(samples$sample)
# > txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

###############################################################
# Differential expression analysis (DESeq2)
###############################################################

# Loaded library
# > library(DESeq2)

# Generated results table:
# > dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
# > dds <- DESeq(dds)
# > res <- results(dds)

# Summarized:
# > summary(res)

# Wrote summary into csv file
# > rescsv <- results(dds)
# > write.csv(rescsv, file="WTvUlt.csv")

# All resulting files in differential_expression/




