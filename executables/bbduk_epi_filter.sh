#! /bin/bash

###################
# Filtering metagenomic samples with BBduk part of the BBtools package
# For use on WEI GLBRC servers running HT Condor
# Elizabeth McDaniel 
##################

# set path where fastp is installed in local home directory bin
BBPATH=/opt/bifxapps/bbmap-38.32/
OUTDIR=/home/GLBRCORG/emcdaniel/Lake-MAGs/troutbog/epi/qced_reads
# queueing r1 r2 metagenomic reads and output folder/file names
reads=$1
sample=$(basename $reads .fastq)
$BBPATH/bbduk.sh in=$reads out=$OUTDIR/$sample.cleaned.fastq qtrim=r trimq=10 maq=10
