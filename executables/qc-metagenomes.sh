#! /bin/bash

# Run QC of raw metagenomic reads

# Read in variables for running mapping
meta=$1
metabase=$(basename $1)

mkdir metagenomes

# Programs
tar -xvzf BBMap_38.07.tar.gz

# Copy over metagenomic read files from Gluster and decompress
cp $1 metagenomes/
cd metagenomes
gunzip $metabase
cd ..
metarun=$(basename $metabase .fastq.gz)

# Filter the reads
bbmap/bbduk.sh in=metagenomes/$metarun.fastq out=metagenomes/$metarun.qced.fastq qtrim=r trimq=10 maq=10

# Move back only the sorted BAM files to Gluster for binning purposes
tar -czvf metagenomes/$metarun.qced.fastq.tar.gz metagenomes/$metarun.qced.fastq
cp metagenomes/$metarun.qced.fastq.tar.gz /mnt/gluster/emcdaniel/GEODES/