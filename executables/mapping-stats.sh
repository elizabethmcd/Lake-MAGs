#! /bin/bash

# Run mapping of short metagenomic reads to a reference

# Read in variables for running mapping
ref=$1
meta=$2
outname=$3

refbase=$(basename $1)
metabase=$(basename $2)

# Setup script for directories

mkdir metagenomes
mkdir refs
mkdir mappingResults

# Programs
tar -xvzf BBMap_38.07.tar.gz
tar -xvzf samtools.tar.gz

# Copy over assembly and metagenomic read files from Gluster and decompress

cp $2 metagenomes/
cd metagenomes
tar -xzvf $metabase
mv metagenomes/*.fastq .
cd ../
cp $1 refs/

metarun=$(basename $metabase .tar.gz)

# Perform mapping
bbmap/bbmap.sh ref=refs/$refbase in=metagenomes/$metarun outm=$outname idtag minid=0.95 nodisk -Xmx48g

# Make sorted BAM files
./samtools/bin/samtools sort $outname -o ${outname%.bam}.sorted.bam

# Move back only the sorted BAM files to Gluster for binning purposes
cp ${outname%.bam}.sorted.bam /mnt/gluster/emcdaniel/GEODES/mendota