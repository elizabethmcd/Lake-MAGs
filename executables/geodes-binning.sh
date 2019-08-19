#! /bin/bash

# Binning with MetaBat

# Read in assembly variable, make a binning directory, and move over corresponding 
ref=$1
refbase=$(basename $1 .scaffolds.fasta)
mkdir $refbase

# Copy over assembly and BAM files from Gluster

cp $1 $refbase
cp /mnt/gluster/emcdaniel/GEODES/$refbase*.sorted.bam $refbase/
cd $refbase/
mkdir $refbase-bins

# Get depth 
jgi_summarize_bam_contig_depths --outputDepth $refbase-depth.txt *.bam

# Run metabat
metabat2 -i $refbase.scaffolds.fasta -a $refbase-depth.txt -o $refbase-bins/bin

# Zip up
tar -czvf $refbase-bins.tar.gz $refbase-bins/
mv $refbase-bins.tar.gz /mnt/gluster/emcdaniel/GEODES/.