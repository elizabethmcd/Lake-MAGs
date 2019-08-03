#! /bin/bash 

# reformat r1 r2 fastqs to interleaved fastq files

# Read in variables for running mapping
r1=$(basename $1)
r2=$(basename $2)
out=$3

# Setup directories
mkdir metagenomes

# Programs
tar -xvzf BBMap_38.07.tar.gz

# copy over files
cp $1 metagenomes/
cp $2 metagenomes/
cd metagenomes/
gunzip $r1
gunzip $r2
cd ../

# run reformatting

bbmap/reformat.sh in1=metagenomes/$r1 in2=metagenomes/$r2 out=metagenomes/$outname

# bring back to gluster
cp metagenomes/$outname /mnt/gluster/emcdaniel/MENDOTA/

