#! /bin/bash

# Run mapping of short metagenomic reads to a reference

# setup

mkdir metagenomes
mkdir refs
mkdir mappingResults

# Programs
tar -xvzf BBMap_38.07.tar.gz
tar -xvzf samtools.tar.gz
tar -xvzf python.tar.gz

# Python
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
chmod u+x *.py

# Copy over assembly and metagenomic read files from Gluster and decompress
# Read in variables for running mapping
ref=$1
meta=$2
outname=$3
refbase=$(basename $1)
metabase=$(basename $2)
refname=$(basename $refbase .fna)
metaname=$(basename $metabase .qced.fastq.tar.gz)
cp $1 refs/
cp $2 metagenomes/
cd metagenomes
tar -xzvf $metabase
mv metagenomes/*.fastq .
cd ../
metarun=$(basename $metabase .tar.gz)

# Perform mapping
bbmap/bbmap.sh ref=refs/$refbase in=metagenomes/$metarun outm=$outname idtag minid=0.95 nodisk -Xmx48g

# Make sorted BAM files
for file in mappingResults/*.bam; do
    outsort=$(basename $file .qced.fastq.tar.bam).sorted.bam;
    ./samtools/bin/samtools sort $file -o mappingResults/$outsort;
done

# get depth
for file in mappingResults/*.sorted.bam; do
    outdepth=$(basename $file .sorted.bam).depth;
    ./samtools/bin/samtools depth $file > mappingResults/$outdepth;
done

# sorted, indexed BAM file
for file in mappingResults/*.sorted.bam; do
    ./samtools/bin/samtools index $file;
done

# refernece lengths file
for file in refs/*.fna; do
    python countBases.py $file;
done
cat refs/*.len > refGenomes.len

# metagenomic reads file
for file in metagenomes/*.fastq; do
    awk '{s++}END{print FILENAME,s/4}' $file >> metaReads.txt;
done

# create stats file
for file in mappingResults/*.depth; do
    python calc-mapping-stats.py $file;
done 

# bring back stats and sorted/indexed BAM files to gluster
mkdir $refname-vs-$metaname
mv *.coverage.txt $refname-vs-$metaname/
mv mappingResults/*.sorted.bam $refname-vs-$metaname/
mv mappingResults/*.sorted.bam.bai $refname-vs-$metaname/
tar -czvf $refname-vs-$metaname.tar.gz $refname-vs-$metaname/
cp $refname-vs-$metaname.tar.gz /mnt/gluster/emcdaniel/GEODES/sparkling

rm *.txt
rm *.len
rm */*.fastq
rm */*.bam
rm */*.gz
rm *.gz