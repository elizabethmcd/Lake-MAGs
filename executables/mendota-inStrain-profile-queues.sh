#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments

genome=$1
mapping=$2

samplename=$(basename $mapping .mum-spRep.sorted.bam)
genomeName=$(basename $genome .fa)
resultName=$genomeName-vs-$samplename.IS

fasta=/home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mendota_finalBins/$genomeName.fa
genes=/home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mendota_finalBins/$genomeName.genes.fna

# cd to mapping folder

cd /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mappingResults

# profile command

inStrain profile --pairing_filter all_reads $mapping $fasta -o /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/inStrain/$resultName -p 8 -g $genes -s /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mendota_finalBins/mendota-scaffolds-to-bins.tsv