#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments

genome=$1
mapping=$2

samplename=$(basename $mapping -spRep.sorted.bam)
genomeName=$(basename $genome .fa)
resultName=$genomeName-vs-$samplename.IS

fasta=/home/GLBRCORG/emcdaniel/Lake-MAGs/ref_MAGs_SAGs/98LQset/$genomeName.fa

# cd to mapping folder

cd /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/LQ98set_mapping

# profile command

inStrain profile $mapping $fasta -o /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/LQ98set_inStrain/$resultName -p 8 -s /home/GLBRCORG/emcdaniel/Lake-MAGs/ref_MAGs_SAGs/98LQset/all98LQset_scaffolds_to_bins.tsv