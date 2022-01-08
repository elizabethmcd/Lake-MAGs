#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''


# arguments
sam=$1
outfolder=$(basename $sam .mum-spRep.sorted.bam)-quick-profile

# cd to mapping results folder

cd /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mappingResults

# inStrain quick profile command

inStrain quick_profile -p 2 -s /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mendota_finalBins/mendota-scaffolds-to-bins.tsv -o /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/inStrain/quick_profiles/$outfolder $sam /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/mendota_finalBins/all-mendota-bins.fasta