#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''


# arguments
sam=$1
outfolder=$(basename $sam -spRep.sorted.bam)-quick-profile

# cd to mapping results folder

cd /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/LQ98set_mapping

# inStrain quick profile command

inStrain quick_profile -p 2 -s /home/GLBRCORG/emcdaniel/Lake-MAGs/ref_MAGs_SAGs/98LQset/all98LQset_scaffolds_to_bins.tsv -o /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/LQ98set_inStrain/quick_profiles/$outfolder $sam /home/GLBRCORG/emcdaniel/Lake-MAGs/ref_MAGs_SAGs/98LQset/all98LQset.fasta