#! /bin/bash

# Append genome name to scaffolds file for each bin, combine into one file
for file in *.fa; do GENOME=`basename ${file%.fa}`; sed -i "s|^>|>${GENOME}~|" $file; done
cat *.fa > all-bins.fasta

# make list of scaffold-to-bins pairings
# for use with InStrain's genome_wide workflow to profile genome wide SNPs and covg
grep '>' bins.fasta | sed 's|[<>,]||g' | awk -F '~' '{print $1"~"$2"\t"$1}' > scaffolds-to-bins.tsv

# index the combined reference file
bowtie2-build bins.fasta bt2/combined-bins.fasta

# queue mapping 
# create file with list of complete paths of each metagenome files, second column the output destination
# If have just a few metagenomes, can write a script of a for loop such as 

for file in ~/EBPR/coassembly/metagenomes/cleaned_fastqs/*.qced.fastq; do
	name=$(basename $file .qced.fastq);
	bowtie2 -p 4 -x bt2/R1R2-EBPR-bins-index.fasta -q $file > mapping/$name-mapping.sam;
done

# or queue from a metadata file so can re-use the queuing script and not have to create new for loop scripts for each dataset


# queue inStrain profiling

# parse inStrain TSVs with python dictionary to create combined output table