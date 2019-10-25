#! /bin/bash 

for file in *.fastq; do
    name=$(basename $file .fastq);
    bash $file $name-F.fastq $name-R.fastq;
done