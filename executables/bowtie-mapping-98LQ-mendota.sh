#! /bin/bash

# setup environment for samtools dependencies to work 

export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

# arguments
reads=$1
samplename=$(basename $reads .mum.QCed.fastq)

cd /home/GLBRCORG/emcdaniel/Lake-MAGs/mendota/LQ98set_mapping

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x /home/GLBRCORG/emcdaniel/Lake-MAGs/ref_MAGs_SAGs/98LQset/bt2/all98LQset.fasta --interleaved $reads > $samplename-spRep.sam


# BAM, sort, index
samtools view -S -b  $samplename-spRep.sam >  $samplename-spRep.bam
samtools sort  $samplename-spRep.bam -o  $samplename-spRep.sorted.bam
samtools index $samplename-spRep.sorted.bam $samplename-spRep.sorted.bam.bai