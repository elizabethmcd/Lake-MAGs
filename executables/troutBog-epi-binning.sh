#! /bin/bash

######################
# MetaBAT Binning on GLBRC
# Elizabeth McDaniel
######################

# Queue assembly files
ASSEMB=$1
NAME=$(basename $ASSEMB .a.fna)
MAPDIR=/home/GLBRCORG/emcdaniel/Lake-MAGs/troutbog/epi/mappingResults
METABATDIR=/opt/bifxapps/metabat-2.12.1/
OUTBIN=/home/GLBRCORG/emcdaniel/Lake-MAGs/troutbog/epi/binningResults

# Get depth 
$METABATDIR/jgi_summarize_bam_contig_depths --outputDepth $OUTBIN/$NAME-depth.txt $MAPDIR/$NAME*.sorted.bam

# Run metabat
$METABATDIR/metabat2 -i $ASSEMB -a $OUTBIN/$NAME-depth.txt -o $OUTBIN/$NAME-bins/$NAME-bin
