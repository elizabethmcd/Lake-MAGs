#!/bin/bash

mkdir input
archive=$1
archbase=$(basename $1)
tar -xvf $archbase
sampleBins=$(basename $archbase .tar.gz)
sample=$(basename $archbase -bins.tar.gz)
mkdir $sample_checkm
for filename in $refbase/*.fa; do cp $filename ${filename%.fa}.fna; done
for filename in $refbase/*.fna; do cp $filename input/; done

# Running lineage specific marker set

checkm lineage_wf input output

tax=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /marker lineage/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
length=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Genome size/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
complete=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Completeness/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
contamination=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Contamination/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')

echo "${1%.*}" $tax $length $complete $contamination > "${sample%.*}".checkm.lineage.txt

mv ${sample%.*}".checkm.lineage.txt $sample_checkm
mv output/ $sample_checkm
tar -czf $sample_checkm.tgz $sample_checkm/