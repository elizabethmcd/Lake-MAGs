#!/bin/bash

mkdir input
mv $1 input

# Running lineage specific marker set

checkm lineage_wf input output

tax=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /marker lineage/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
length=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Genome size/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
complete=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Completeness/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')
contamination=$(awk -F', ' '{for(i=1;i<=NF;i++){if ($i ~ /Contamination/){print $i}}}' output/storage/bin_stats_ext.tsv | awk -F': ' '{print $2}')

echo "${1%.*}" $tax $length $complete $contamination > "${sample%.*}".checkm.lineage.txt

rm -r output/