library(tidyverse)

# Mendota clusters check and quality tables
setwd('/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs')

# Quality stats
stats = read_delim('results/mendota_historical/mendota_checkm_stats.tsv', delim='\t')
colnames(stats) = c('bin', 'lineage', 'completeness', 'contamination', 'genome_size', 'contigs', 'percent_gc', 'unknown')
mendota_stats = stats %>% select(-unknown)

# GTDB classifications
gtdb = read_delim('results/mendota_historical/mendota_gtdb_classifications.txt', delim="\t")
colnames(gtdb) = c('bin','classification')

# ANI comparisons
ani = read_delim('results/mendota_historical/mendota.all.ani.cleaned', delim='\t')
colnames(ani) = c('bin1','bin2','ANI1','ANI2','AF1','AF2')
species_clusters = ani %>% filter(ANI1 > 95)
species_clusters$bin1 <- gsub(".fna","", species_clusters$bin1)
species_clusters$bin2 <- gsub(".fna","", species_clusters$bin2)

# merge preliminary stats and GTDB classfs
info = left_join(gtdb, mendota_stats)
mendota_info = info %>% select(-lineage)
# write out preliminary mendota information
write.csv(mendota_info, 'results/mendota_historical/mendota-mags-preliminary-dereplicated-metadata.csv', quote=FALSE, row.names = FALSE)

# merge cluster information with stats
cluster_stats = mendota_stats %>% select(-lineage, -percent_gc)
ani_stats1 = left_join(species_clusters, cluster_stats, by=c("bin1" = "bin"))
colnames(ani_stats1) <- c("bin1", "bin2", "ANI1", "ANI2", "AF1", "AF2", "comp1", "contam1", "size1", "contigs1")
ani_stats2 <- left_join(ani_stats1, cluster_stats, by=c("bin2" = "bin"))
colnames(ani_stats2) <- c("bin1", "bin2", "ANI1", "ANI2", "AF1", "AF2", "comp1", "contam1", "size1", "contigs1", "comp2", "contam2", "size2", "contigs2")
ani_clusters = left_join(ani_stats2, gtdb, by=c('bin1'='bin'))
cluster_table <- ani_clusters %>% select(bin1, bin2, classification, ANI1, ANI2, AF1, AF2, comp1, comp2, contam1, contam2, size1, size2, contigs1, contigs2)
# write out cluster table
write.csv(cluster_table, 'results/mendota_historical/mendota-ani-cluster-table.csv', quote = FALSE, row.names = FALSE)

# cyanobacteria information for metox/patricia group
cyanobacteria = mendota_info %>% filter(grepl('Cyanobacteria', classification))
write.csv(cyanobacteria, 'results/mendota_historical/mendota-historical-cyanobacteria-stats.csv', quote=FALSE, row.names = FALSE)
