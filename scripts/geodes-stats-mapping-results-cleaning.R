library(tidyverse)
library(reshape2)

# merge stats and gtdbtk metadata

stats = read.delim("results/GEODES/geodes_bin_stats/all-stats.txt", sep="\t")
gtdb = read.delim("results/GEODES/geodes_bin_stats/gtdbtk.bac120.summary.tsv", sep="\t")
colnames(gtdb)[1] = c("bin")
names = gtdb[,c(1,2)]
merged = left_join(names, stats)
write.csv(merged, 'results/GEODES/geodes_bin_stats/all-geodes-final-stats.csv', quote=FALSE, row.names = FALSE)

# sparkling stats
sparkling_stats = merged %>% filter(grepl('GEODES005|GEODES006', bin))
write.csv(sparkling_stats, 'results/GEODES/geodes_bin_stats/sparkling-bins-stats-final.csv', quote=FALSE, row.names = FALSE)
# tb stats
tb_stats = merged %>% filter(grepl('GEODES057|GEODES058', bin))
write.csv(tb_stats, 'results/GEODES/geodes_bin_stats/tb-bins-stats-final.csv', quote= FALSE, row.names = FALSE)
# mendota stats
mendota_stats = merged %>% filter(grepl('GEODES117|GEODES118', bin))
write.csv(mendota_stats, 'results/GEODES/geodes_bin_stats/mendota-bins-stats-final.csv', quote=FALSE, row.names = FALSE)

# mapping files
mendota = read.delim("results/GEODES/GEODES_mapping_results/mendota-mapping-results-v2.txt", sep="\t")
troutBog = read.delim("results/GEODES/GEODES_mapping_results/troutBog-mapping-results.txt", sep="\t")
sparkling = read.delim("results/GEODES/GEODES_mapping_results/sparkling-mapping-results-v2.txt", sep="\t")

# mendota
men = mendota[,c(2,3,8)]
men.s = spread(men, meta, AvgCov)
colnames(men.s)[1] = c("bin")
men_mapping_stats = left_join(men.s, merged)
write.csv(men_mapping_stats, "metadata/mapping/mendota-mapping-stats-metadata.txt", quote=FALSE, row.names=FALSE)

# tb
trout = troutBog[,c(2,3,8)]
trout.s = spread(trout, meta, AvgCov)
colnames(trout.s)[1] = c("bin")
trout_mapping_stats = left_join(trout.s, merged)
write.csv(trout_mapping_stats, "metadata/mapping/troutBog-mapping-stats-metadata.txt", quote=FALSE, row.names=FALSE)

# spark
spark = sparkling[,c(2,3,8)]
spark.s = spread(spark, meta, AvgCov)
colnames(spark.s)[1] = c("bin")
spark_mapping_stats = left_join(spark.s, merged)
write.csv(spark_mapping_stats, "metadata/mapping/sparkling-mapping-stats-metadata.txt", quote=FALSE, row.names=FALSE)

# combined stats and mapping file for all three lakes
men.c = men.s
colnames(men.c) = c("bin", "metaG1", "metaG2")
trout.c = trout.s
colnames(trout.c) = c("bin", "metaG1", "metaG2")
spark.c = spark.s
colnames(spark.c) = c("bin", "metaG1", "metaG2")

combined = rbind(men.c, trout.c, spark.c)
joined_stats_mapping = left_join(merged, combined)

write.csv(joined_stats_mapping, "metadata/all-geodes-stats-mapping.csv", quote=FALSE, row.names=FALSE)
