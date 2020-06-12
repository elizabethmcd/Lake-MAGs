# Organizing re-assembly of Lake MAGs

library(tidyverse)

# Trout Bog epi 

tb_epi_checkm <- read.delim("results/troutBog_epi_historical/troutBog-epi-checkm-stats.tsv", sep = "\t", header = FALSE)
tb_epi_gtdb <- read.delim("results/troutBog_epi_historical/troutBog-epi-gtdbtk.tsv", sep="\t")
tb_epi_bins <- read.table("results/troutBog_epi_historical/troutBog-epi-finalBins.txt")
colnames(tb_epi_checkm) <- c("user_genome", "lineage", "completeness", "contamination", "size_mbp", "contigs", "percent_gc", "un")
colnames(tb_epi_bins) <- c("user_genome")
tb_epi_merged <- left_join(tb_epi_bins, tb_epi_gtdb)
tb_epi_table <- left_join(tb_epi_merged, tb_epi_checkm) %>% select(-un, -lineage)

# Trout Bog hypo 

tb_hypo_checkm <- read.delim("results/troutBog_hypo_historical/troutBog-hypo-checkM.tsv", sep="\t", header = FALSE)
tb_hypo_gtdb <- read.delim("results/troutBog_hypo_historical/troutBog-hypo-gtdb.tsv", sep="\t")
tb_hypo_bins <- read.table("results/troutBog_hypo_historical/troutBog-hypo-finalBins.txt")
colnames(tb_hypo_checkm) <- c("user_genome", "lineage", "completeness", "contamination", "size_mbp", "contigs", "percent_gc", "un")
colnames(tb_hypo_bins) <- c("user_genome")
tb_hypo_merged <- left_join(tb_hypo_bins, tb_hypo_gtdb)
tb_hypo_table <- left_join(tb_hypo_merged, tb_hypo_checkm) %>% select(-un, -lineage)

# Crystal Bog

crystal_checkm <- read.delim("results/crystalBog_historical/new_spades_assemb/crystalBog-new-checkm-stats.tsv", sep="\t", header = FALSE)
crystal_gtdb <- read.delim("results/crystalBog_historical/new_spades_assemb/crystalBog-new-gtdb.tsv", sep="\t")
crystal_bins <- read.table("results/crystalBog_historical/new_spades_assemb/crystalBog-new-finalBins.txt")
colnames(crystal_checkm) <- c("user_genome", "lineage", "completeness", "contamination", "size_mbp", "contigs", "percent_gc", "un")
colnames(crystal_bins) <- c("user_genome")
crystal_merged <- left_join(crystal_bins, crystal_gtdb)
crystal_table <- left_join(crystal_merged, crystal_checkm) %>% select(-un, -lineage)

# Mary

mary_checkm <- read.delim("results/mary_historical/new_spades_assemb/mary-checkm-stats.tsv", sep="\t", header=FALSE)
mary_gtdb <- read.delim("results/mary_historical/new_spades_assemb/mary-gtdb-classifications.tsv", sep="\t")
mary_bins <- read.table("results/mary_historical/new_spades_assemb/mary-new-finalBins.txt")
colnames(mary_checkm) <- c("user_genome", "lineage", "completeness", "contamination", "size_mbp", "contigs", "percent_gc", "un")
colnames(mary_bins) <- c("user_genome")
mary_merged <- left_join(mary_bins, mary_gtdb)
mary_table <- left_join(mary_merged, mary_checkm) %>% select(-un, -lineage)

write.csv(tb_epi_table, "results/troutBog_epi_historical/troutBog_epi_historical_metadata_final.csv", quote=FALSE, row.names = FALSE)
write.csv(tb_hypo_table, "results/troutBog_hypo_historical/troutBog_hypo_historical_metadata_final.csv", quote=FALSE, row.names=FALSE)
write.csv(crystal_table, "results/crystalBog_historical/new_spades_assemb/crystalBog_historical_metadata_final.csv", quote=FALSE, row.names = FALSE)
write.csv(mary_table, "results/mary_historical/new_spades_assemb/mary_historical_metadata_final.csv", quote=FALSE, row.names = FALSE)
