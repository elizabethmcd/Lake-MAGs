library(tidyverse)

# Historical Lake MAGs stats organization for Mendota, Trout (epi/hypo), Mary, and Crystal Bog

# Mendota
mend_checkm = read_delim('results/mendota_historical/mendota_checkm_stats.tsv', delim="\t", col_names=TRUE)
colnames(mend_checkm) = c('Bin', 'lineage', 'Completeness', 'Contamination', "size_bp", "n_contigs", 'percent_gc', "misc")
mend_classf = read_delim('results/mendota_historical/mendota_gtdb_classifications.txt', delim="\t")
colnames(mend_classf) = c('Bin', 'Classification')
mend_merged = left_join(mend_classf, mend_checkm) %>% select(Bin, Classification, Completeness, Contamination, size_bp, n_contigs, percent_gc)
write.csv(mend_classf, 'results/mendota_historical/mendota-dRep-bins-table.csv', quote = FALSE, row.names = FALSE)

# Trout epi
trout_epi_checkm = read_delim('results/troutBog_epi_historical/troutBog-epi-checkm-stats.tsv', delim="\t", col_names=FALSE)
colnames(trout_epi_checkm) = c('Bin', 'lineage', 'Completeness', 'Contamination', "size_bp", "n_contigs", 'percent_gc', "misc")
trout_epi_classf = read_delim('results/troutBog_epi_historical/troutBog-epi-gtdbtk.tsv', delim="\t")
colnames(trout_epi_classf) = c('Bin', 'Classification')
trout_epi_merged = left_join(trout_epi_classf, trout_epi_checkm) %>% select(Bin, Classification, Completeness, Contamination, size_bp, n_contigs, percent_gc)
write.csv(trout_epi_merged, 'results/troutBog_epi_historical/troutEpi-dRep-bins-table.csv', quote=FALSE, row.names = FALSE)

# Trout hypo
trout_hypo_checkm = read_delim('results/troutBog_hypo_historical/troutBog-hypo-checkM.tsv', delim="\t", col_names=FALSE)
colnames(trout_hypo_checkm) = c('Bin', 'lineage', 'Completeness', 'Contamination', "size_bp", "n_contigs", 'percent_gc', "misc")
trout_hypo_classf = read_delim('results/troutBog_hypo_historical/troutBog-hypo-gtdb.tsv', delim="\t")
colnames(trout_hypo_classf) = c('Bin', 'Classification')
trout_hypo_merged = left_join(trout_hypo_classf, trout_hypo_checkm) %>% select(Bin, Classification, Completeness, Contamination, size_bp, n_contigs, percent_gc)
write.csv(trout_hypo_merged, 'results/troutBog_hypo_historical/troutHypo-dRep-bins-table.csv', quote=FALSE, row.names = FALSE)

# Mary
mary_checkm = read_delim('results/mary_historical/mary-checkm-stats.tsv', delim="\t", col_names = FALSE)
colnames(mary_checkm) = c('Bin', 'lineage', 'Completeness', 'Contamination', "size_bp", "n_contigs", 'percent_gc', "misc")
mary_classf = read_delim('results/mary_historical/mary-gtdb-classfs.tsv', delim="\t")
colnames(mary_classf) = c('Bin', 'Classification')
mary_merged = left_join(mary_classf, mary_checkm) %>% select(Bin, Classification, Completeness, Contamination, size_bp, n_contigs, percent_gc)
write.csv(mary_merged, 'results/mary_historical/mary-dRep-bins-table.csv', quote=FALSE, row.names = FALSE)

# Crystal Bog
crystal_checkm = read_delim('results/crystalBog_historical/crystalBog_checkm_stats.tsv', delim="\t")
colnames(crystal_checkm) = c('Bin', 'lineage', 'Completeness', 'Contamination', "size_bp", "n_contigs", 'percent_gc', "misc")
crystal_classf = read_delim('results/crystalBog_historical/crystalBog_gtdb_classf.tsv', delim="\t")
colnames(crystal_classf) = c('Bin', 'Classification')
crystal_merged = left_join(crystal_classf, crystal_checkm) %>% select(Bin, Classification, Completeness, Contamination, size_bp, n_contigs, percent_gc)
write.csv(crystal_merged, 'results/crystalBog_historical/crystalBog-dRep-bins-table.csv', quote=FALSE, row.names = FALSE)

# Combining lakes together into one dataset for redone historical metagenomes binning
mendota = data.frame(lake = 'Mendota', mend_merged)
troutBog_epi = data.frame(lake = "TroutBog_epi", trout_epi_merged)
troutBog_hypo = data.frame(lake = "TroutBog_hypo", trout_hypo_merged)
mary = data.frame(lake = "Mary", mary_merged)
crystal = data.frame(lake = "CrystalBog", crystal_merged)
allLakes = rbind(mendota, troutBog_epi, troutBog_hypo, mary, crystal)
write.csv(allLakes, 'results/allLakes-historical-rebinned-dRep-table.csv', row.names = FALSE, quote=FALSE)
