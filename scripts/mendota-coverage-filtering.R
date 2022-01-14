library(tidyverse)

# Mendota coverage checks for genomes 
# read in all genome info files to one dataframe

genome_info_path <- "results/mendota_historical/inStrain/quick_profiles/"
files <- dir(genome_info_path, pattern="*.csv")
mendota_info <- data_frame(filename = files) %>%
  mutate(file_contents = map(filename, ~ read.csv(file.path(genome_info_path, .)))
  ) %>%
  unnest()

mendota_table <- mendota_info %>% 
  mutate(sample = gsub("-quick-profile.csv", "", filename)) %>% 
  select(genome, sample, coverage, breadth) %>% 
  filter(coverage > 10 & breadth > 0.9)

# combine with metadata 
mendota_metadata <- read.csv("results/mendota_historical/mendota_historical_finalBins_metadata.csv")
colnames(mendota_metadata)[1] <- c("genome")

mendota_table_info <- left_join(mendota_table, mendota_metadata) %>% 
  filter(completeness > 80 & contamination < 10) %>% 
  select(genome, gtdb_classification, completeness, contamination, sample, coverage, breadth)

mendota_table_info %>% 
  count(sample) %>% 
  arrange(desc(n))

genome_sample_counts <- mendota_table_info %>% 
  count(genome, gtdb_classification) %>% 
  arrange(desc(n))

mendota_table_info %>% 
  count(sample) %>% 
  arrange(desc(n))

write.csv(genome_sample_counts, "results/mendota_historical/inStrain/mendota_mag_coverage_samples.csv", quote = FALSE, row.names = FALSE)

mendota_queues <- mendota_table_info %>% 
  select(genome, sample) %>% 
  mutate(genome = paste(genome, "fa", sep=".")) %>% 
  mutate(sample = paste(sample, "mum-spRep.sorted.bam", sep="."))

write_tsv(mendota_queues, "metadata/queues/mendota-inStrain-queues.txt", col_names = FALSE)

# Mendota metagenomes metadata (library codes to Julian date)

mendota_metagenome_info <- read.csv("metadata/ref_mags_sags_data/mendota-tb-historical-metaGs-supp-updated-EAM-2020-05-13.csv") %>% 
  filter(lake == "Mendota") %>% 
  select(library_code, date) %>% 
  mutate(sample = library_code) %>% 
  select(sample, date)

mendota_coverage_info <- left_join(mendota_table_info, mendota_metagenome_info)

mendota_coverage_info %>% ggplot(aes(x=date, y=coverage)) + geom_point() + facet_wrap(~ genome)

top_lineages <- genome_sample_counts %>% 
  top_n(5, n) %>% 
  pull(genome)

top_lineages_coverage_dynamics_plot <- mendota_coverage_info %>% 
  filter(genome %in% top_lineages) %>% 
  drop_na() %>% 
  ggplot(aes(x=date, y=coverage)) + geom_point() + facet_wrap(~ gtdb_classification, ncol=1) + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

mendota_coverage_info$gtdb_classification <- as.character(mendota_coverage_info$gtdb_classification)

coverage_dynamics_plot <- mendota_coverage_info %>%
  drop_na() %>% 
  separate(gtdb_classification, into=c("taxonomy", "phylum"), sep="p__") %>% 
  mutate(phylum = gsub(";.*", "", phylum)) %>%
  ggplot(aes(x=date, y=coverage)) + geom_point(aes(color=phylum)) + theme_bw() + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

mendota_coverage_info %>% 
  drop_na() %>% 
  count(date) %>% 
  arrange(desc(n))

ggsave("figs/coverage-dynamics.png", coverage_dynamics_plot, width=15, height=5, units=c("in"))

ggsave("figs/top5-lineages-coverage-dynamics.png", top_lineages_coverage_dynamics_plot, width=15, height=7, units=c("in")) 
