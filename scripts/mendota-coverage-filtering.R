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

write.csv(genome_sample_counts, "results/mendota_historical/inStrain/mendota_mag_coverage_samples.csv", quote = FALSE, row.names = FALSE)
