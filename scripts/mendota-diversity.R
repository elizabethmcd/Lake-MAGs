library(tidyverse)
library(patchwork)
library(cowplot)

# Mendota diversity results

diversity_path <- "results/mendota_historical/inStrain/genome_info/"
diversity_files <- dir(diversity_path, pattern=".tsv")
mendota_div <- data_frame(filename = diversity_files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(diversity_path, .)))
  ) %>%
  unnest()

mendota_div_table <- mendota_div %>% 
  filter(coverage > 10 & breadth > 0.9) %>% 
  separate(filename, into=c("org", "sample"), sep="-vs-") %>% 
  mutate(sample = gsub(".IS_genome_info.tsv", "", sample)) %>% 
  select(genome, sample, coverage, breadth, nucl_diversity, r2_mean, d_prime_mean)

mendota_div_info <- left_join(mendota_div_table, mendota_metadata) %>% 
  select(genome, gtdb_classification, sample, coverage, breadth, nucl_diversity, r2_mean, d_prime_mean)

mendota_div_df <- left_join(mendota_div_info, mendota_metagenome_info) %>% 
  drop_na()

top_lineages_diversity_dynamics <- mendota_div_df %>% 
  filter(genome %in% top_lineages) %>% 
  ggplot(aes(x=date, y=nucl_diversity)) + geom_point() + facet_wrap(~ gtdb_classification, ncol=1) + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

top_lineages_coverage_dynamics <- mendota_div_df %>% 
  filter(genome %in% top_lineages) %>% 
  ggplot(aes(x=date, y=coverage)) + geom_point() + facet_wrap(~ gtdb_classification, ncol=1) + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

top_lineages_recombination_dynamics <- mendota_div_df %>% 
  filter(genome %in% top_lineages) %>% 
  ggplot(aes(x=date, y=r2_mean)) + geom_point() + facet_wrap(~ gtdb_classification, ncol=1) + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

top_lineages_recombination_dynamics

mendota_top_lineages_grid <- plot_grid(top_lineages_coverage_dynamics, top_lineages_diversity_dynamics, ncol=2)

ggsave("figs/top_lineages_diversity_dynamics.png", top_lineages_diversity_dynamics, width=15, height=5, units=c("in"))

ggsave("figs/top5-lineages-coverage-dynamics.png", top_lineages_coverage_dynamics, width=15, height=5, units=c("in"))

ggsave("figs/top_lineages_grid.png", mendota_top_lineages_grid, width=20, height=6, units=c("in"))


mendota_div_df_clean <- mendota_div_df %>%
  select(genome, gtdb_classification, sample, date, coverage, breadth, nucl_diversity, r2_mean, d_prime_mean)

write.csv(mendota_div_df_clean, "results/mendota_historical/inStrain/mendota-diversity-table.csv", quote=FALSE, row.names = FALSE)

# LQ98set diversity results 
lq98set_path <- "results/mendota_historical/LQ98set_inStrain/"
lq98set_diversity_files <- dir(lq98set_path, pattern=".tsv")
lq98set_div <- data_frame(filename = lq98set_diversity_files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(lq98set_path, .)))
  ) %>%
  unnest()

lq98set_div_table <- lq98set_div %>% 
  filter(coverage > 10 & breadth > 0.9) %>% 
  separate(filename, into=c("org", "sample"), sep="-vs-") %>% 
  mutate(sample = gsub(".IS_genome_info.tsv", "", sample)) %>% 
  select(genome, sample, coverage, breadth, nucl_diversity)
  # only 22 genome : sample pairs met the coverage/breadth threshold when using only paired-end reads mapping in exact matches 

# LQ98set greedy using all mapped reads
lq98set_greedy_path <- "results/mendota_historical/LQ98set_inStrain_greedy/"
lq98set_greedy_diversity_files <- dir(lq98set_greedy_path, pattern=".tsv")
lq98set_greedy_div <- data_frame(filename = lq98set_greedy_diversity_files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(lq98set_greedy_path, .)))
  ) %>%
  unnest()

lq98set_div_table_greedy <- lq98set_greedy_div %>% 
  filter(coverage > 10 & breadth > 0.9) %>% 
  separate(filename, into=c("org", "sample"), sep="-vs-") %>% 
  mutate(sample = gsub(".IS_genome_info.tsv", "", sample)) %>% 
  select(genome, sample, coverage, breadth, nucl_diversity)
  # greedy mapping has 519 genome : sample pairs

lq98set_div_info <- left_join(lq98set_div_table_greedy, mendota_metagenome_info) %>% drop_na()

lq98set_div_table_greedy %>% 
  count(genome) %>% 
  arrange(desc(n))

top10_lq98set <- lq98set_div_info %>% 
  count(genome) %>% 
  top_n(10, n) %>% 
  pull(genome)

lq98set_covg_plot <- lq98set_div_info %>% 
  filter(genome %in% top10_lq98set) %>% 
  ggplot(aes(x=date, y=coverage)) + geom_point() + facet_wrap(~ genome, ncol=1) + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

lq98set_nucldiv_plot <- lq98set_div_info %>%  
  filter(genome %in% top10_lq98set) %>% 
  ggplot(aes(x=date, y=nucl_diversity)) + geom_point() + facet_wrap(~ genome, ncol=1) + theme(axis.text.x=element_text(angle=85, vjust=0.5, hjust=0.5))

lq98set_grid <- plot_grid(lq98set_covg_plot, lq98set_nucldiv_plot, ncol=2)

ggsave("figs/mendota-lq98set-covg-nucldiv-plot.png", width=20, height=10, units=c("in"))

write.csv(lq98set_div_info, "results/mendota_historical/lq98set_mendota_greedy_diversity_table.csv", quote=FALSE, row.names = FALSE)
