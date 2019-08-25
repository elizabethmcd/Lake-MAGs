# Preliminary quality checks of first round of GEODES and MENDH bins

library(tidyverse)

# stats
geodes_stats <- read.delim("metadata/bin_stats/geodes-preliminary-noS006-bins-qual-stats.txt", sep=" ", header=FALSE)
mendota_stats <- read.delim("metadata/bin_stats/mendota-bins-checkm-stats.txt", sep=" ", header=FALSE)

colnames(geodes_stats) <- c("bin", "classf", "size", "completeness", "redundancy")
colnames(mendota_stats) <- c("bin", "classf", "size", "completeness", "redundancy")

# only work with bins above 50% complete and less than 10% redundant
geodes_medium <- geodes_stats %>% filter(completeness > 50 & redundancy < 10)
mendota_medium <- mendota_stats %>% filter(completeness > 50 & redundancy < 10)

# for curiosity sake, lists of high qual = abofe 90% complete and less than 5% redundant to keep in mind for transcriptional mapping
geodes_highqual <- geodes_stats %>% filter(completeness > 90 & redundancy < 5)
# can't map RNA reads to these because don't have hypo samples from GEODES, but good to know regardless
mendota_highqual <- mendota_stats %>% filter(completeness > 90 & redundancy < 5)

# write out medium stats
write_delim(geodes_medium, 'metadata/bin_stats/geodes-medium-qual-bin-stats.txt', delim="\t")
write_delim(mendota_medium, 'metadata/bin_stats/mendota-medium-qual-bin-stats.txt', delim="\t")
