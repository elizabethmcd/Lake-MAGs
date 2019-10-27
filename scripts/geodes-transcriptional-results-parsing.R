library(tximport)
library(tibble)
library(tidyverse)

# Merging kallisto count files per lake GEODES transcriptional time-series experiment
dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/transcriptomes"
# sample names
mendota <- read.table(file.path(dir, "mendota-samples.txt"), header=TRUE)
sparkling <- read.table(file.path(dir, "sparkling-samples.txt"), header=TRUE)
trout <- read.table(file.path(dir, "troutBog-samples.txt"), header=TRUE)
# filenames
mendota_files <- file.path(dir, mendota$sample, "abundance.tsv")
names(mendota_files) <- mendota$sample
sparkling_files <- file.path(dir, sparkling$sample, "abundance.tsv")
names(sparkling_files) <- sparkling$sample
troutBog_files <- file.path(dir, trout$sample, "abundance.tsv")
names(troutBog_files) <- trout$sample
# kallisto counts
mendota.txi.tsv <- tximport(mendota_files, type="kallisto", txOut=TRUE)
sparkling.txi.tsv <- tximport(sparkling_files, type="kallisto", txOut=TRUE)
troutBog.txi.tsv <- tximport(troutBog_files, type="kallisto", txOut=TRUE)
# save the normalized count (tpm) dfs
mendota.df <- as.data.frame(mendota.txi.tsv$abundance)
sparkling.df <- as.data.frame(sparkling.txi.tsv$abundance)
troutBog.df <- as.data.frame(troutBog.txi.tsv$abundance)
write.csv(mendota.df, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/mendota-normalized-tpm-counts.csv", quote=FALSE, row.names=FALSE)
write.csv(sparkling.df, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/sparkling-normalized-tpm-counts.csv", quote=FALSE, row.names=FALSE)
write.csv(troutBog.df, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/troutBog-normalized-tpm-counts.csv", quote=FALSE, row.names=FALSE)
# raw counts
mendota.raw.counts <- as.data.frame(mendota.txi.tsv$counts)
troutBog.raw.counts <- as.data.frame(troutBog.txi.tsv$counts)
sparkling.raw.counts <- as.data.frame(sparkling.txi.tsv$counts)

# with genome name as extra column for total aggregations
#raw
mendota.raw.names <- rownames_to_column(mendota.raw.counts, var="locus_tag")
sparkling.raw.names <- rownames_to_column(sparkling.raw.counts, var="locus_tag")
trout.raw.names <- rownames_to_column(troutBog.raw.counts, var="locus_tag")
#normalized
mendota.names <- rownames_to_column(mendota.df, var="locus_tag")
sparkling.names <- rownames_to_column(sparkling.df, var="locus_tag")
trout.names <- rownames_to_column(troutBog.df, var="locus_tag")
    # mendota raw
mendota.raw.split <- mendota.raw.names %>% separate(locus_tag, c("genome"), sep="_") %>% cbind(mendota.raw.names$locus_tag)
colnames(mendota.raw.split)[41] <- c("locus_tag")
mendota.raw.final <- mendota.raw.split[,c(1,41,2:40)]
    # mendota norm
mendota.split <- mendota.names %>% separate(locus_tag, c("genome"), sep='_') %>% cbind(mendota.names$locus_tag)
colnames(mendota.split)[41] <- c("locus_tag")
mendota.final <- mendota.split[,c(1,41,2:40)]
    # sparkling raw
sparkling.raw.split <- sparkling.raw.names %>% separate(locus_tag, c("genome"), sep="_") %>% cbind(sparkling.raw.names$locus_tag)
colnames(sparkling.raw.split)[41] < -c("locus_tag")
sparkling.raw.final <- sparkling.raw.split[,c(1,41,2:40)]
    # sparkling norm
sparkling.split <- sparkling.names %>% separate(locus_tag, c("genome"), sep='_') %>% cbind(sparkling.names$locus_tag)
colnames(sparkling.split)[41] <- c("locus_tag")
sparkling.final <- sparkling.split[,c(1,41,2:40)]
    # trout Bog raw
trout.raw.split <- trout.raw.names %>% separate(locus_tag, c("genome"), sep='_') %>% cbind(trout.raw.names$locus_tag)
colnames(trout.raw.split)[33] <- c("locus_tag")
trout.raw.final <- trout.raw.split[,c(1,33,2:32)]
    # trout Bog norm
troutBog.split <- trout.names %>% separate(locus_tag, c("genome"), sep='_') %>% cbind(trout.names$locus_tag)
colnames(troutBog.split)[33] <- c("locus_tag")
troutBog.final <- troutBog.split[,c(1,33,2:32)]

# aggregate raw counts for each dataset
  # mendota
mendota.raw.bins <- aggregate(mendota.raw.final[3:41], list(mendota.raw.final$genome), sum)
mendota.raw.bins$sum <- rowSums(mendota.raw.bins[2:40])
mendota.raw.bins$avg <- rowMeans(mendota.raw.bins[2:40])
  # sparkling
sparkling.raw.bins <- aggregate(sparkling.raw.final[3:41], list(sparkling.raw.final$genome), sum)
sparkling.raw.bins$sum <- rowSums(sparkling.raw.bins[2:40])
sparkling.raw.bins$avg <- rowMeans(sparkling.raw.bins[2:40])
  # tb and take out NaN sample columns
trout.raw.bins <- aggregate(trout.raw.final[3:33], list(trout.raw.final$genome), sum) %>% select(-c(GEODES068, GEODES072))
trout.raw.bins$sum <- rowSums(trout.raw.bins[2:30])
trout.raw.bins$avg <- rowMeans(trout.raw.bins[2:30])

# aggregate normalized counts for each dataset by bin (avg and sum of counts)
  # mendota
mendota.bins <- aggregate(mendota.final[3:41], list(mendota.final$genome), sum)
mendota.bins$sum <- rowSums(mendota.bins[2:40])
mendota.bins$avg <- rowMeans(mendota.bins[2:40])
  # sparkling
sparkling.bins <- aggregate(sparkling.final[3:41], list(sparkling.final$genome), sum)
sparkling.bins$sum <- rowSums(sparkling.bins[2:40])
sparkling.bins$avg <- rowMeans(sparkling.bins[2:40])
  # trout
      # take out sames with weird NaN mapping & double check
trout.bins <- aggregate(troutBog.final[3:33], list(troutBog.final$genome), sum) %>% select(-c(GEODES068, GEODES072))
trout.bins$sum <- rowSums(trout.bins[2:30])
trout.bins$avg <- rowMeans(trout.bins[2:30])
# save to select certain bins
write.csv(mendota.bins, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/mendota-bin-norm-counts.csv", quote=FALSE, row.names=FALSE)
write.csv(sparkling.bins, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/sparkling-bin-norm-counts.csv", quote=FALSE, row.names=FALSE)
write.csv(trout.bins, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/troutBog-bin-norm-counts.csv", quote=FALSE, row.names=FALSE)

# combine raw and normalized counts columns of totals and averages for comparison
mendota.combined <- mendota.bins %>% select(Group.1, sum, avg) %>% cbind(mendota.raw.bins$sum) %>% cbind(mendota.raw.bins$avg)
sparkling.combined <- sparkling.bins %>% select(Group.1, sum, avg) %>% cbind(sparkling.raw.bins$sum) %>% cbind(sparkling.raw.bins$avg)
trout.combined <- trout.bins %>% select(Group.1, sum, avg) %>% cbind(trout.raw.bins$sum) %>% cbind(trout.raw.bins$avg)
comb.names <- c("genome", "norm_sum", "norm_avg", "raw_sum", "raw_avg")
colnames(mendota.combined) <- comb.names
colnames(sparkling.combined) <- comb.names
colnames(trout.combined) <- comb.names

# combine with metadata and metaG mapping
metadata <- read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/metadata/all-geodes-stats-mapping.csv")
colnames(metadata)[1] <- c("genome")
all.bins.norm.counts <- rbind(mendota.combined, sparkling.combined, trout.combined)
merged.metadata.counts <- left_join(metadata, all.bins.norm.counts)
# save
write.csv(merged.metadata.counts, "/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/metadata/all-geodes-metaT-metadata.csv", row.names=FALSE, quote=FALSE)

# geodes 116 thermo bin
thermo <- troutBog.final %>% filter(genome == 'GEODES057-bin.68') %>% select(-c(GEODES068, GEODES072))
write.csv(thermo, "/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/thermoleo/metadata/annotations/geodes-thermo-all-counts.csv", quote=FALSE, row.names=FALSE)
thermo$sum <- rowSums(thermo[3:31])
thermo$avg <- rowMeans(thermo[3:31])
thermo.no.zeros <- thermo %>% filter(sum > 0)
write.csv(thermo.no.zeros, "/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/thermoleo/metadata/annotations/geodes-thermo-no-zeros.csv", quote=FALSE, row.names=FALSE)
colSums(thermo[3:33])
colMeans(thermo[3:33])

# plot of avg expression by bin for each lake
  # mendota
mendota.classf <- merged.metadata.counts %>% filter(grepl('GEODES11', genome)) %>% select(genome, classification, norm_avg)
mendota.cleaned <- separate(mendota.classf, classification, into=c("Domain", "Phyla", "Class", "Order", "Family", "Genus", "Species"), sep=';') %>% select(genome, Phyla, norm_avg)
men.plt <- ggplot(mendota.cleaned, aes(x=reorder(genome,-norm_avg), y=norm_avg, fill=Phyla)) + geom_col() + labs(x="Bin", y="Average Normalized Expressoin (TPM)") + theme_classic() + theme(axis.text.x=element_text(angle=85, vjust=0.5)) + scale_fill_brewer(palette="Paired")
men.plt
ggsave(filename="/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/mendota-avg-expression-bins.png", plot=men.plt, height=15, width=25, units=c("cm"))
  # trout
trout.classf <- merged.metadata.counts %>% filter(grepl('GEODES05', genome)) %>% select(genome, classification, norm_avg)
trout.cleaned <- separate(trout.classf, classification, into=c("Domain", "Phyla", "Class", "Order", "Family", "Genus", "Species"), sep=';') %>% select(genome, Phyla, norm_avg)
trout.plt <- ggplot(trout.cleaned, aes(x=reorder(genome,-norm_avg), y=norm_avg, fill=Phyla)) + geom_col() + labs(x="Bin", y="Average Normalized Expressoin (TPM)") + theme_classic() + theme(axis.text.x=element_text(angle=85, vjust=0.5)) + scale_fill_brewer(palette="Paired")
trout.plt
ggsave(filename="/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/trout-avg-expression-bins.png", plot=trout.plt, height=15, width=25, units=c("cm"))
  # sparkling
sparkling.classf <- merged.metadata.counts %>% filter(grepl('GEODES00', genome)) %>% select(genome, classification, norm_avg)
sparkling.cleaned <- separate(sparkling.classf, classification, into=c("Domain", "Phyla", "Class", "Order", "Family", "Genus", "Species"), sep=';') %>% select(genome, Phyla, norm_avg)
sparkling.plt <- ggplot(sparkling.cleaned, aes(x=reorder(genome,-norm_avg), y=norm_avg, fill=Phyla)) + geom_col() + labs(x="Bin", y="Average Normalized Expressoin (TPM)") + theme_classic() + theme(axis.text.x=element_text(angle=85, vjust=0.5)) + scale_fill_brewer(palette="Paired")
sparkling.plt
ggsave(filename="/Users/emcdaniel/Desktop/McMahon-Lab/Lake-MAGs/results/transcriptomes/sparkling-avg-expression-bins.png", plot=sparkling.plt, height=15, width=25, units=c("cm"))

# merge with Prokka annotations 
annots <- read_delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/thermoleo/metadata/annotations/GEODES057-bin.68-annots.txt", delim="\t", col_names=TRUE)
colnames(annots) <- c("locus_tag", "annotation")
thermo.counts.annots <- left_join(thermo.no.zeros, annots)
thermo.no.ribosomal <- thermo.counts.annots %>% filter(!grepl('ribosomal', annotation))
thermo.no.hypothetical <- thermo.no.ribosomal %>% filter(!grepl('hypothetical', annotation))
# automatically pulls out tRNAs
thermo.no.missing <- thermo.no.hypothetical %>% filter(!is.na(annotation))
