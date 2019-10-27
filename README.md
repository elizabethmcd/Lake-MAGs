# Assembly and Analysis of Lake MAGs

This repository contains scripts and workflows for reassembling and binning metagenomes from Trout Bog, Sparkling Lake, and the epilimnion/hypolimnion of Lake Mendota from the GEODES and MEHG metagenomics/metatranscriptomics projects. 

## Dependencies

- BBTools (bbuk, bbmap)
- MetaBAT
- CheckM
- dRep
- GTDB-tk
- Prokka
- Kallisto
- KofamKOALA

## Metagenomic Binning, Refinement, and Curation Workflow 

1. Map reads to assemblies with `map-metas-to-refs.sub` with a queue file of the names of the metagenomic assemblies
2. Bin with `geodes-binning.sub`
3. Run CheckM with Sarah Stevens' [checkm-chtc-pipeline](https://github.com/sstevens2/checkm-chtc-pipeline)
4. Provide the CheckM quality information to `drep` to dereplicate bins between samples of the same lake. 
- Install `drep` by creating a conda environment with `conda create -n drep drep`
- From the combined stats files of all bins, get the ones that are only medium quality and put into a CSV file for input into `drep` with `cat all-geodes-stats.txt | awk -F " " '{if (($4 > 50) && ($5 < 10)) print $1".fna,"$4","$5}' > med-geodes-input-drep.txt`. Add the header `genome,completeness,contamination`
- Run `drep` with `dRep dereplicate outputDir -g inputDir/*.fna --genomeInfo bin-stats.csv
5. Map metagenomic reads back to assembled bins, check uniform coverage with anvi'o, run final quality statistics
6. Classify automatically with GTDB-tk and manual ribosomal protein classification workflow
7. Run functional annotation and file formatting with Prokka, other databases of interest

## Metatranscriptomics Workflow 

Before using the kallisto metatranscriptomics pipeline, each RNAseq experiment must be quality filtered, rRNA removed, and files deinterleaved to work with kallisto. The latter can be done with `queue-deinterleave.sh` to queue the set of fastq files to be deinterleaved with the `deainterleave_fastq.sh` script. 

1. Install kallisto with `conda create -n kallisto kallisto`
2. Create reference mapping database with `kallisto index`
3. Map RNAseq reads to the reference mapping database with `kallisto quant`
4. Quantify and normalize reads (R package tximport)
5. Functional annotation (KofamKOALA), or BlastKOALA of the collection of genes

# IMPORTANT NOTES

I mixed up with metagenomic samples came from Mendota and Sparkling, therefore the queue metadata files and some scripts for calling the metagenomes don't match what the samples actually are. The naming scheme doesn't mess up the actual assemblies/bins because they were all kept clustered together by sample names (such as GEODES005 and GEODES006 would still be the same whether or not I assigned them to Mendota or Sparkling), but this will cause issues with mapping the metatranscriptomes, since those don't have the same 5/6, 57/58, 117/118 naming scheme. The metaG mapping is correct since those are queued by the filenames themselves (which refer to the coded sample names and not the names of the lake). However, equating lakes with different codes for mapping the metaTs to the bins will mess you up, so make sure you are mapping to the correct things (or you will get next to zero mapping except for really conserved genes if you mix up assembled bins vs metaT sample). For reference, the coded metaGs and bins correspond to the following lakes: 

- GEODES005/006 : Sparkling Lake
- GEODES057/058 : Trout Bog Lake
- GEODES117/118 : Lake Mendota