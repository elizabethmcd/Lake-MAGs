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