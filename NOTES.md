# IMPORTANT NOTES


## Transcriptomics notes for GEODES

I mixed up with metagenomic samples came from Mendota and Sparkling, therefore the queue metadata files and some scripts for calling the metagenomes don't match what the samples actually are. The naming scheme doesn't mess up the actual assemblies/bins because they were all kept clustered together by sample names (such as GEODES005 and GEODES006 would still be the same whether or not I assigned them to Mendota or Sparkling), but this will cause issues with mapping the metatranscriptomes, since those don't have the same 5/6, 57/58, 117/118 naming scheme. The metaG mapping is correct since those are queued by the filenames themselves (which refer to the coded sample names and not the names of the lake). However, equating lakes with different codes for mapping the metaTs to the bins will mess you up, so make sure you are mapping to the correct things (or you will get next to zero mapping except for really conserved genes if you mix up assembled bins vs metaT sample). For reference, the coded metaGs and bins correspond to the following lakes: 

- GEODES005/006 : Sparkling Lake
- GEODES057/058 : Trout Bog Lake
- GEODES117/118 : Lake Mendota