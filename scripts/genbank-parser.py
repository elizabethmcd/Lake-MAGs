#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

# Arguments

parser = argparse.ArgumentParser(description = "Parse Prokka annotations from a genbank file to add external gene calls and functions to Anvi'o")
parser.add_argument('gbk_file', metavar='GBK', help='Annotation file from Prokka in Genbank format')
parser.add_argument('--gene-calls', default='gene_calls.txt', help='Output: External gene calls (Default: gene_calls.txt)')
parser.add_argument('--annotation', default='gene_annot.txt', help="Output: Functional annotation for external gene calls (Default: gene_annot.txt)")

args = parser.parse_args()

# Input and output files
GBK = args.gbk_file
OUT_CDS = open(args.gene_calls, "w")
OUT_ANNO = open(args.annotation, "w")

# Headers for Anvi'o output
OUT_CDS.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n")
OUT_ANNO.write("locus_tag\taccession\tfunction\n")

# Gene ID and e-value
gene_id = 1
e_value = "0"

# Parse the GBK file
for record in SeqIO.parse(GBK, "genbank"):
    contig = record.name
    for f in record.features:
        if f.type == 'CDS':
            span, inference, function, locus_tag = (f.location, f.qualifiers["inference"][0], f.qualifiers["product"][0], f.qualifiers["protein_id"][0])
            protein_id = locus_tag.replace("X:", "")
            source = inference.split(":")[1]
            version = inference.split(":")[2]
            beg = span.start # biopython already reads in the genbank file to start counting from 0, so don't have to subtract 1
            end = span.end
            strand = str(span.strand)
            r = strand.replace("-1", "r")
            direction = r.replace("1", "f")
            if (float(beg - end)/float(3)).is_integer() == True:
                partial = str(0)
            else:
                partial = str(1)
            try:
                gene_acc = f.qualifiers["gene"][0]
            except KeyError:
                gene_acc = ""
            # Write out to gene calls and annotation files
            OUT_CDS.write('%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n' %(gene_id, contig, beg, end, direction, partial, source, version))
            OUT_ANNO.write('%s\t%s\t%s\n' % (protein_id, gene_acc, function))
            gene_id = gene_id + 1