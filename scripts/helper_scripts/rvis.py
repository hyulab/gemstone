# Title: gene_expression.py
# Description: This script is to annotate gene RVIS (Residual Variation Intolerance Score) for genes.
# Author: Siwei Chen
# Date created: 14/03/2016
# Date last modified: 14/03/2016
# Python version: 2.7.5

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="input file")
args = parser.parse_args()

# Store RVIS scores for genes
gene_rvis = dict([line.strip().split("\t")[0],line.strip().split("\t")[1:]] for line in open("../../GeMSTONE-data/RVIS_score_full.txt").readlines())
# Read in variants for annotation
# Write out variants with annotated RVIS
fh2rvis = open(args.f).readlines()
fh_rvis = open(args.f[:-3] + 'rvis.txt','w')
fh_rvis.write(fh2rvis[0][:-1]+'\t'+'\t'.join(gene_rvis['GENE']) + '\n')
for line in fh2rvis[1:]:
    gene = line.strip().split("\t")[0].upper()
    rvis = gene_rvis[gene] if gene in gene_rvis else ['.']*len(gene_rvis['GENE'])
    fh_rvis.write(line[:-1]+'\t'+'\t'.join(rvis) + '\n')
fh_rvis.close()
        