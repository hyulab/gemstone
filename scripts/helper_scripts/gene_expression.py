# Title: gene_expression.py
# Description: This script is to annotate gene expression level (from GTEx and HPA) in different tissues for genes.
# Author: Siwei Chen
# Date created: 13/03/2016
# Date last modified: 13/03/2016
# Python version: 2.7.5

import imp
import argparse
import os
# import python scripts with pre-populated dictionaries with gene expression values 
from gtex_gene_expr import gtex_tissue_genes
from gtex_gene_expr import tissue_subregions
from gtex_gene_expr import gene_pref_tissue
from gtex_gene_expr import gene_sex_tissue
from gtex_gene_expr import gene_eth_tissue
from gtex_gene_expr import gene_age
from hpa_gene_expr import hpa_gene_tissue
from hpa_gene_expr import gene_subcellular
from hpa_gene_expr import gene_enri_tissue
from hpa_gene_expr import gene_enri_group
parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
args = parser.parse_args()

c = imp.load_source('cache', args.c).cache

# Annotate all file and then subset for recurrent file
fh=open("../query/output/{0}/Cross-sample_all_genes.txt".format(c['qid'])).readlines() 
of=open("../query/output/{0}/Cross-sample_all_expr.txt".format(c['qid']),'w')
hdr=fh[0][:-1]
if c['expression']['hpa'] != []:
    hdr += '\t'+'\t'.join(["TISSUE_EXPR_HPA","TISSUE_SUBCELLULAR_HPA","TISSUE_ENRICHMENT_HPA","TISSUE_GROUP_ENRICHMENT_HPA"])
if c['expression']['gtex'] != []:
    hdr += '\t'+'\t'.join(["TISSUE_EXPR_GTEx(sub-regions)",
                           "TISSUE_PREFERENCE_GTEx",
                           "TISSUE_SEX_DIFF_GTEx(log2FoldChange(females/males)_FDR)",
                           "TISSUE_ETHNICITY_DIFF_GTEx(log2FoldChange(AA/EA)_FDR)",
                           "GENE_EXPR_AGE_DIFF_GTEx(coefficient_FDR)"])
of.write(hdr+'\n')
for line in fh[1:]:
    content= line[:-1]
    gene= line.split("\t")[2]
    # Annotate gene expression(level) in selected tissues
    if c['expression']['hpa'] != []:
        expr=[tissue + '('+hpa_gene_tissue[gene][tissue]+')' for tissue in c['expression']['hpa'] if gene in hpa_gene_tissue and tissue in hpa_gene_tissue[gene]]
        content += '\t' + (', '.join(expr) if expr != [] else '.')
        # Annotate sub-cellular locations of the gene if expressed in selected tissues
        subcellular= gene_subcellular[gene] if gene in gene_subcellular and expr != [] else '.'
        content += '\t' + subcellular
        # Annotate tissue enrichment and tissue group enrichment
        enri_tissue= gene_enri_tissue[gene] if gene in gene_enri_tissue else '.'
        enri_group= gene_enri_group[gene] if gene in gene_enri_group else '.'
        content += '\t' + '\t'.join([enri_tissue,enri_group])
    if c['expression']['gtex'] != []:
        # Annotate gene expression(subregions) in selected tissues
        expr= []
        for tissue in c['expression']['gtex']:
            if tissue_subregions[tissue] != []:
                subregions=[sub for sub in tissue_subregions[tissue] if gene in gtex_tissue_genes[sub]]
                if subregions !=[]:
                    expr.append(tissue + '(' + ';'.join(subregions) + ')') 
            else:
                if gene in gtex_tissue_genes[tissue]:
                    expr.append(tissue)
        content += '\t' + (", ".join(expr) if expr != [] else '.')
        # Annotate tissue preferentially expressed genes, sex, enth, age diff 
        pref= ", ".join(gene_pref_tissue[gene]) if gene in gene_pref_tissue else '.'
        diffsex= ", ".join(gene_sex_tissue[gene]) if gene in gene_sex_tissue else '.'
        diffeth= ", ".join(gene_eth_tissue[gene]) if gene in gene_eth_tissue else '.'
        diffage= gene_age[gene] if gene in gene_age else '.'
        content += '\t' + '\t'.join([pref,diffsex,diffeth,diffage])
       
    of.write(content+'\n')
of.close()
os.system("""mv ../query/output/{0}/Cross-sample_all_expr.txt ../query/output/{0}/Cross-sample_all_genes.txt""".format(c['qid']))













