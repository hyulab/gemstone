# Title: gene_expression.py
# Description: This script is to calculate gene set enrichment in pathways (GO, KEGG, Reactome, BioCarta).
# Author: Siwei Chen
# Date created: 14/09/2015
# Date last modified: 15/03/2016
# Python version: 2.7.5


import operator
import MySQLdb
import imp
import argparse
from fisher import pvalue
import subprocess as sub

parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
parser.add_argument("-out", help="out file")
parser.add_argument("-target", help="target gene list")
parser.add_argument("-background", help="background gene list")
args = parser.parse_args()
c = imp.load_source('cache', args.c).cache
mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

fh_pea = open(args.out,'w')

# Get target gene list from filtered gene file,
# backgroud gene list option selected by user
target = list(set([g.strip() for g in open(args.target).readlines()]))
background = list(set([g.strip() for g in open(args.background).readlines()]))

# Format output
fh_pea.write('\t'.join(['DATABASE', 'PATHWAY', 'GENE_ID','P-VALUE', 'Q-VALUE']) + '\n')

# Define a function to calculate gene set enrichment (p-value and q-value)
# given pathway database, target gene list and backgroud gene list
def pathwayEnrichment(database,target,background):
    db_enrichment = {}
    cur.execute('select * from {0};'.format(database))
    rows = cur.fetchall()
    for row in rows:
        pathway = row[0]
        pathway_genes = row[2].split(",")
        target_hits = len([i for i in target if i in pathway_genes])
        target_hits_genes = [i for i in target if i in pathway_genes]
        background_hits = len([i for i in background if i in pathway_genes])
        pval = pvalue(target_hits, len(target)-target_hits, background_hits, len(background)-background_hits).right_tail
        db_enrichment[tuple([database,pathway])] = [",".join(target_hits_genes),pval,pval*len(rows)]
    return db_enrichment

# Map quries to MySQL tables
cgi_mysql = {'GO_bp':'GO_biological_process','GO_cc':'GO_cellular_component','GO_mf':'GO_molecular_function','KEGG':'KEGG_EntrezID','BIOCARTA':'BIOCARTA_EntrezID','REACTOME':'REACTOME_EntrezID'}

# Calculate gene set enrichment in pathways and write to result file
for db in c['pea']:
    results = []
    for result in pathwayEnrichment(cgi_mysql[db],target,background).items():
        results.append(list(result[0])+result[1])
    ranked_qval = sorted(results, key=operator.itemgetter(1))
    for result in ranked_qval:
        if result[-1] <= float(c['qval']):
            fh_pea.write("\t".join(result[:3]+["%.2E" % result[3],"%.2E" % result[4]]) + '\n')
            
fh_pea.close()




