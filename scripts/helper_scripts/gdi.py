# Title: gdi.py
# Description: This script is to annotate GDI (Gene Damage Index) for genes.
# Author: Siwei Chen
# Date created: 14/03/2016
# Date last modified: 14/03/2016
# Python version: 2.7.5

import argparse
import os
import imp
import MySQLdb
mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
parser.add_argument("-f", help="input file")
args = parser.parse_args()
c = imp.load_source('cache', args.c).cache

# Format output 
query_mapping = {'All diseases': 'All_diseases', 'Mendelian (general model)':'Mendelian', 'Mendelian (autosomal dominant)': 'Mendelian_AD', 'Mendelian (autosomal recessive)':'Mendelian_AR', 'Primary immunodeficiency (general model)':'PID', 'Primary immunodeficiency (autosomal dominant)':'PID_AD', 'Primary immunodeficiency (autosomal recessive)':'PID_AR', 'Cancer (general model)':'Cancer','Cancer (autosomal recessive)':'Cancer_AR','Cancer (autosomal dominant)':'Cancer_AD'}
diseases = ["All_diseases","Mendelian","Mendelian_AD","Mendelian_AR","PID","PID_AD","PID_AR","Cancer","Cancer_AR","Cancer_AD"]
query = ['GDI','GDI_Phred'] + [i for i in diseases if i in [query_mapping[j] for j in c['gdi']]]

# Read in variants for annotation
# Write out variants with annotated GDI scores
fh2gdi = open(args.f).readlines()
fh_gdi = open(args.f[:-3]+'gdi.txt','w')
fh_gdi.write(fh2gdi[0][:-1] + '\t' + '\t'.join(query) + '\n')
# Get GDI for genes from MySQL table 
for line in fh2gdi[1:]:
    gene = line.strip().split("\t")[0]
    cur.execute('select {0} from Gene_Damage_Index where Gene = %s;'.format(','.join(query)),gene)
    rows= cur.fetchall()
    if rows == ():
        gdi = ["."]*len(query)
    else:
        gdi = list(rows[0])
    fh_gdi.write(line[:-1] + '\t' + '\t'.join(gdi) + '\n')
fh_gdi.close()