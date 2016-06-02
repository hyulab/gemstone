# Title: corss_family.py
# Description: This script is to filter on recurrences of variants across families/sporadic samples.
# Author: Siwei Chen, Juan Felipe Beltran
# Date created: 30/08/2015
# Date last modified: 11/03/2016
# Python version: 2.7.5

import argparse
import os
import imp
from glob import glob
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
args = parser.parse_args()
c = imp.load_source('cache', args.c).cache

# Collect all variant files of individual family/sporadic sample 
files = glob("../query/output/{0}/*_variants.txt".format(c['qid']))
header = open(files[0], 'r').readline()

# Sort out occurences of each variant cross samples
fh_tmpSort_all = open('../query/output/{0}/tmpSort_all.txt'
                      .format(c['qid']), 'w')
if c['inheritance'] != 'Rec_compound':
    variants_to_individuals = defaultdict(list)
    [[variants_to_individuals[tuple(row[:6] + row[7:])].append(row[6]) for row in [line.split('\t') for line in open(f).read().strip().split('\n')[1:]]] for f in files]
    for variants in variants_to_individuals:
        identifiers = ','.join(['[' + individual + ']' for individual in variants_to_individuals[variants]])
        variants = list(variants)
        fh_tmpSort_all.write('\t'.join(variants[:6] + [identifiers] + variants[6:]) + '\n')
else:
    variants_to_gt_individuals = defaultdict(list)
    [[variants_to_gt_individuals[tuple(row[:5] + row[7:])].append(row[5:7]) for row in [line.split('\t') for line in open(f).read().strip().split('\n')[1:]]] for f in files]

    for variants in variants_to_gt_individuals:
        #n_individuals = len(variants_to_individuals[variants])
        genotypes= ','.join(['[' + identifiers[0] + ']' for identifiers in variants_to_gt_individuals[variants]])
        samples = ','.join(['[' + identifiers[1] + ']' for identifiers in variants_to_gt_individuals[variants]])
        variants = list(variants)
        fh_tmpSort_all.write('\t'.join(variants[:5] + [genotypes,samples] + variants[7:]) + '\n')     
fh_tmpSort_all.close()

# Generate an "all_variant" result file containing all variants occurred at least once cross all samples
fh_all = open('../query/output/{0}/Cross-sample_all_variants.txt'
              .format(c['qid']), 'w')
fh_all.write(header)
fh_all.close()
# Use sort to format result file
os.system(	"sort -V -k1,1 -k2,2 ../query/output/{0}/tmpSort_all.txt".format(c['qid']) +
           ">> ../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid']))
os.system(	"rm -f ../query/output/{0}/tmpSort_all.txt".format(c['qid']))

# Filter on recurrences of each variant from "all_variant" file to generate "recurrent_variant" result file
fh_2ocr= open("../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid'])).readlines()
fh_ocr= open("../query/output/{0}/Cross-sample_recurrent_variants.txt".format(c['qid']),'w')
fh_ocr.write(fh_2ocr[0])
for line in fh_2ocr[1:]:
    samples= line[:-1].split("\t")[6][1:-1].split("],[")
    fm= 0
    ss= 0
    for i in samples:
        if len(i.split(",")) > 1:
            fm += 1
        else:
            ss += 1
    if fm+ss >= c['rc_F_lower'] and fm+ss <= c['rc_F_upper'] and ss >= c['rc_I_lower'] and ss <= c['rc_I_upper'] and fm >= c['family']:
        fh_ocr.write(line)
            
fh_ocr.close()




