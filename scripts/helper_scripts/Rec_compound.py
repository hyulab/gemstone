# Title: Rec_compound.py
# Description: This script is to run co-segregation analysis with recessive compound heterozygous inheritance model.
# Author: Siwei Chen
# Date created: 16/01/2016
# Date last modified: 26/04/2016
# Python version: 2.7.5

import os
import subprocess as sub
import imp
from itertools import groupby
import argparse

import MySQLdb
mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

vcftools = 'vcftools'

parser = argparse.ArgumentParser()
parser.add_argument("-para", help="parameter file")
args = parser.parse_args()


def bash_pipe(command):
    print
    print "=" * 80
    print command
    print "-" * 80
    os.system(command)
    print "=" * 80

filename = args.para
c = imp.load_source('cache', filename).cache

ped_contents = [line.split('\t') for line in open(
    '../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n')]
individual_info = dict([(vector[1], vector[2:])
                        for vector in ped_contents])
family_members = dict([(key, list([g[1] for g in group])) for key, group in groupby(
    sorted([(p[0], p[1]) for p in ped_contents]), lambda x: x[0])])

samples = sub.Popen(["""bcftools query -l ../query/output/{0}/{0}.alts.filter.eff.vcf""".format(c['qid'])], shell=True, stdout=sub.PIPE, stderr=sub.PIPE).communicate()[0].strip().split("\n")
ped_samples = [line.strip() for line in open("../query/output/{0}/{0}.indv.txt".format(c['qid'])).readlines()]

joint_control = 'Y' if [i for i in samples if i not in ped_samples] else ''
control_vcf = 'Y' if c['control_filename'][-7:] == '.vcf.gz' or c['control_filename'][-4:] == '.vcf' else ''
control_txt = 'Y' if c['control_filename'][-4:] == '.txt' else ''

if c['control_filename'] == '':
    os.system('cp ../query/output/{0}/{0}.alts.filter.eff.vcf ../query/output/{0}/{0}.alts.filter.2gemini.vcf'.format(c['qid']))
elif control_txt:
    bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} !(($1$2$4$5) in vals)' \
              ../query/data/{0}/{1} ../query/output/{0}/{0}.alts.filter.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.2gemini.vcf""".format(c['qid'], c['vcf_filename']))
elif control_vcf:
    # Merge with additional control vcf
    # Slim control vcf only keeping variants occur in the main vcf (to reduce runnig time on processing big VCF file)
    bash_pipe("""bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ../query/output/{0}/{0}.alts.filter.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.recode.pos.txt""".format(c['qid']))
    bash_pipe("""grep '^#' ../query/output/{0}/{0}.ctrl.alts.vcf > ../query/output/{0}/{0}.ctrl.alts.2merge.vcf""".format(c['qid']))
    bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.alts.filter.recode.pos.txt ../query/output/{0}/{0}.ctrl.alts.vcf \
              >> ../query/output/{0}/{0}.ctrl.alts.2merge.vcf""".format(c['qid'], c['vcf_filename']))
    bash_pipe("""bgzip ../query/output/{0}/{0}.alts.filter.eff.vcf""".format(c['qid']))
    bash_pipe("""tabix ../query/output/{0}/{0}.alts.filter.eff.vcf.gz""".format(c['qid']))
    bash_pipe("""bgzip ../query/output/{0}/{0}.ctrl.alts.2merge.vcf""".format(c['qid']))
    bash_pipe("""tabix ../query/output/{0}/{0}.ctrl.alts.2merge.vcf.gz""".format(c['qid']))
    bash_pipe("""bcftools merge ../query/output/{0}/{0}.alts.filter.eff.vcf.gz ../query/output/{0}/{0}.ctrl.alts.2merge.vcf.gz \
              > ../query/output/{0}/{0}.alts.filter.2gemini.vcf""".format(c['qid']))

#  Map hg38 position to hg19 using (positions extracted from dbNSFP -- compound heterozygous model only focuses on nonsynonymous coding mutaitns)
if c['build'] == 'GRCh38.79':
    os.rename('../query/output/{0}/{0}.alts.filter.2gemini.vcf ../query/output/{0}/{0}.alts.filter.2gemini.2hg38.vcf'.format(c['qid']))
    fh_2hg38 = open("../query/output/{0}/{0}.alts.filter.2gemini.2hg38.vcf".format(c['qid'])).readlines()
    fh_hg38 = open("../query/output/{0}/{0}.alts.filter.2gemini.vcf".format(c['qid']), 'w')
    for line in fh_2hg38:
        if line[0] == "#":
            fh_hg38.write(line)
        else:
            content = line.strip().split("\t")
            hg38_pos = tuple(content[:2])
            cur.execute('select hg_19 from hg38_hg19_mapping where (CHROM,hg_38)= (%s, %s);', hg38_pos)
            rows = cur.fetchall()
            if rows == ():
                continue
            else:
                content[1] == rows[0][0]
            fh_hg38.write("\t".join(content) + '\n')
    fh_hg38.close()

# Create PED file of the merged VCF for GEMINI
samples_all = sub.Popen(["""bcftools query -l ../query/output/{0}/{0}.alts.filter.2gemini.vcf.gz""".format(c['qid'])], shell=True, stdout=sub.PIPE, stderr=sub.PIPE).communicate()[0][:-1].split('\n')
samples_ped = individual_info.keys()
comphets_ped = [line.split('\t')[:6] for line in open('../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n')]
samples_ctrl = [i for i in samples_all if i not in samples_ped and i != '']
# Controls: set phenotype = 1 in ped
if samples_ctrl != []:
    for indv in samples_ctrl:
        comphets_ped.append([indv, indv, '0', '0', '0', '1'])
of_comphets_ped = open('../query/output/{0}/{0}.comphets.ped'.format(c['qid']), 'w')
for i in comphets_ped:
    of_comphets_ped.write('\t'.join(i) + '\n')
of_comphets_ped.close()

# Run gemini
bash_pipe("python helper_scripts/convert_hap2dip.py \
          -i ../query/output/{0}/{0}.alts.filter.2gemini.vcf \
          -o ../query/output/{0}/{0}.alts.filter.2gemini.vcf".format(c['qid']))
bash_pipe('bgzip ../query/output/{0}/{0}.alts.filter.2gemini.vcf'.format(c['qid']))
bash_pipe('tabix ../query/output/{0}/{0}.alts.filter.2gemini.vcf.gz'.format(c['qid']))
bash_pipe("""gemini load --cores 4 -v ../query/output/{0}/{0}.alts.filter.2gemini.vcf.gz \
          -p ../query/output/{0}/{0}.comphets.ped \
          -t snpEff --skip-gerp-bp --skip-cadd --skip-gene-tables \
          ../query/output/{0}/{0}.alts.filter.comphets.db""".format(c['qid']))
bash_pipe('''gemini comp_hets --max-priority 1 \
          --columns "chrom,start,rs_ids,ref,alt,filter,impact_so,impact_severity,gene,transcript,biotype,codon_change,aa_change,\
          aaf_adj_exac_all,aaf_adj_exac_afr,aaf_adj_exac_amr,aaf_adj_exac_eas,aaf_adj_exac_fin,aaf_adj_exac_nfe,aaf_adj_exac_sas,aaf_adj_exac_oth,\
          aaf_esp_all,aaf_esp_ea,aaf_esp_aa,\
          aaf_1kg_all,aaf_1kg_amr,aaf_1kg_afr,aaf_1kg_eas,aaf_1kg_eur,aaf_1kg_sas" \
          ../query/output/{0}/{0}.alts.filter.comphets.db \
          > ../query/output/{0}/{0}.alts.filter.eff.spc.cs.txt'''.format(c['qid']))
