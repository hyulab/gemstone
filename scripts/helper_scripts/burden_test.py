# Title: burden_test.py
# Description: This script is to run burden tests using PLINK/SEQ package.
# Author: Siwei Chen
# Date created: 15/09/2015
# Date last modified: 09/05/2016
# Python version: 2.7.5

from glob import glob
import time
from multiprocessing.pool import ThreadPool
import os
import subprocess as sub
import imp
from itertools import groupby
import shutil
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-c", help="parameter file")
args = parser.parse_args()

filename= args.c
c = imp.load_source('cache', filename).cache

os.system("/usr/local/bin/pseq ../query/output/{0}/burden/{0} new-project --resources ../../GeMSTONE-data/pseq/hg19".format(c['qid']))

# Cases: unrelated samples only (first to be selected in families)
# Controls: ONLY from separate control vcf to make sure there are no related samples
# Pheno cases 2, control 1

if os.path.isfile("../query/output/{0}/{0}.alts.filter.recode.vcf.gz".format(c['qid'])):
    os.system("gunzip ../query/output/{0}/{0}.alts.filter.recode.vcf.gz".format(c['qid']))

# Parse PED file and store pedigree structures
ped_contents = [line.split('\t') for line in open(
    '../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n') if line.split('\t')[5] == '2'] # only include affected indv in each family
aff_family_members = dict([(key, list([g[1] for g in group])) for key, group in groupby([(p[0], p[1]) for p in ped_contents], lambda x: x[0])]) # not sorted keep order
# Extract the first member in each family as proband
cases= [i[0] for i in aff_family_members.values()]
# Extract controls in VCF file
controls = sub.Popen(["""bcftools query -l ../query/data/{0}/{1}""".format(c['qid'],c['burden_control'])],shell = True,stdout=sub.PIPE,stderr=sub.PIPE).communicate()[0].strip().split("\n")

if c['burden_control'][-3:] == '.gz':
    os.system("""gunzip ../query/data/{0}/{1}""".format(c['qid'],c['burden_control']))
    c['burden_control'] = c['burden_control'][:-3]

if c['build'] == 'GRCh38.79':
    os.system("""bcftools convert --haploid2diploid ../query/data/{0}/{1} | \
              sed 's/ID=AD,Number=./ID=AD,Number=R/' | \
              sed 's/chr//g' | \
              vt decompose -s - | \
              vt decompose_blocksub -a - | \
              vt normalize -n -r ../../GeMSTONE-data/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa - | \
              vt uniq - > ../query/output/{0}/{0}.ctrl_burden.alts.vcf """.format(c['qid'], c['burden_control']))
if c['build'] == 'GRCh37.75':
    os.system("""bcftools convert --haploid2diploid ../query/data/{0}/{1} |\
              sed 's/ID=AD,Number=./ID=AD,Number=R/' | \
              sed 's/chr//g' | \
              vt decompose -s - | \
              vt decompose_blocksub -a - | \
              vt normalize -n -r ../../GeMSTONE-data/human_g1k_v37.fasta - | \
              vt uniq - > ../query/output/{0}/{0}.ctrl_burden.alts.vcf """.format(c['qid'], c['burden_control']))

# QC on control VCF
os.system("""vcftools --vcf ../query/output/{0}/{0}.ctrl_burden.alts.vcf --minQ {1} {2} --minGQ {3} --minDP {4} """
          .format(c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
          """--recode --out ../query/output/{0}/{0}.ctrl_burden.alts.filter""".format(c['qid']))

# Index VCFs
# liftOver from hg38 (if the input is of hg38) to hg19 for PLINK/SEQ input
# Use CorssMap (http://crossmap.sourceforge.net/#convert-vcf-format-files)
if c['build'] == 'GRCh38.79':
    os.rename('../query/output/{0}/{0}.alts.filter.recode.vcf ../query/output/{0}/{0}.alts.filter.recode.2hg19.vcf'.format(c['qid']))
    bash_pipe("python /bin/CrossMap.py \
              vcf ../../GeMSTONE-data/CrossMap/hg38ToHg19.over.chain.gz \
              ../query/output/{0}/{0}.alts.filter.recode.2hg19.vcf \
              ../../GeMSTONE-data/hg19.fa \
              ../query/output/{0}/{0}.alts.filter.recode.vcf".format(c['qid']))
    os.rename('../query/output/{0}/{0}.ctrl_burden.alts.filter.recode.vcf ../query/output/{0}/{0}.ctrl_burden.alts.filter.recode.hg19.vcf'.format(c['qid']))
    bash_pipe("python /bin/CrossMap.py \
              vcf ../../GeMSTONE-data/CrossMap/hg38ToHg19.over.chain.gz \
              ../query/output/{0}/{0}.ctrl_burden.alts.filter.recode.2hg19.vcf \
              ../../GeMSTONE-data/hg19.fa \
              ../query/output/{0}/{0}.ctrl_burden.alts.filter.recode.vcf".format(c['qid']))


os.system("bgzip -c ../query/output/{0}/{0}.alts.filter.recode.vcf > ../query/output/{0}/{0}.alts.filter.recode.vcf.gz".format(c['qid']))
os.system("""/usr/local/bin/pseq ../query/output/{0}/burden/{0}.pseq index-vcf --vcf ../query/output/{0}/{0}.alts.filter.recode.vcf.gz""".format(c['qid']))
os.system("bgzip ../query/output/{0}/{0}.ctrl_burden.alts.filter.recode.vcf".format(c['qid']))
os.system("""/usr/local/bin/pseq ../query/output/{0}/burden/{0}.pseq index-vcf --vcf ../query/output/{0}/{0}.ctrl_burden.alts.filter.recode.vcf.gz""".format(c['qid']))


# Create masking indvidual list
open("../query/output/{0}/burden/indv_list.txt".format(c['qid']),'w').write('\n'.join(cases+controls) + '\n')
# Create and load phe file
phe_file=open("../query/output/{0}/burden/{0}.phe".format(c['qid']),'w')
phe_file.write('''##PHENO,Integer,0,"Primary disease phenotype"\n#ID\tPHENO\n''')
for i in cases:
    phe_file.write(i+'\t2\n')
for i in controls:
    phe_file.write(i+'\t1\n')
phe_file.close()
os.system("""/usr/local/bin/pseq ../query/output/{0}/burden/{0}.pseq load-pheno --file ../query/output/{0}/burden/{0}.phe""".format(c['qid']))

# Create masking loci list from the all_ file
os.system("""sed "/^GENE_NAME/d" ../query/output/{0}/Cross-sample_all_genes.txt | \
          awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $4}}' \
          > ../query/output/{0}/burden/gene_list.txt""".format(c['qid']))


# Run burden tests in parallel -- not multithreading as PINK/SEQ database would be blocked,
# instead copy the database for multiple tests.
def runBurden(arg_tuple):
    qid, burden_test = arg_tuple
    os.system("""cp ../query/output/{0}/burden/{0}.pseq ../query/output/{0}/burden/{1}.pseq""".format(qid,burden_test))
    os.system("""/usr/local/bin/pseq ../query/output/{0}/burden/{1}.pseq assoc \
              --mask loc.group=ensembl --mask loc.subset=ensembl,@../query/output/{0}/burden/gene_list.txt \
              indiv=@../query/output/{0}/burden/indv_list.txt \
              --tests {1} --phenotype PHENO \
              > ../query/output/{0}/burden/{1}.txt""".format(qid,burden_test))
    
for t in c['burden']:
    runBurden((c['qid'],t.lower()))


# Create a dict for assigning values after burden tests
transcript_list= [line[:-1] for line in open("../query/output/{0}/burden/gene_list.txt".format(c['qid'])).readlines()]
assoc_pval={}
for i in transcript_list:
    assoc_pval[i] = []
fh=open("../query/output/{0}/Cross-sample_all_genes.txt".format(c['qid'])).readlines() 
of=open("../query/output/{0}/Cross-sample_all_burden.txt".format(c['qid']),'w')
of.write(fh[0].strip() + '\t' + '\t'.join([i.upper() for i in c['burden']]) + '\n')
for test in c['burden']:
    scores= dict([line[:-1].split("\t")[0],line[:-1].split("\t")[5]] for
        line in open("../query/output/{0}/burden/{1}.txt".format(c['qid'],test.lower())).readlines()[1:]
        if len(line[:-1].split("\t")) > 5)
    for t in transcript_list:
        if t in scores:
            assoc_pval[t].append(scores[t] if scores[t] != 'nan' else '.')
        else:
            assoc_pval[t].append('.')
for line in fh[1:]:
    of.write(line.strip()+'\t'+'\t'.join(assoc_pval[line[:-1].split("\t")[3]]) + '\n')
of.close()

# Keep file name unchanged for next piping
os.system('''mv ../query/output/{0}/Cross-sample_all_burden.txt ../query/output/{0}/Cross-sample_all_genes.txt'''.format(c['qid']))
















