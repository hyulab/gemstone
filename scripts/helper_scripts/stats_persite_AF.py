# Title: stats_persite_AF.py
# Description: This script is to calculate statistics (of variant density, variant quality, mean depth, allele frequency, ts/tv ratio, insertion/deletion lengths, variant type) on filtered VCF for visualization.
# Author: Siwei Chen
# Date created: 14/10/2015
# Date last modified: 26/04/2016
# Python version: 2.7.5

import argparse
import os
import imp
import MySQLdb

mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

parser = argparse.ArgumentParser()
parser.add_argument("-v", help="vcf file")
args = parser.parse_args()

vcf_file = args.v
gendir = '../query/D3_output/'+vcf_file.split('/')[3]+'/filtered/general'

chrom= ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

for c in chrom:
    if not os.path.exists("""{0}/{1}""".format(gendir,c)):
	os.mkdir("""{0}/{1}""".format(gendir,c))


# Define a function that takes in data as a list and generates input for plotting a histogram
def listToHistogram(L, b_bins):
    import numpy as np
    return '\n'.join(['{0},{1}'.format(x, height) for height,x in sorted(zip(*np.histogram(L,bins=b_bins)))])

# Chromosome length rounded to 10k
chr_pos4= {'20': '2718570000',
           '21': '2781600000',
           '22': '2829730000',
           '1': '0',
           '3': '492450000',
           '2': '249250000',
           '5': '881620000',
           '4': '690470000',
           '7': '1233660000',
           '6': '1062540000',
           '9': '1539160000',
           '8': '1392800000',
           '11': '1815900000',
           '10': '1680370000',
           '13': '2084760000',
           '12': '1950910000',
           '15': '2307280000',
           '14': '2199930000',
           '17': '2500160000',
           '16': '2409810000',
           '19': '2659440000',
           '18': '2581360000',
           'Y': '3036300000',
           'X': '2881030000',
           'M': '3095670000'}

# Variant density
os.system('''vcftools --vcf {0} --SNPdensity 1000000 --out {1}/tmp'''.format(vcf_file,gendir))
of_density= open('{0}/variant_density_1m.csv'.format(gendir),'w')
of_density.write('BIN_START,SNP_COUNT,ABSOLUTE_POS_4\n')
if len(open('{0}/tmp.snpden'.format(gendir)).readlines()) > 1:
    for i in open('{0}/tmp.snpden'.format(gendir)).readlines()[1:]:
	of_density.write(','.join([str(int(i[:-1].split('\t')[1])/10000),i[:-1].split('\t')[2],str((int(i[:-1].split('\t')[1])+int(chr_pos4[i[:-1].split('\t')[0]]))/10000)]) + '\n')
of_density.close()
THRESHOLD = 100
SOFTNESS  = 25
variant_density = [line.split(',') for line in open('{0}/variant_density_1m.csv'.format(gendir)).read().strip().split('\n')]
header = ','.join(variant_density[0])
open('{0}/variant_density.csv'.format(gendir),'w').write(header+'\n')
if len(variant_density) > 1:
    data = [[int(x) for x in line] for line in variant_density[1:]]
    for index in range(len(data)-1):
	    start, snp, absolute    = data[index]
	    start2, snp2, absolute2 = data[index+1]
	    if abs(absolute - absolute2) > THRESHOLD:
		    data.append([start +SOFTNESS, 0, absolute +SOFTNESS])
		    data.append([start2-SOFTNESS, 0, absolute2-SOFTNESS])
		    #print 'find gap!'
    open('{0}/variant_density.csv'.format(gendir),'w').write(header+'\n'+'\n'.join([','.join([str(x) for x in line]) for line in sorted(data,key=lambda x:x[-1])]))
os.system("""rm {0}/variant_density_1m.csv""".format(gendir))


# Variant quality
os.system('''vcftools --vcf {0} --site-quality --out {1}/tmp'''.format(vcf_file,gendir))
quals= [float(i.split('\t')[2]) for i in open('{0}/tmp.lqual'.format(gendir)).readlines()[1:] if len(open('{0}/tmp.lqual'.format(gendir)).readlines()) > 1 and 'nan' not in i.split('\t')[2]]
open('{0}/variant_quality_hist.csv'.format(gendir),'w').write('POSITION,HEIGHT\n' + listToHistogram(quals,100))
# Generate data for each chromosome
for c in chrom:
    quals= [float(i.split('\t')[2]) for i in open('{0}/tmp.lqual'.format(gendir)).readlines()[1:] if len(open('{0}/tmp.lqual'.format(gendir)).readlines()) > 1 and i.split('\t')[0] == c and 'nan' not in i.split('\t')[2]]
    open('{0}/{1}/variant_quality_hist.csv'.format(gendir,c),'w').write('POSITION,HEIGHT\n' + listToHistogram(quals,100))
os.system("""rm {0}/tmp.lqual""".format(gendir))    

# Extract positions of varaints of interest  
os.system('''sed '/^#/d' {0} | awk 'BEGIN {{FS= "\t"; OFS= "\t"}} {{print $1,$2,$4,$5}}' > {1}/tmp_pos.txt'''.format(vcf_file,gendir))

# Define a function takes in a chromosome number (1,2 instead of chr1,chr2) and return a list of MAF (in ExAC database) of variants on this chromosome 
def afSpectrum(chromosome):
    af= []
    for var in [line for line in open('{0}/tmp_pos.txt'.format(gendir)).readlines() if line.strip().split("\t")[0] == chromosome]:
	query= tuple(var[:-1].split('\t'))
	cur.execute('select AF_Adj from ExAC_hg19 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);',query)
	rows= cur.fetchall()
	if rows == ():
	    #noval += 1
	    af.append(0.0)
	elif round(rows[0][0],3) <= 50.0:
	    af.append(round(rows[0][0],3))
	else:
	    af.append(100.0 - round(rows[0][0],3))
    open('{0}/{1}/af_spectrum_hist.csv'.format(gendir,chromosome),'w').write('POSITION,HEIGHT\n' + listToHistogram(af,100))
    return af

# Compile chromosomes data to generate the whole spectrum
whole_spectrum = []
for c in chrom:
    whole_spectrum += afSpectrum(c)
open('{0}/af_spectrum_hist.csv'.format(gendir),'w').write('POSITION,HEIGHT\n' + listToHistogram(whole_spectrum,100))


# DP: the mean depth per site averaged across all individuals
os.system('''vcftools --vcf {0} --site-mean-depth --out {1}/tmp'''.format(vcf_file,gendir))
open('{0}/mean_depth_hist.csv'.format(gendir),'w').write('POSITION,HEIGHT\n' + listToHistogram([float(i[:-1].split('\t')[2]) for i in open('{0}/tmp.ldepth.mean'.format(gendir)).readlines()[1:] if len(open('{0}/tmp.ldepth.mean'.format(gendir)).readlines()) > 1 and 'nan' not in i.split('\t')[2]],100))
# Generate data for each chromosome
for c in chrom:
    depths= [float(i.split('\t')[2]) for i in open('{0}/tmp.ldepth.mean'.format(gendir)).readlines()[1:] if len(open('{0}/tmp.ldepth.mean'.format(gendir)).readlines()) > 1 and i.split('\t')[0] == c and 'nan' not in i.split('\t')[2]]
    open('{0}/{1}/mean_depth_hist.csv'.format(gendir,c),'w').write('POSITION,HEIGHT\n' + listToHistogram(depths,100))
os.system("""rm {0}/tmp.ldepth.mean""".format(gendir)) 

# TsTv ratio
os.system('''vcftools --vcf {0} --FILTER-summary --out {1}/tmp'''.format(vcf_file,gendir))
os.system('''sed 's/\t/,/g' {0}/tmp.FILTER.summary > {0}/tstv_ratio.csv'''.format(gendir))
# Generate data for each chromosome
for c in chrom:
    os.system('''vcftools --vcf {0} --chr {1} --FILTER-summary --out {2}/{1}/tmp'''.format(vcf_file,c,gendir))
    os.system('''sed 's/\t/,/g' {0}/{1}/tmp.FILTER.summary > {0}/{1}/tstv_ratio.csv'''.format(gendir,c))
    os.system('''rm {0}/{1}/tmp.FILTER.summary'''.format(gendir,c))


# Variant Types
# INDEL Lengths
os.system('''vcftools --vcf {0} --hist-indel-len --out {1}/tmp'''.format(vcf_file,gendir))
SNP= 0
Ins= 0
Del= 0
All= [pos for pos in open('{0}/tmp_pos.txt'.format(gendir)).readlines()]
Ttl= len(All) if len(All) != 0 else 1
for indel in open('{0}/tmp.indel.hist'.format(gendir)).readlines()[1:]:
    if int(indel.split('\t')[0]) < 0:
        Del += int(indel.split('\t')[1])
    if int(indel.split('\t')[0]) == 0:
        SNP += int(indel.split('\t')[1])
    if int(indel.split('\t')[0]) > 0:
        Ins += int(indel.split('\t')[1])
Oth= Ttl - SNP - Ins - Del
SNPPRCT= round(float(SNP)/Ttl*100, 2) if Ttl != 0 else 0.0
InsPRCT= round(float(Ins)/Ttl*100, 2) if Ttl != 0 else 0.0
DelPRCT= round(float(Del)/Ttl*100, 2) if Ttl != 0 else 0.0
OthPRCT= round(float(Oth)/Ttl*100, 2) if Ttl != 0 else 0.0
of_type= open('{0}/variant_type.csv'.format(gendir),'w')
of_type.write('TYPE,COUNT,PRCT\n')
of_type.write('\n'.join(['SNP,'+str(SNP)+','+str(SNPPRCT)+'%','Ins,'+str(Ins)+','+str(InsPRCT)+'%','Del,'+str(Del)+','+str(DelPRCT)+'%','Other,'+str(Oth)+','+str(OthPRCT)+'%']) + '\n')
of_type.close()
os.system('''awk 'BEGIN {{FS= "\t"; OFS= "\t"}};{{print $1,$2}}' {0}/tmp.indel.hist | sed 's/\t/,/g' > {0}/indel_length.csv'''.format(gendir))
# Generate data for each chromosome
for c in chrom:
    os.system('''vcftools --vcf {0} --chr {1} --hist-indel-len --out {2}/{1}/tmp'''.format(vcf_file,c,gendir))
    SNP= 0
    Ins= 0
    Del= 0
    All= [pos for pos in open('{0}/tmp_pos.txt'.format(gendir)).readlines() if pos.split('\t')[0] == c]
    Ttl= len(All) if len(All) != 0 else 1
    for indel in open('{0}/{1}/tmp.indel.hist'.format(gendir,c)).readlines()[1:]:
	if int(indel.split('\t')[0]) < 0:
	    Del += int(indel.split('\t')[1])
	if int(indel.split('\t')[0]) == 0:
	    SNP += int(indel.split('\t')[1])
	if int(indel.split('\t')[0]) > 0:
	    Ins += int(indel.split('\t')[1])
    Oth= Ttl - SNP - Ins - Del
    SNPPRCT= round(float(SNP)/Ttl*100, 2) if Ttl != 0 else 0.0
    InsPRCT= round(float(Ins)/Ttl*100, 2) if Ttl != 0 else 0.0
    DelPRCT= round(float(Del)/Ttl*100, 2) if Ttl != 0 else 0.0
    OthPRCT= round(float(Oth)/Ttl*100, 2) if Ttl != 0 else 0.0
    of_type= open('{0}/{1}/variant_type.csv'.format(gendir,c),'w')
    of_type.write('TYPE,COUNT,PRCT\n')
    of_type.write('\n'.join(['SNP,'+str(SNP)+','+str(SNPPRCT)+'%','Ins,'+str(Ins)+','+str(InsPRCT)+'%','Del,'+str(Del)+','+str(DelPRCT)+'%','Other,'+str(Oth)+','+str(OthPRCT)+'%']) + '\n')
    of_type.close()
    os.system('''awk 'BEGIN {{FS= "\t"; OFS= "\t"}};{{print $1,$2}}' {0}/{1}/tmp.indel.hist | sed 's/\t/,/g' > {0}/{1}/indel_length.csv'''.format(gendir,c))
    os.system("""rm {0}/{1}/tmp.*""".format(gendir,c))
os.system("""rm {0}/tmp*""".format(gendir))

