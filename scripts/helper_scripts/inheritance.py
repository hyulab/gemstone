# Title: inheritance.py
# Description: This script is to run co-segragation analysis with any one of the following inheritance models:
#               Autosomal dominant
#               Autosomal recessive homozygous
#               X-linked dominant
#               X-linked recessive
#               Y-linked
#               OR
#               No inheritance model screening variants with no co-segregating restraints
# Author: Siwei Chen
# Date created: 18/08/2015
# Date last modified: 29/02/2016
# Python version: 2.7.5


import argparse
import os
import imp


parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
parser.add_argument("-f", help="family id")
args = parser.parse_args()

c = imp.load_source('cache', args.c).cache

pedigrees = dict([(line.split('\t')[1], line.split('\t')[2:])
                  for line in open('../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n')
                  if line.split('\t')[0] == args.f])

# For each inheritance model,
# use CRHOM,POS,REF,ALT as unique identifier for each variant,
# use VCFtools to select/remove samples,
# use BCFtools to extract genotypes.

if c['inheritance'] == 'Dominant':
    # Categorize sample(s) within each family from PED file
    # for differnt requirements on genotypes
    # Unaffected sample(s)
    indv_homo = [individual for individual in pedigrees if pedigrees[individual][3] == '1']
    # Affected male sample(s)
    indv_nonpar = [individual for individual in pedigrees if pedigrees[individual][3] == '2' and pedigrees[individual][2] == '1']
    # Affected samples are required to be heterozygous
    # Remove homozygous variants, except for chromosome X non-pseudoautosomal regions (PAR) in male samples.
    if c['nonPAR'] != 'Y' and len(indv_nonpar) == len(pedigrees):
        open('../query/output/{0}/{0}.fm_nonpar.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_nonpar) + '\n')
        os.system(	"""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """			.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_nonpar.txt --recode --stdout | """			.format(c['qid'], c['ped_filename']) + 
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | sed '/1[\|\/]1/d' |"""
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.nonpar.txt"""						.format(c['qid'], c['vcf_filename']))
        os.system(	"""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """			.format(c['qid'], c['vcf_filename']) +
                   """--chr X --chr Y """ + 
                   """--exclude-positions ../../GeMSTONE-data/PAR.pos.{0}.txt """.format(c['build']) + 
                   """--keep ../query/output/{0}/{0}.fm_nonpar.txt --recode --stdout | """			.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | """
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' """ +
                   """>> ../query/output/{0}/{0}.cs.pos.heter.nonpar.txt"""						.format(c['qid'], c['vcf_filename']))
        os.system("""sort ../query/output/{0}/{0}.cs.pos.heter.nonpar.txt | \
                  uniq > ../query/output/{0}/{0}.cs.pos.heter.txt""".format(c['qid'], c['vcf_filename']))   
    else:
        indv_heter = [individual for individual in pedigrees if pedigrees[individual][3] == '2']        
        open('../query/output/{0}/{0}.fm_heter.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_heter) + '\n')    
        indv_aff = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
        open('../query/output/{0}/{0}.fm_aff.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_aff) + '\n')    
        os.system(	"""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """			.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_heter.txt --recode --stdout | """			.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | sed '/1[\|\/]1/d' |"""
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.txt"""						.format(c['qid'], c['vcf_filename']))
    # Unaffected family members (if there are any) are required to be homozygous wildtype
    if  indv_homo != []: 
        open('../query/output/{0}/{0}.fm_homo.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """		.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_homo.txt --recode --stdout | """    	.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/][\.01]/d' | sed '/[\.0][\|\/]1/d' | sed '/\.[\|\/]\./d' | """ + 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """
                   """> ../query/output/{0}/{0}.cs.pos.homo.txt"""			                .format(c['qid'], c['vcf_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.heter.txt ../query/output/{0}/{0}.cs.pos.homo.txt \
                  > ../query/output/{0}/{0}.cs.pos.txt"""                                           .format(c['qid'], c['vcf_filename']))
    else:
        os.system("""mv ../query/output/{0}/{0}.cs.pos.heter.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    # Use variant identifilers from position file to filter VCF file
    os.system("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.alts.filter.spc.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt""".format(c['qid'], c['vcf_filename']))


elif c['inheritance'] == 'Recessive':
    # Affected samples are required to be homozygous mutant type
    # Unaffected parents (if there are any) of affected child are required to be heterozygous
    indv_homo = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
    indv_heter = list(set([pedigrees[pedigrees[individual][0]] for individual in indv_homo if pedigrees[individual][0] in pedigrees and pedigrees[pedigrees[individual][0]][3] == '1'] +
                          [pedigrees[pedigrees[individual][1]] for individual in indv_homo if pedigrees[individual][1] in pedigrees and pedigrees[pedigrees[individual][1]][3] == '1']))
    indv_oth = [i for i in pedigrees if i not in indv_homo and i not in indv_heter]
    open('../query/output/{0}/{0}.fm_aff.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo) + '\n')
    os.system(	"""vcftools """
               """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """	.format(c['qid'], c['vcf_filename']) +
               """--keep ../query/output/{0}/{0}.fm_aff.txt --recode --stdout | """	.format(c['qid'], c['ped_filename']) +
               """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.01]/d' | sed '/1[\|\/][\.0]/d' | sed '/\.[\|\/][\.01]/d' | """ 
               """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
               """> ../query/output/{0}/{0}.cs.pos.txt"""			.format(c['qid'], c['vcf_filename']))
    # Keep common positions from the intersection of different categories
    if indv_heter != []:
        open('../query/output/{0}/{0}.fm_heter.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_heter) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """	.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_heter.txt --recode --stdout | """	.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/]1/d' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.txt"""			.format(c['qid'], c['vcf_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.heter.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    if indv_oth != []:
        open('../query/output/{0}/{0}.fm_oth.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_oth) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """	.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_oth.txt --recode --stdout | """	.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/]1/d' | sed '/\.[\|\/]\./d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.oth.txt"""			.format(c['qid'], c['vcf_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.oth.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    
    os.system("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.alts.filter.spc.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt""".format(c['qid'], c['vcf_filename']))


elif c['inheritance'] == 'XR':
    # Affected M and F = homo1
    # Unaffected F whose father is affected = heter
    # Unaffected M = homo0
    # Other unaffected = homo0, heter
    indv_homo1 = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
    indv_heter = [individual for individual in pedigrees if pedigrees[individual][2] == '2' and pedigrees[individual][3] == '1' and pedigrees[individual][0] in indv_homo1]
    indv_homo0 = [individual for individual in pedigrees if pedigrees[individual][2] == '1' and pedigrees[individual][3] == '1']
    indv_oth = [individual for individual in set([individual for individual in pedigrees if pedigrees[individual][3] == '1']) - set(indv_heter + indv_homo0)]
    open('../query/output/{0}/{0}.fm_homo1.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo1) + '\n')
    indv_aff = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
    open('../query/output/{0}/{0}.fm_aff.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_aff) + '\n')
    os.system(	"""vcftools """
               """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
               """--keep ../query/output/{0}/{0}.fm_homo1.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
               """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.01]/d' | sed '/1[\|\/][\.0]/d' | sed '/\.[\|\/][\.01]/d' | """ 
               """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
               """> ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    if indv_heter != []:
        open('../query/output/{0}/{0}.fm_heter.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_heter) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_heter.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/]1/d' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.heter.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    if indv_homo0 != []:
        open('../query/output/{0}/{0}.fm_homo0.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo0) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_homo0.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/][\.10]/d' | sed '/0[\|\/]1/d' | sed '/\.[\|\/][\.1]/d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.homo0.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.homo0.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))        
    if indv_oth != []:
        open('../query/output/{0}/{0}.fm_oth.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_oth) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_oth.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/]1/d' | sed '/\.[\|\/]\./d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.oth.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.oth.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
        
    os.system("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.alts.filter.spc.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt""".format(c['qid'], c['vcf_filename']))

elif c['inheritance'] == 'XD':
    # Affected M = homo1
    # Affected F with one and only one of her parents is affected = heter
    # Other affected F = heter + homo1
    # Unaffected M and F = homo0
    indv_aff = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
    open('../query/output/{0}/{0}.fm_aff.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_aff) + '\n')
    indv_homo1 = [individual for individual in pedigrees if pedigrees[individual][2] == '1' and pedigrees[individual][3] == '2']
    indv_heter = [individual for individual in pedigrees if pedigrees[individual][2] == '2' and pedigrees[individual][3] == '2' and set(
        pedigrees[individual][0:2]) < set(pedigrees) and pedigrees[pedigrees[individual][0]][3] == '2' and pedigrees[pedigrees[individual][1]][3] == '1']
    indv_heter += [individual for individual in pedigrees if pedigrees[individual][2] == '2' and pedigrees[individual][3] ==
                   '2' and set(pedigrees[individual][0:2]) < set(pedigrees) and pedigrees[pedigrees[individual][0]][3] == '1' and pedigrees[pedigrees[individual][1]][3] == '2']
    indv_oth = [individual for individual in set([individual for individual in pedigrees if pedigrees[individual][3] == '2']) - set(indv_homo1 + indv_heter)]
    indv_homo0 = [individual for individual in pedigrees if pedigrees[individual][3] == '1']      
    open('../query/output/{0}/{0}.fm_homo1.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo1) + '\n')
    os.system(	"""vcftools """
               """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
               """--keep ../query/output/{0}/{0}.fm_homo1.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
               """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.01]/d' | sed '/1[\|\/][\.0]/d' | sed '/\.[\|\/][\.01]/d' | """ 
               """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
               """> ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    if indv_heter != []:
        open('../query/output/{0}/{0}.fm_heter.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_heter) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_heter.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/]1/d' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.heter.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    if indv_homo0 != []:
        open('../query/output/{0}/{0}.fm_homo0.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo0) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_homo0.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/][\.10]/d' | sed '/0[\|\/]1/d' | sed '/\.[\|\/][\.1]/d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.homo0.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.homo0.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))        
    if indv_oth != []:
        open('../query/output/{0}/{0}.fm_oth.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_oth) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_oth.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/]1/d' | sed '/\.[\|\/]\./d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.oth.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.oth.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))        
    os.system("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.alts.filter.spc.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt""".format(c['qid'], c['vcf_filename']))
    
    
elif c['inheritance'] == 'Y':
    # Affected M and their fathers and sons= homo1
    # Unaffected M = homo0
    indv_aff = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
    open('../query/output/{0}/{0}.fm_aff.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_aff) + '\n')    
    indv_homo1 = [individual for individual in pedigrees if pedigrees[individual][2] == '1' and pedigrees[individual][3] == '2']
    indv_homo1 += [individual for individual in indv_homo1 if pedigrees[individual][0] in pedigrees and pedigrees[pedigrees[individual][0]][3] == '2']
    indv_homo1 += [individual for individual in pedigrees if pedigrees[individual][0] in indv_homo1]
    indv_homo1 = [individual for individual in set(indv_homo1)]
    indv_homo0 = [individual for individual in pedigrees if pedigrees[individual][2] == '1' and pedigrees[individual][3] == '1']

    open('../query/output/{0}/{0}.fm_homo1.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo1) + '\n')
    os.system(	"""vcftools """
               """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
               """--keep ../query/output/{0}/{0}.fm_homo1.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
               """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.01]/d' | sed '/1[\|\/][\.0]/d' | sed '/\.[\|\/][\.01]/d' | """ 
               """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
               """> ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    if indv_homo0 != []:
        open('../query/output/{0}/{0}.fm_homo0.txt'.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo0) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_homo0.txt --recode --stdout | """.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/][\.10]/d' | sed '/0[\|\/]1/d' | sed '/\.[\|\/][\.1]/d' | """ 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.homo0.txt""".format(c['qid'], c['ped_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.cs.pos.homo0.txt | sort """.format(c['qid'], c['vcf_filename']) +
                   """> ../query/output/{0}/{0}.cs.pos.comm.txt"""	.format(c['qid'], c['vcf_filename']))
        os.system("""mv ../query/output/{0}/{0}.cs.pos.comm.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    
    os.system("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.alts.filter.spc.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt""".format(c['qid'], c['vcf_filename']))
    
elif c['inheritance'] == 'No':    
    indv_homo = [individual for individual in pedigrees if pedigrees[individual][3] == '1']
    indv_nonpar = [individual for individual in pedigrees if pedigrees[individual][3] == '2' and pedigrees[individual][2] == '1']
    if c['nonPAR'] == 'Y' and len(indv_nonpar) == len(pedigrees):
        open('../query/output/{0}/{0}.fm_nonpar.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_nonpar) + '\n')
        os.system(	"""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """			.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_nonpar.txt --recode --stdout | """			.format(c['qid'], c['ped_filename']) + 
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' |"""
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.nonpar.txt"""						.format(c['qid'], c['vcf_filename']))
        os.system(	"""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """			.format(c['qid'], c['vcf_filename']) +
                   """--chr X --chr Y """ + 
                   """--exclude-positions ../../GeMSTONE-data/PAR.pos.{0}.txt """.format(c['build']) + 
                   """--keep ../query/output/{0}/{0}.fm_nonpar.txt --recode --stdout | """			.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' | """
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' """ +
                   """>> ../query/output/{0}/{0}.cs.pos.heter.nonpar.txt"""						.format(c['qid'], c['vcf_filename']))
        os.system("""sort ../query/output/{0}/{0}.cs.pos.heter.nonpar.txt | \
                  uniq > ../query/output/{0}/{0}.cs.pos.heter.txt""".format(c['qid'], c['vcf_filename']))    
    else:
        indv_heter = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
        
        open('../query/output/{0}/{0}.fm_heter.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_heter) + '\n')
    
        indv_aff = [individual for individual in pedigrees if pedigrees[individual][3] == '2']
        open('../query/output/{0}/{0}.fm_aff.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_aff) + '\n')
    
        os.system(	"""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """			.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_heter.txt --recode --stdout | """			.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/0[\|\/][\.0]/d' | sed '/\.[\|\/][\.0]/d' |"""
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """ +
                   """> ../query/output/{0}/{0}.cs.pos.heter.txt"""						.format(c['qid'], c['vcf_filename']))         
    if  indv_homo != []: 
        open('../query/output/{0}/{0}.fm_homo.txt'	.format(c['qid'], c['ped_filename'])	, 'w').write('\n'.join(indv_homo) + '\n')
        os.system("""vcftools """
                   """--vcf ../query/output/{0}/{0}.alts.filter.spc.eff.vcf """		.format(c['qid'], c['vcf_filename']) +
                   """--keep ../query/output/{0}/{0}.fm_homo.txt --recode --stdout | """    	.format(c['qid'], c['ped_filename']) +
                   """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed '/1[\|\/][\.01]/d' | sed '/[\.0][\|\/]1/d' | sed '/\.[\|\/]\./d' | """ + 
                   """awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort """
                   """> ../query/output/{0}/{0}.cs.pos.homo.txt"""			                .format(c['qid'], c['vcf_filename']))
        os.system("""comm -12 ../query/output/{0}/{0}.cs.pos.heter.txt ../query/output/{0}/{0}.cs.pos.homo.txt \
                  > ../query/output/{0}/{0}.cs.pos.txt"""                                           .format(c['qid'], c['vcf_filename']))
    else:
        os.system("""mv ../query/output/{0}/{0}.cs.pos.heter.txt ../query/output/{0}/{0}.cs.pos.txt""".format(c['qid'], c['vcf_filename']))
    
    os.system("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/{0}.cs.pos.txt ../query/output/{0}/{0}.alts.filter.spc.eff.vcf \
              > ../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt""".format(c['qid'], c['vcf_filename']))

