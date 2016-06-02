# Title: af_sc_ANN_comphets.py
# Description: This script is to annotate and/or filter variants for allele frequency, variant consequence and location (particularly on output from recessive compound heterozygous inheritance model).
# Author: Siwei Chen
# Date created: 16/03/2016
# Date last modified: 17/04/2016
# Python version: 2.7.5

import argparse
import os
import imp
import MySQLdb
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
args = parser.parse_args()
c = imp.load_source('cache', args.c).cache

mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

# First convert chromosome name from chr1 to 1
os.system("""sed -i 's/^chr//g' ../query/output/{0}/{0}.alts.filter.eff.spc.cs.txt""".format(c['qid'], c['vcf_filename']))

ped_contents = [line.split('\t') for line in open(
    '../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n')]
sample_pop = dict([(vector[1], vector[-1]) for vector in ped_contents])

# Population code for ethnicity-specific AF filter
pop_code = ['esp_AA', 'esp_ALL', 'esp_EA', 'exac_AFR', 'exac_ALL', 'exac_AMR',
            'exac_Adj', 'exac_EAS', 'exac_FIN', 'exac_NFE', 'exac_OTH', 'exac_SAS',
            'kg_ALL','kg_AFR', 'kg_AMR', 'kg_EAS', 'kg_EUR', 'kg_SAS', 'tagc_AJ',
            'exac','esp','kg','tagc']

# Storing by comphets id may NOT satisfy ceratin requirement,
# insdead using tupe([family_id,comp_het_id] with candidate pairs)
comphets_id2rm = [] 

fh_FM = open('../query/output/{0}/{0}.alts.filter.eff.spc.cs.txt'.format(c['qid'], c['vcf_filename'])).readlines()
fh_ANN = open('../query/output/{0}/{0}.alts.filter.eff.spc.cs.ANN.txt'.format(c['qid'], c['vcf_filename']), 'w')
out_header = "CHROM\tPOS\tID\tREF\tALT\tFILTER" + \
             "\tVARIANT_ID\tFAMILY_ID\tFAMILY_MEMBERS\tFAMILY_GENOTYPES\tSAMPLES\tCOMP_HET_ID\tPRIORITY" + \
             "\tCONSEQUENCE\tPUTATIVE_IMPACT" + \
             "\tGENE_NAME\tENTREZ_ID\tENSEMBL_ID\tFEATURE_ID\tTRANSCRIPT_BIOTYPE\tHGVS.c\tHGVS.p\tUNIPROT\tPFAM_ID\tDOMAIN" \
             if c['domains'] != 'N' else \
             "CHROM\tPOS\tID\tREF\tALT\tFILTER" + \
             "\tVARIANT_ID\tFAMILY_ID\tFAMILY_MEMBERS\tFAMILY_GENOTYPES\tSAMPLES\tCOMP_HET_ID\tPRIORITY" + \
             "\tCONSEQUENCE\tPUTATIVE_IMPACT" + \
             "\tGENE_NAME\tENTREZ_ID\tENSEMBL_ID\tFEATURE_ID\tTRANSCRIPT_BIOTYPE\tHGVS.c\tHGVS.p\tUNIPROT"
if 'ExAC' in c['ann_af']: # control coverage for hg38
    out_header += "\tExAC_MEAN_COVERAGE\tExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_SAS\tExAC_OTH" if c['build'] == 'GRCh37.75' else "\tExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_SAS\tExAC_OTH"
if 'ESP6500' in c['ann_af']:
    out_header += "\tESP6500_ALL\tESP6500_EA\tESP6500_AA"
if '1000Genome' in c['ann_af']:
    out_header += "\t1000G_ALL\t1000G_AMR\t1000G_AFR\t1000G_EAS\t1000G_EUR\t1000G_SAS"
if 'TAGC' in c['ann_af']:
    out_header += "\tTAGC_AF\tTAGC_AC"
fh_ANN.write(out_header + '\n')

for line in fh_FM[1:]:
    content = line.strip().split("\t")
    cols6 = content[:6]
    cols6[1] = str(int(cols6[1])+1)
    cols6[2] = cols6[2] if cols6[2] != 'None' else '.'
    cols6[-1] = cols6[-1] if cols6[-1] != 'None' else 'PASS'
    variant_pair = content[30:35] + content[36:]
    comphets_id = tuple([variant_pair[1],variant_pair[5]])

    # Map hg19 position back to hg38 using (positions extracted from dbNSFP -- compound heterozygous model only focuses on nonsynonymous coding mutaitns)
    if c['build'] == 'GRCh38.79':
        hg19_pos = tuple(cols6[:2])
        cur.execute('select hg_38 from hg19_hg38_mapping where (CHROM,hg_38)= (%s, %s);', hg19_pos)
        rows = cur.fetchall()
        if rows == ():
            comphets_id2rm.append(comphets_id)
            continue
        else:
            cols6[1] = rows[0][0]
            
    so = content[6]
    # GEMINI sometimes does not follow the order of so impact annotated by snpeff, so need to further confirm the so is in c['consequences'] list
    if not [i for i in c['consequences'] if i in so]:
        comphets_id2rm.append(comphets_id)
        continue
    impact = content[6:8]
    gname = content[8]
    entrez = '.'
    cur.execute('select ID from gID_gSymbol where Symbol= %s;', gname)
    rows = cur.fetchall()
    if rows != ():
        entrez = rows[0][0]
    else:
        cur.execute('select ID from gSymbol_gID_biomart where Symbol= %s;', gname)
        rows = cur.fetchall()
        if rows != ():
            entrez = [rows[0][0]]
    enst = content[9]
    cur.execute('select ENSEMBL_GENE_ID from enst_ensg where ENSEMBL_TRANSCRIPT_ID= %s;', enst)
    rows = cur.fetchall()
    ensg = rows[0][0] if rows != () else '.'
    biotype = content[10]
    if biotype not in c['transcript_biotype']:
        comphets_id2rm.append(comphets_id)
        continue
    hgvs = content[11:13]
    uniprot = '.'
    if entrez != '.':
        cur.execute('select Uniprot from gID_Uniprot where Gene_ID = %s;', entrez)
        rows = cur.fetchall()
        if rows != ():
            uniprot = rows[0][0]
    
    # Pfam annotate&/filter
    pfam = [] if c['domains'] == 'N' else ['.','.']
    hgvsp = hgvs[1]
    if (c['domains'] == 'All' or c['ProDomains'] != ['None']) and len(hgvsp) > 5 and hgvsp[0:2] == 'p.' and uniprot != '.':
        aa_pos= int("".join(itertools.takewhile(str.isdigit, hgvsp[5:])))
        cur.execute('''select PFAM,DOMAIN_NAME from UNIPROT_PFAM_DOMAINS where UNIPROT= '%s' and START <= '%d' and END >= '%d';''' % (uniprot, aa_pos, aa_pos))
        rows = cur.fetchall()
        if rows != ():
            if (c['ProDomains'] != ['None'] and rows[0][1] in c['ProDomains']) or c['ProDomains'] == ['None']:
                pfam = [rows[0][0],rows[0][1]]
    if pfam == ['.','.'] and c['domains_filter']  == 'Y':
        comphets_id2rm.append(comphets_id)
        continue
    
    # AF, coverage
    # For AF, change to percentage, "-1" to "."
    def formatAF(list_of_af):
        if "-1" in list_of_af:
            return ["."]*len(list_of_af)
        else:
            return ["%.3f" % af for af in [float(i)*100 for i in list_of_af]]
    # AF from comphets output
    exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH = formatAF(content[13:21])
    esp_All_EA_AA = formatAF(content[21:24])
    kg_All_AMR_AFR_EAS_EUR_SAS= formatAF(content[24:30])
    # AF,sc ann from mysql
    query = tuple(cols6[0:2] + cols6[3:5])
    if c['build'] == 'GRCh37.75':
        cur.execute('select AF,AC from TAGC128_hg19 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);', query)
        rows = cur.fetchall()
        if rows == ():
            tagc = ['.','.']
        else:
            tagc = ("%.3f,%s" % (rows[0])).split(',')
        cur.execute('select * from ExAC_r03_coverage where (CHROM,POS)= (%s, %s);', query[0:2])
        rows = cur.fetchall()
        if rows == ():
            exac_coverage = '.'
        else:
            exac_coverage = "%.3f" % (rows[0][2])
    if c['build'] == 'GRCh38.79':
        cur.execute('select AF,AC from TAGC128_hg38 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);', query)
        rows = cur.fetchall()
        if rows == ():
            tagc = ['.','.']
        else:
            tagc = ("%.3f,%s" % (rows[0])).split(',')
        
    # Write out candidates
        
    pop_code_mapping = {'esp_AA':float(esp_All_EA_AA[2] if esp_All_EA_AA[2] != '.' else -1),
                        'esp_ALL':float(esp_All_EA_AA[0] if esp_All_EA_AA[0] != '.' else -1),
                        'esp_EA':float(esp_All_EA_AA[1] if esp_All_EA_AA[1] != '.' else -1),
                        'exac_AFR':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[1] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[1] != '.' else -1),
                        'exac_ALL':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[0] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[0] != '.' else -1),
                        'exac_AMR':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[2] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[2] != '.' else -1),
                        'exac_EAS':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[3] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[3] != '.' else -1),
                        'exac_FIN':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[4] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[4] != '.' else -1),
                        'exac_NFE':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[5] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[5] != '.' else -1),
                        'exac_OTH':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[7] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[7] != '.' else -1),
                        'exac_SAS':float(exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[6] if exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH[6] != '.' else -1),
                        'kg_ALL':float(kg_All_AMR_AFR_EAS_EUR_SAS[0] if kg_All_AMR_AFR_EAS_EUR_SAS[0] != '.' else -1),
                        'kg_AFR':float(kg_All_AMR_AFR_EAS_EUR_SAS[2] if kg_All_AMR_AFR_EAS_EUR_SAS[2] != '.' else -1),
                        'kg_AMR':float(kg_All_AMR_AFR_EAS_EUR_SAS[1] if kg_All_AMR_AFR_EAS_EUR_SAS[1] != '.' else -1),
                        'kg_EAS':float(kg_All_AMR_AFR_EAS_EUR_SAS[3] if kg_All_AMR_AFR_EAS_EUR_SAS[3] != '.' else -1),
                        'kg_EUR':float(kg_All_AMR_AFR_EAS_EUR_SAS[4] if kg_All_AMR_AFR_EAS_EUR_SAS[4] != '.' else -1),
                        'kg_SAS':float(kg_All_AMR_AFR_EAS_EUR_SAS[5] if kg_All_AMR_AFR_EAS_EUR_SAS[5] != '.' else -1),
                        'tagc_AJ':float(tagc[0]) if tagc[0] != '.' else -1,
                        'exac':max([float(i) for i in [-1 if j == '.' else j for j in exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH]]),
                        'esp':max([float(i) for i in [-1 if j == '.' else j for j in esp_All_EA_AA]]),
                        'kg':max([float(i) for i in [-1 if j == '.' else j for j in kg_All_AMR_AFR_EAS_EUR_SAS]]),
                        'tagc':float(tagc[0]) if tagc[0] != '.' else -1}
    
    # AF filter
    samples = variant_pair[4].split(",")
    samples_pop_code = [sample_pop[i] if sample_pop[i] in pop_code else 'exac_ALL' for i in samples]
    if max([pop_code_mapping[code] for code in samples_pop_code]) > float(c['af']):
        comphets_id2rm.append(comphets_id)
        continue
    # Output formatting
    out_content = cols6 + variant_pair + impact + [gname,entrez,ensg,enst,biotype]+hgvs+[uniprot] + pfam
    exac_col = ([exac_coverage] + exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH) if c['build'] == 'GRCh37.75' else (exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH)
    esp_col = esp_All_EA_AA
    kg_col = kg_All_AMR_AFR_EAS_EUR_SAS
    if 'ExAC' in c['ann_af']:
        out_content += exac_col
    if 'ESP6500' in c['ann_af']:
        out_content += esp_col
    if '1000Genome' in c['ann_af']:
        out_content += kg_col
    if 'TAGC' in c['ann_af']:
        out_content += tagc
    fh_ANN.write("\t".join(out_content) + '\n')

fh_ANN.close()


fh_ANN2rc = open('../query/output/{0}/{0}.alts.filter.eff.spc.cs.ANN.txt'.format(c['qid'])).readlines()
fh_rcall = open('../query/output/{0}/Cross-sample_all_variants.txt'.format(c['qid']),'w')

# Store #fm and ss for each compound pair
#devide by 2 for recurrent file
compid_fm = {} 
compid_ss = {}
compid_fmss = {}
fh_rcall.write(fh_ANN2rc[0])
for line in fh_ANN2rc[1:]:
    identifiler = tuple([line.strip().split("\t")[7],line.strip().split("\t")[11]])
    if identifiler in comphets_id2rm:
        continue
    fh_rcall.write(line)    
    compid = identifiler[1]
    n_ss = len(line.strip().split("\t")[10].split(","))
    if compid not in compid_ss:
        compid_ss[compid] = n_ss
    else:
        compid_ss[compid] += n_ss
    if compid not in compid_fm:
        compid_fm[compid] = 1 if n_ss > 1 else 0
    else:
        compid_fm[compid] += 1 if n_ss > 1 else 0
    if compid not in compid_fmss:
        compid_fmss[compid] = 1
    else:
        compid_fmss[compid] += 1
        
fh_rcall.close()


fh_rcall = open('../query/output/{0}/Cross-sample_all_variants.txt'.format(c['qid'])).readlines()
fh_rc = open("../query/output/{0}/Cross-sample_recurrent_variants.txt".format(c['qid']), 'w')
fh_rc.write(fh_rcall[0])
for line in fh_rcall[1:]:
    compid = line.strip().split("\t")[11]
    ss = compid_ss[compid]/2
    fm = compid_fm[compid]/2
    fmss = compid_fmss[compid]/2
    
    if fmss >= c['rc_F_lower'] and fmss <= c['rc_F_upper'] and ss >= c['rc_I_lower'] and ss <= c['rc_I_upper'] and fm >= c['family']:  
        fh_rc.write(line)
        
fh_rc.close()



