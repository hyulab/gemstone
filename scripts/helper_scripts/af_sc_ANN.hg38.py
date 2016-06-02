# Title: af_sc_ANN.hg38.py
# Description: This script is to annotate and/or filter variants for allele frequency, variant consequence and location on genome build GRCh38/hg38.
# Author: Siwei Chen
# Date created: 18/03/2016
# Date last modified: 25/04/2016
# Python version: 2.7.5

import argparse
import os
import imp
import MySQLdb
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
parser.add_argument("-p", help="population")
args = parser.parse_args()
c = imp.load_source('cache', args.c).cache



mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

fh_FM = open('../query/output/{0}/{0}.alts.filter.spc.eff.cs.txt'.format(c['qid'], c['vcf_filename'])).readlines()
fh_ANN = open('../query/output/{0}/{0}.alts.filter.spc.eff.cs.af2ann.txt'.format(c['qid'], c['vcf_filename']), 'w')
fh_ANN.write('CHROM\tPOS\tID\tREF\tALT\tFILTER\tINDIVIDUAL_ID\tCONSEQUENCE\tPUTATIVE_IMPACT\tGENE_NAME\tENTREZ_ID\tENSEMBL_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_BIOTYPE\tHGVS.c\tHGVS.p\t' +
             'UNIPROT\tPFAM_ID\tDOMAIN\t' + 
             'ExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_SAS\tExAC_OTH\t' +  
             'ESP6500_ALL\tESP6500_EA\tESP6500_AA\t' + 
             '1000G_ALL\t1000G_AMR\t1000G_AFR\t1000G_EAS\t1000G_EUR\t1000G_SAS\t' +  
             'TAGC_AF\tTAGC_AC\n') #ExAC_MEAN_COVERAGE not available for hg38
for var in fh_FM:
    # Parsing
    if var.split('\t')[7].split('|')[7] in c['transcript_biotype']:
        chr_pos_id_ref_alt = var.split('\t')[0:5]
        fil = [var.split('\t')[6]]
        indv_id = [','.join(open('../query/output/{0}/{0}.fm_aff.txt'.format(c['qid'], c['ped_filename'])).read()[:-1].split('\n'))]
        so = var.split('\t')[7].split('|')[1]
        if not [i for i in c['consequences'] if i in so]: # string matching
            continue
        cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp = var.split('\t')[7].split('|')[1:4] + ['.'] + var.split('\t')[7].split('|')[4:8] + var.split('\t')[7].split('|')[9:11] + ['.','.','.'] 
        # Map gname to gid (new version of snpEff does not give gid yet but giving dup cols of gname!)
        gname = cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[2]
        cur.execute('select ID from gID_gSymbol where Symbol= %s;', gname)
        rows = cur.fetchall()
        if rows == ():
            cur.execute('select ID from gSymbol_gID_biomart where Symbol= %s;', gname)
            rows = cur.fetchall()
            if rows != ():
                cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[3] = rows[0][0] 
        else:
            cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[3] = rows[0][0]
        if cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[3] != '.':
            cur.execute('select Uniprot from gID_Uniprot where Gene_ID = %s;', cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[3])
            rows = cur.fetchall()
            if rows != ():
                cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[10]= rows[0][0]
                hgvsp= cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[9]
                if len(hgvsp) > 5 and hgvsp[0:2] == 'p.':
                    aa_pos= int("".join(itertools.takewhile(str.isdigit, hgvsp[5:])))
                    uniprot= cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[10]
                    if c['domains'] == 'All' or c['ProDomains'] != ['None']:
                        cur.execute('''select PFAM,DOMAIN_NAME from UNIPROT_PFAM_DOMAINS where UNIPROT= '%s' and START <= '%d' and END >= '%d';''' % (uniprot, aa_pos, aa_pos))
                        rows = cur.fetchall()
                        if rows != ():
                            if (c['ProDomains'] != ['None'] and rows[0][1] in c['ProDomains']) or c['ProDomains'] == ['None']:
                                cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[11]= rows[0][0]
                                cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[12]= rows[0][1]

        # Annotate af, coverage
        query = tuple(var.split('\t')[0:2] + var.split('\t')[3:5])
        cur.execute('select AF_Adj,AF_AFR,AF_AMR,AF_EAS,AF_FIN,AF_NFE,AF_SAS, AF_OTH from ExAC_hg38 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);', query)
        rows = cur.fetchall()
        if rows == ():
            exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH = ['.', '.', '.', '.', '.', '.', '.', '.']
        else:
            exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH = ("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % (rows[0])).split(',')
    
        cur.execute('select * from ESP6500_hg19 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);', query)
        rows = cur.fetchall()
        if rows == ():
            esp_All_EA_AA = ['.', '.', '.']
        else:
            esp_All_EA_AA = ["%.3f" % (rows[0][6]), "%.3f" % (rows[0][4]), "%.3f" % (rows[0][5])]
    
        cur.execute('select * from 1000G_hg38 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);', query)
        rows = cur.fetchall()
        if rows == ():
            kg_All_AMR_AFR_EAS_EUR_SAS = ['.','.', '.', '.', '.', '.']
        else:
            kg_All_AMR_AFR_EAS_EUR_SAS = ("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % (rows[0][5:11])).split(',')
    
        cur.execute('select AF,AC from TAGC128_hg38 where (CHROM,POS,REF,ALT)= (%s, %s, %s, %s);', query)
        rows = cur.fetchall()
        if rows == ():
            tagc = ['.','.']
        else:
            tagc = ("%.3f,%s" % (rows[0])).split(',')
    
    
        if c['domains_filter'] == ['Y'] and cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp[11] == '.':
            continue
        fh_ANN.write('\t'.join(chr_pos_id_ref_alt + fil + indv_id + cons_imp_gname_entr_ensg_ftype_fid_trans_hgvscp +
                               exac_All_AFR_AMR_EAS_FIN_NFE_SAS_OTH + esp_All_EA_AA + kg_All_AMR_AFR_EAS_EUR_SAS + tagc) + '\n')

fh_ANN.close()

# Map population code to column index for filtering by awk
ann1_col = {'ExAC': '$21,$22,$23,$24,$25,$26,$27,$28,', 'ESP6500': '$29,$30,$31,', '1000Genome': '$32,$33,$34,$35,$36,$37,', 'TAGC': '$38,$39,'}
pop_col = {'kg_EAS': '$35', 'kg_ALL': '$32', 'exac_AMR': '$23', 'exac_NFE': '$26', 'kg_AFR': '$34', 'kg_AMR': '$33', 'exac_AFR': '$22', 'esp_ALL': '$29', 'kg_EUR': '$36', 'exac_FIN': '$25', 'tagc_AJ_AC': '$39', 'exac_SAS': '$27', 'esp_EA': '$30', 'tagc_AJ': '$38', 'esp_AA': '$31', 'exac_EAS': '$24', 'exac_OTH': '$28', 'exac_ALL': '$21', 'kg_SAS': '$37'}
pop_sub = {'exac':['exac_ALL','exac_AFR','exac_AMR','exac_EAS','exac_FIN','exac_NFE','exac_SAS','exac_OTH'],
            'esp':['esp_ALL','esp_EA','esp_AA'],
            'kg':['kg_ALL','kg_AMR','kg_AFR','kg_EAS','kg_EUR','kg_SAS']}
pop = []
for i in  args.p.split(','):
    if i in pop_sub:
        pop += pop_sub[i]
    else:
        pop.append(i)

cols = '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,' if c['domains'] != 'N' else '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,'

ann_af = [i for i in ['ExAC', 'ESP6500', '1000Genome','TAGC'] if i in c["ann_af"]]

for i in ann_af:
    cols += ann1_col[i]

# Filter and annotate variants using awk
os.system('''awk 'BEGIN {FS= "\t"; OFS= "\t"}; ''' + ''' ($1 == "CHROM") {print ''' + cols[:-1] + ''' }' ''' + '../query/output/{0}/{0}.alts.filter.spc.eff.cs.af2ann.txt'.format(
    c['qid'], c['vcf_filename']) + ' > ' + '../query/output/{0}/{0}.alts.filter.spc.eff.cs.ANN.txt'.format(c['qid'], c['vcf_filename']))
os.system('''awk 'BEGIN {FS= "\t"; OFS= "\t"}; (''' + ' && '.join([pop_col[p] + "<=" + c['af'] for p in pop]) + ''') {print ''' + cols[:-1] + '''}' ''' +
          '../query/output/{0}/{0}.alts.filter.spc.eff.cs.af2ann.txt '''.format(c['qid'], c['vcf_filename']) + '''>> ../query/output/{0}/{0}.alts.filter.spc.eff.cs.ANN.txt'''.format(c['qid'], c['vcf_filename']))














