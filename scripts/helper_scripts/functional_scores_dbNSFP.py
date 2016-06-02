# Title: functional_scores_dbNSFP.py
# Description: This script is to annotate and/or filter variants on functional prediciton scores from dbNSFP.
# Author: Siwei Chen
# Date created: 27/08/2015
# Date last modified: 20/03/2016
# Python version: 2.7.5

import argparse
import os
import imp


parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
parser.add_argument("-f", help="input file")
#!! input file path can NOT be absolute!!
args = parser.parse_args()

c = imp.load_source('cache', args.c).cache

# Prediction algorisms
function = ['SIFT', 'PROVEAN', 'PPH2_HDIV','PPH2_HVAR', 'LRT', 'MT', 'MA', 'FATHMM','FATHMMMKL','VEST3','CADD_RAW','CADD_PHRED','DANN', 'MetaSVM','MetaLR', 'fitCons']
conservation = ['GERP', 'phyloP_vertebrate','phyloP_mammalian','phastCons_vertebrate','phastCons_mammalian', 'SiPhy']

# Get user's selections
function = [i for i in function if i in c['ann_fip']]
conservation = [i for i in conservation if i in c['ann_fip']]

# Map prediction algorisms to column number for querying dbNSFP database
dbNSFP_col = {'SIFT':'24,26,', 'PROVEAN': '53,55,', 'PPH2_HDIV': '30,32,','PPH2_HVAR':'33,35,', 'LRT':'36,38,', 'MT':'40,42,', 'MA':'47,49,', 'FATHMM':'50,52,','FATHMMMKL':'65,67,','VEST3':'58,',
              'CADD_RAW':'60,', 'CADD_PHRED':'62,','DANN':'63,', 'MetaSVM':'69,71,','MetaLR':'72,74,', 'fitCons':'76,','GERP':'89,', 'phyloP_vertebrate':'91,','phyloP_mammalian':'93,','phastCons_vertebrate':'95,','phastCons_mammalian':'97,', 'SiPhy':'100,'}

query_fip = False
query_cons = False

# Variant functional predictions
if function != []:
    os.system("""awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$4,$5}}' {0} | sed '/^CHROM/d' | sed '/synonymous_variant/d' \
              > {1}/dbNSFP_function.in""".format(args.f,os.path.dirname(args.f)))
    os.chdir(os.path.dirname(args.f))
    indir=os.getcwd()
    #open("{1}/dbNSFP_function.out".format(args.f,os.path.dirname(args.f)),'a').close()

    searchdb_cols = ''
    for i in function:
        searchdb_cols += dbNSFP_col[i]
    # Use dbNSFP search program to query prediction scores for variants
    # (have to be under the directory of the program)
    os.chdir('../../GeMSTONE-data/dbNSFPv3.1a/') 
    if c['build'] == 'GRCh38.79':
        os.system("""java -Xmx5g search_dbNSFP31a \
                  -i {0}/dbNSFP_function.in \
                  -o {0}/dbNSFP_function.out \
                  -w {1}""".format(indir,'1,2,3,4,'+searchdb_cols[:-1]))      
    if c['build'] == 'GRCh37.75':
        os.system("""java -Xmx5g search_dbNSFP31a -v hg19 \
                  -i {0}/dbNSFP_function.in \
                  -o {0}/dbNSFP_function.out \
                  -w {1}""".format(indir,'8,9,3,4,'+searchdb_cols[:-1]))         
    os.chdir('/data/germline-web/scripts/')

    out_header= open('{0}/dbNSFP_function.out'.format(indir)).readlines()[0][:-1].split('\t')
    # Index entries with more than one returned values (score and binned prediction) for formatting purpose
    sift_idx = -1
    provean_idx= -1
    pph2_hdiv_idx= -1
    pph2_hvar_idx= -1
    mt_idx= -1
    vest_idx= -1
    fathmm_idx= -1
    if 'SIFT_score' in out_header:
        sift_idx = out_header.index('SIFT_score')
    if 'PROVEAN_score' in out_header:
        provean_idx = out_header.index('PROVEAN_score')
    if 'Polyphen2_HDIV_score' in out_header:
        pph2_hdiv_idx = out_header.index('Polyphen2_HDIV_score')
    if 'Polyphen2_HVAR_score' in out_header:
        pph2_hvar_idx = out_header.index('Polyphen2_HVAR_score')    
    if 'MutationTaster_score' in out_header:
        mt_idx = out_header.index('MutationTaster_score')
    if 'FATHMM_score' in out_header:
        fathmm_idx= out_header.index('FATHMM_score')    
    if 'VEST3_score' in out_header:
        vest_idx= out_header.index('VEST3_score')
    query_fip= {}
    # Parse and format output from dbNSFP
    for preds in open('{0}/dbNSFP_function.out'.format(indir)).readlines()[1:]:
        k= tuple(preds.split('\t')[:4]) # CHROM,POS,REF,ALT as identifier/key for each variant
        v= preds[:-1].split('\t')       # predictions as value 
        if sift_idx != -1:
            sift_score= v[sift_idx].split(';')
            sift_pred= v[sift_idx+1].split(';')
            if set(sift_score) == set(['.']):
               sift_score= sift_pred= ['.']
            elif len(sift_score) != 1:
                sift_score= [str(min([float(i) for i in sift_score if i != '.']))]
                sift_pred= ['D' if 'D' in sift_pred else 'T']
            v[sift_idx]= sift_score[0]
            v[sift_idx+1]= sift_pred[0]            
        if provean_idx != -1:    
            provean_score= v[provean_idx].split(';')
            provean_pred= v[provean_idx+1].split(';')
            if set(provean_score) == set(['.']):
               provean_score= provean_pred= ['.']
            elif len(provean_score) != 1:
                provean_score= [str(min([float(i) for i in provean_score if i != '.']))]
                provean_pred= ['D' if 'D' in provean_pred else 'N']
            v[provean_idx]= provean_score[0]
            v[provean_idx+1]= provean_pred[0]        
        if pph2_hdiv_idx != -1:
            pph2_hdiv_score= v[pph2_hdiv_idx].split(';')
            pph2_hdiv_pred= v[pph2_hdiv_idx+1].split(';')
            if set(pph2_hdiv_score) == set(['.']):
                pph2_hdiv_score= pph2_hdiv_pred= ['.']
            elif len(pph2_hdiv_score) != 1:
                pph2_hdiv_score= [str(max([float(i) for i in pph2_hdiv_score if i != '.']))]
                pph2_hdiv_pred= ['D' if 'D' in pph2_hdiv_pred else 'P' if 'P' in pph2_hdiv_pred else 'B']
            v[pph2_hdiv_idx]= pph2_hdiv_score[0]
            v[pph2_hdiv_idx+1]= pph2_hdiv_pred[0]
        if pph2_hvar_idx != -1:
            pph2_hvar_score= v[pph2_hvar_idx].split(';')
            pph2_hvar_pred= v[pph2_hvar_idx+1].split(';')
            if set(pph2_hvar_score) == set(['.']):
                pph2_hvar_score= pph2_hvar_pred= ['.']
            elif len(pph2_hvar_score) != 1:
                pph2_hvar_score= [str(max([float(i) for i in pph2_hvar_score if i != '.']))]
                pph2_hvar_pred= ['D' if 'D' in pph2_hvar_pred else 'P' if 'P' in pph2_hvar_pred else 'B']
            v[pph2_hvar_idx]= pph2_hvar_score[0]
            v[pph2_hvar_idx+1]= pph2_hvar_pred[0]
        if mt_idx != -1:
            mt_score= v[mt_idx].split(';')
            mt_pred= v[mt_idx+1].split(';')
            if set(mt_score) == set(['.']):
                mt_score = mt_pred = ['.']
            elif len(mt_score) != 1:
                mt_score= [str(max([float(i) for i in mt_score if i != '.']))]
                mt_pred= ['A' if 'A' in mt_pred else 'D' if 'D' in mt_pred else 'N' if 'N' in mt_pred else 'P']
            v[mt_idx]= mt_score[0]
            v[mt_idx+1]= mt_pred[0]
        if fathmm_idx != -1:
            fathmm_score= v[fathmm_idx].split(';')
            fathmm_pred= v[fathmm_idx+1].split(';')
            if set(fathmm_score) == set(['.']):
                fathmm_score= fathmm_pred = ['.']
            elif len(fathmm_score) != 1:
                fathmm_score= [str(min([float(i) for i in fathmm_score if i != '.']))]
                fathmm_pred= ['D' if 'D' in fathmm_pred else 'T']
            v[fathmm_idx]= fathmm_score[0]
            v[fathmm_idx+1]= fathmm_pred[0]
        if vest_idx != -1:
            vest_score= v[vest_idx].split(';')
            if set(vest_score) == set(['.']):
                vest_score= ['.']
            elif len(vest_score) != 1:
                vest_score= [str(max([float(i) for i in vest_score if i != '.']))]
            v[vest_idx]= vest_score[0]
        # Store prediction scores of variants    
        query_fip[k]= v[4:]
        
# Variant conservation scores
if conservation != []:
    os.system("""awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2}}' {0} | sed '/^CHROM/d' \
              > {1}/dbNSFP_conservation.in""".format(args.f,os.path.dirname(args.f)))
    os.chdir(os.path.dirname(args.f))
    indir=os.getcwd()
    #open("{1}/dbNSFP_function.out".format(args.f,os.path.dirname(args.f)),'a').close()   
    searchdb_cols = ''
    for i in conservation:
        searchdb_cols += dbNSFP_col[i]          
    os.chdir('../../GeMSTONE-data/dbNSFPv3.1a/') 
    if c['build'] == 'GRCh38.79':
        os.system("""java -Xmx5g search_dbNSFP31a \
                  -i {0}/dbNSFP_conservation.in \
                  -o {0}/dbNSFP_conservation.out \
                  -w {1}""".format(indir,'1,2,'+searchdb_cols[:-1]))       
    if c['build'] == 'GRCh37.75':
        os.system("""java -Xmx5g search_dbNSFP31a -v hg19 \
                  -i {0}/dbNSFP_conservation.in \
                  -o {0}/dbNSFP_conservation.out \
                  -w {1}""".format(indir,'8,9,'+searchdb_cols[:-1]))        
    os.chdir('/data/germline-web/scripts/')
    
    # Store conservation scores of variants
    # CHROM,POS as identifier/key for each variant,
    # conservation score as value
    query_cons = dict([tuple(line.strip().split("\t")[:2]),line.strip().split("\t")[2:]] for line in open('{0}/dbNSFP_conservation.out'.format(indir)).readlines()[1:])

# Read in variant file for annotation and/or filtering
# Write out annotated variant file
fh_2fip = open(args.f).readlines()
fh_fip = open(args.f[:-3]+'fip.txt', 'w')
fip_header = {'SIFT': '\tSIFT_SCORE\tSIFT_PREDICTION',
              'PROVEAN': '\tPROVEAN_SCORE\tPROVEAN_PREDICTION',
              'PPH2_HDIV': '\tPPH2_HDIV_SCORE\tPPH2_HDIV_PREDICTION',
              'PPH2_HVAR': '\tPPH2_HVAR_SCORE\tPPH2_HVAR_PREDICTION',
              'LRT': '\tLRT_SCORE\tLRT_PREDICTION',
              'MT': '\tMT_SCORE\tMT_PREDICTION',
              'MA': '\tMA_SCORE\tMA_PREDICTION',
              'FATHMM': '\tFATHMM_SCORE\tFATHMM_PREDICTION',
              'FATHMMMKL': '\tFATHMM-MKL_SCORE\tFATHMM-MKL_PREDICTION',
              'VEST3': '\tVEST3',
              'CADD_RAW':'\tCADD_RAW_SCORE',
              'CADD_PHRED':'\tCADD_PHRED_SCORE',
              'DANN':'\tDANN',
              'MetaSVM': '\tMetaSVM_SCORE\tMetaSVM_PREDICTION',
              'MetaLR':'\tMetaLR_SCORE\tMetaLR_PREDICTION',
              'fitCons': '\tfitCons',
              'GERP': '\tGERP++',
              'phyloP_vertebrate': '\tphyloP7way_vertebrate',
              'phyloP_mammalian':'\tphyloP20way_mammalian',
              'phastCons_vertebrate':'\tphastCons7way_vertebrate',
              'phastCons_mammalian':'\tphastCons20way_mammalian',
              'SiPhy': '\tSiPhy'}

# Format output header and contents
score_header = ''.join([fip_header[i] for i in function+conservation]).split("\t")[1:]
fh_fip.write(fh_2fip[0][:-1] + '\t' + '\t'.join(score_header) + '\t(Deleterious_prediction_counts)\n')
filter_handle = {'SIFT': 'SIFT_SCORE',
              'PROVEAN': 'PROVEAN_SCORE',
              'PPH2_HDIV': 'PPH2_HDIV_SCORE',
              'PPH2_HVAR': 'PPH2_HVAR_SCORE',
              'LRT': 'LRT_SCORE',
              'MT': 'MT_SCORE',
              'MA': 'MA_SCORE',
              'FATHMM': 'FATHMM_SCORE',
              'FATHMMMKL': 'FATHMM-MKL_SCORE',
              'VEST3': 'VEST3',
              'CADD_RAW':'CADD_RAW_SCORE',
              'CADD_PHRED':'CADD_PHRED_SCORE',
              'DANN':'DANN',
              'MetaSVM': 'MetaSVM_SCORE',
              'MetaLR':'MetaLR_SCORE',
              'fitCons': 'fitCons',
              'GERP': 'GERP++',
              'phyloP_vertebrate': 'phyloP7way_vertebrate',
              'phyloP_mammalian':'phyloP20way_mammalian',
              'phastCons_vertebrate':'phastCons7way_vertebrate',
              'phastCons_mammalian':'phastCons20way_mammalian',
              'SiPhy': 'SiPhy'}
filter_idx = {} # Index of score in the score_header for filtering(counting del)
for tool in function+conservation:
    filter_idx[tool] = score_header.index(filter_handle[tool])

if c['inheritance'] != 'Rec_compound':
    for var in fh_2fip[1:]:
        content= var[:-1].split('\t')
        fquery= tuple(content[0:2] + content[3:5])
        cquery = tuple(content[:2])
        fip= (query_fip[fquery] if fquery in query_fip else ['.']*len(query_fip.values()[0])) if function else []
        cons = (query_cons[cquery] if cquery in query_cons else ['.']*len(query_cons.values()[0])) if conservation else []
        score_content = fip+cons
        if 'missense' not in content[7]:
            #del_count = content[7]
            fh_fip.write('\t'.join(content+score_content) + '\t.\n')
        else:
            del_count= 0
            for tool in function+conservation:
                if score_content[filter_idx[tool]] != '.':
                    score = float(score_content[filter_idx[tool]])
                    if score >= float(c['fip_cutoffs'][tool][0]) and score <= float(c['fip_cutoffs'][tool][1]):
                        del_count += 1
            if del_count >= c["del_count"]:
                fh_fip.write('\t'.join(content+score_content) + '\t' + str(del_count) +'\n')    
    fh_fip.close()

else:
    out_content = []
    comphets_id2rm = []
    for var in fh_2fip[1:]:
        content= var[:-1].split('\t')
        fquery= tuple(content[0:2] + content[3:5])
        cquery = tuple(content[:2])
        fip= (query_fip[fquery] if fquery in query_fip else ['.']*len(query_fip.values()[0])) if query_fip else []
        cons = (query_cons[cquery] if cquery in query_cons else ['.']*len(query_cons.values()[0])) if query_cons else []
        score_content = fip+cons
        if 'missense' not in content[13]:
            #fh_fip.write('\t'.join(content+score_content) + '\t.\n')
            out_content.append('\t'.join(content+score_content) + '\t.\n')
        else:
            del_count= 0
            for tool in function+conservation:
                if score_content[filter_idx[tool]] != '.':
                    score = float(score_content[filter_idx[tool]])
                    if score >= float(c['fip_cutoffs'][tool][0]) and score <= float(c['fip_cutoffs'][tool][1]):
                        del_count += 1
            if del_count >= c["del_count"]:
                #fh_fip.write('\t'.join(content+score_content) + '\t' + str(del_count) +'\n')
                out_content.append('\t'.join(content+score_content) + '\t' + str(del_count) +'\n')
            else:
                comphets_id2rm.append(tuple([content[7],content[11]]))
    for line in out_content:
        if tuple([line.strip().split("\t")[7],line.strip().split("\t")[11]]) not in comphets_id2rm:
            fh_fip.write(line)
    fh_fip.close()

