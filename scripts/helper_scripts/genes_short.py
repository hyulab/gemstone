# Title: genes_short.py
# Description: This script is to annotate and/or filter genes for interaction partners, Gene Ontology terms, known associated diseases (from HGMD, Clinver, OMIM), pathways (KEGG, Reactome, BioCarta), and user's list of interested genes.
# Author: Siwei Chen
# Date created: 05/09/2015
# Date last modified: 25/04/2016
# Python version: 2.7.5

import imp
import MySQLdb
mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()
import argparse
import subprocess as sub


parser = argparse.ArgumentParser()
parser.add_argument("-c", help="cache file")
parser.add_argument("-f", help="variant file/gene list")
parser.add_argument("-o", help="output file")
parser.add_argument("-l", help="risk gene list (optional)")
args = parser.parse_args()

c = imp.load_source('cache', args.c).cache

fh_2dp= open(args.f).readlines()
fh_dp= open(args.o,'w')

# Output header according to user's selections
header = ['INTERACTIONS_IntAct',
    'INTERACTIONS_BioGRID',
    'INTERACTIONS_ConsensusPathDB',
    'INTERACTIONS_HINT',
    'GO_BP\t(interaction_partners)',
    'GO_MF\t(interaction_partners)',
    'GO_CC\t(interaction_partners)',
    'HGMD\t(interaction_partners)',
    'OMIM\t(interaction_partners)',
    'CLINVAR\t(interaction_partners)',
    'MY_CANDIDATES'] if c['go_interactors'] == ['y'] else [
    'INTERACTIONS_IntAct',
    'INTERACTIONS_BioGRID',
    'INTERACTIONS_ConsensusPathDB',
    'INTERACTIONS_HINT',
    'GO_BP',
    'GO_MF',
    'GO_CC',
    'HGMD\t(interaction_partners)',
    'OMIM\t(interaction_partners)',
    'CLINVAR\t(interaction_partners)',
    'MY_CANDIDATES']
header2 = ['KEGG\t(interaction_partners)',
    'BIOCARTA\t(interaction_partners)',
    'REACTOME\t(interaction_partners)'] if c["pathway_interactors"] == ['y'] else ['KEGG','BIOCARTA','REACTOME']
header += header2

# Map header items with parameter dict entries
header_argv = {'GO_BP\t(interaction_partners)': c['go_bp'],
         'GO_MF\t(interaction_partners)': c['go_mf'],
         'GO_CC\t(interaction_partners)': c['go_cc'],
         'HGMD\t(interaction_partners)': c['dg_hgmd'],
         'OMIM\t(interaction_partners)': c['dg_omim'],
         'CLINVAR\t(interaction_partners)': c['dg_clinvar'],
         'MY_CANDIDATES': c['gl_user']} if c['go_interactors'] == ['y'] else {
         'GO_BP': c['go_bp'],
         'GO_MF': c['go_mf'],
         'GO_CC': c['go_cc'],
         'HGMD\t(interaction_partners)': c['dg_hgmd'],
         'OMIM\t(interaction_partners)': c['dg_omim'],
         'CLINVAR\t(interaction_partners)': c['dg_clinvar'],
         'MY_CANDIDATES': c['gl_user']}

header_argv2 = {'KEGG\t(interaction_partners)': c['pw_kegg'],
         'BIOCARTA\t(interaction_partners)': c['pw_biocarta'],
         'REACTOME\t(interaction_partners)': c['pw_reactome']} if c["pathway_interactors"] == ['y'] else {
         'KEGG': c['pw_kegg'],
         'BIOCARTA': c['pw_biocarta'],
         'REACTOME': c['pw_reactome']}
header_argv.update(header_argv2)

# Header to output file 
dp_hd= fh_2dp[0][:-1]
fh_dp.write(dp_hd + '\t' + '\t'.join([h for h in header if (h in header_argv and header_argv[h] != ['None'] and header_argv[h] != '' and header_argv[h] != 'None') or h.split('_')[-1] in c['proteinprotein']]) + '\n')
# Annotations and/or filterings 
for var in fh_2dp[1:]:
    line = var[:-1]
    # ppi
    dp_interaction_partners = []
    if c['proteinprotein'] != ['None']:
        intDB= [i for i in ['IntAct','BioGRID','ConsensusPathDB'] if i in c['proteinprotein']]
        if intDB != []:
            cur.execute('select {0} from Gene_InteractionDB where Ensembl_gene = %s;'.format(','.join(intDB)),var[:-1].split('\t')[2])
            rows= cur.fetchall() # by gene names
            if rows == ():
                interaction_partners = ["."]*len(intDB)
            else:
                interaction_partners= [i.strip(";") for i in rows[0]]
            for i in interaction_partners:
                if len(i.split(";")) <=20:
                    line += '\t' + i
                else:
                    line += '\t' + str(len(i.split(";")))
            dp_interaction_partners += ";".join(interaction_partners).split(";")
        if 'HINT' in c['proteinprotein']:
            cur.execute('select Gene_B from HINT_interactors where Gene_A = %s;', var[:-1].split('\t')[0])
            row1 = cur.fetchall()
            cur.execute('select Gene_A from HINT_interactors where Gene_B = %s;', var[:-1].split('\t')[0])
            row2 = cur.fetchall()
            if row1 + row2 != ():
                hint = ';'.join(list(set([i for i in [i[0] for i in row1 + row2]])))
            else:
                hint = '.'
            line += '\t'+ hint
            dp_interaction_partners += hint.split(";")
        dp_interaction_partners = list(set([i for i in dp_interaction_partners if i != '.'])) 
    
    # Gene Ontology (GO)
    go_count = 0
    # GO biological process
    go_bp = []
    if c['go_bp'] != ['None']:
        go_bp_genes=[]
        bp_ann = '.'
        for bp in c['go_bp']:
            cur.execute('select EntrezID from GO_biological_process where BP = %s;', bp)
            rows = cur.fetchall()
            go_bp_genes += rows[0][0].split(',')
            if var[:-1].split('\t')[1] in rows[0][0].split(','):
              go_bp.append(bp.replace('_', ' ').title())
        if go_bp != []:
            bp_ann = ', '.join(go_bp)
            go_count += 1
        if c['go_interactors'] == ['y']:
            bp_int_ann = '.'
            if go_bp_genes != []:
                cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(go_bp_genes))) 
                rows = cur.fetchall()
                go_bp_genes_names = list(set([i[0] for i in rows]))
                go_bp_int = [i for i in go_bp_genes_names if i in dp_interaction_partners]
                if go_bp_int != []:
                    bp_int_ann = ';'.join(go_bp_int)
            line += '\t' + '\t'.join([bp_ann,bp_int_ann])
        else:
            line += '\t' + bp_ann
    # GO molecullar function
    go_mf = []
    if c['go_mf'] != ['None']:
        go_mf_genes=[]
        mf_ann = '.'
        
        for mf in c['go_mf']:
            cur.execute('select EntrezID from GO_molecular_function where MF = %s;', mf)
            rows = cur.fetchall()
            go_mf_genes += rows[0][0].split(',')
            if var[:-1].split('\t')[1] in rows[0][0].split(','):
              go_mf.append(mf.replace('_', ' ').title())
        if go_mf != []:
            mf_ann = ', '.join(go_mf)
            go_count += 1
        if c['go_interactors'] == ['y']:
            mf_int_ann = '.'
            if go_mf_genes != []:
                cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(go_mf_genes))) 
                rows = cur.fetchall()
                go_mf_genes_names = list(set([i[0] for i in rows]))
                go_mf_int = [i for i in go_mf_genes_names if i in dp_interaction_partners]
                if go_mf_int != []:
                    mf_int_ann = ';'.join(go_mf_int)
            line += '\t' + '\t'.join([mf_ann,mf_int_ann])
        else:
            line += '\t' + mf_ann
    # GO cellular component
    go_cc = []
    if c['go_cc'] != ['None']:
        go_cc_genes=[]
        cc_ann = '.'
        
        for cc in c['go_cc']:
            cur.execute('select EntrezID from GO_cellular_component where CC = %s;', cc)
            rows = cur.fetchall()
            go_cc_genes += rows[0][0].split(',')
            if var[:-1].split('\t')[1] in rows[0][0].split(','):
              go_cc.append(cc.replace('_', ' ').title())
        if go_cc != []:
            cc_ann = ', '.join(go_cc)
            go_count += 1
        if c['go_interactors'] == ['y']:
            cc_int_ann = '.'
            if go_cc_genes != []:
                cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(go_cc_genes)))
                rows = cur.fetchall()
                go_cc_genes_names = list(set([i[0] for i in rows]))
                go_cc_int = [i for i in go_cc_genes_names if i in dp_interaction_partners]
                if go_cc_int != []:
                    cc_int_ann = ';'.join(go_cc_int)
            line += '\t' + '\t'.join([cc_ann,cc_int_ann])
        else:
            line += '\t' + cc_ann
            
    if c['go_filter'] == ['y'] and go_count == 0:
        continue
    
    # HGMD
    hgmd = []
    if c['dg_hgmd'] != ['None']:
        hgmd_genes=[]
        hgmd_ann = '.'
        hgmd_int_ann = '.'
        for disease in c['dg_hgmd']:
            cur.execute('select Gene_ID from HGMD_disease_gene where Disease = %s;', disease)
            rows = cur.fetchall()
            hgmd_genes += [i[0] for i in rows]
            if var[:-1].split('\t')[1] in [i[0] for i in rows]:
                hgmd.append(disease)
        if hgmd != []:
            hgmd_ann = ','.join(hgmd)
        if hgmd_genes != []:
            cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join([i for i in hgmd_genes if i != '.'])))
            rows = cur.fetchall()
            hgmd_genes_names = list(set([i[0] for i in rows]))
            hgmd_int = [i for i in hgmd_genes_names if i in dp_interaction_partners]
            if hgmd_int != []:
                hgmd_int_ann = ';'.join(hgmd_int)
        line += '\t' + '\t'.join([hgmd_ann,hgmd_int_ann])
        
    #OMIM
    omim = []
    if c['dg_omim'] != ['None']:
        omim_genes=[] 
        omim_ann = '.'
        omim_int_ann = '.'
        for disease in c['dg_omim']:
            cur.execute('select Gene_Name from OMIM_std_all where Disease = %s;', disease)
            rows = cur.fetchall()
            omim_genes += [i[0] for i in rows]
            if var[:-1].split('\t')[0] in [i[0] for i in rows]:
                omim.append(disease)
        if omim != []:
            omim_ann = ','.join(omim)
        if omim_genes != []:
            omim_int = [i for i in list(set(omim_genes)) if i in dp_interaction_partners]
            if omim_int != []:
                omim_int_ann = ';'.join(omim_int)
        line += '\t' + '\t'.join([omim_ann,omim_int_ann])
        
    #Clinvar
    clinvar = []    
    if c['dg_clinvar'] != ['None']:
        clinvar_genes=[]
        clinvar_ann = '.'
        clinvar_int_ann = '.'
        for disease in c['dg_clinvar']:
            cur.execute('select UniProt from ClinVar where Disease = %s;', disease)
            rows = cur.fetchall()
            uniprots = list(set([i[0] for i in rows]))
            cur.execute("select Gene_ID from gID_Uniprot where Uniprot in ( '{0}' );".format("','".join(uniprots)))
            rows = cur.fetchall()
            clinvar_genes += [i[0] for i in rows]

            if var[:-1].split('\t')[4] in uniprots:
                clinvar.append(disease)
        if clinvar != []:
            clinvar_ann = ','.join(clinvar)
        if clinvar_genes != []:
            cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(clinvar_genes)))
            rows = cur.fetchall()
            clinvar_genes_names = list(set([i[0] for i in rows]))
            clinvar_int = [i for i in clinvar_genes_names if i in dp_interaction_partners]
            if clinvar_int != []:
                clinvar_int_ann = ';'.join(clinvar_int)
        line += '\t' + '\t'.join([clinvar_ann,clinvar_int_ann])

    # User's uploaded gene list of interest 
    mc = '.'
    if c['gl_user'] != 'None':
        risk_genes_ann = {}
        for l in open(args.l).readlines():
            risk_genes_ann[l.strip().split("\t")[0]] = (':'+l.strip().split("\t")[1]) if len(l.strip().split("\t")) > 1 else ''
        for identifier in line.split("\t")[0:3]:
            if identifier in risk_genes_ann:
                mc = identifier + risk_genes_ann[identifier] 
        line += '\t' + mc
    
    # KEGG pathway
    pathway_count = 0    
    kegg = [] 
    if c['pw_kegg'] != ['None']:
        kegg_genes=[]
        kegg_ann = '.'        
        for pw in c['pw_kegg']:
            cur.execute('select EntrezID from KEGG_EntrezID where Pathway = %s;', pw)
            rows = cur.fetchall()
            kegg_genes += rows[0][0].split(',')
            if var[:-1].split('\t')[1] in rows[0][0].split(','):
              kegg.append(pw.replace('_', ' ').title())
        if kegg != []:
            kegg_ann = ', '.join(kegg)
            pathway_count += 1
        if c["pathway_interactors"] == ['y']:
            kegg_int_ann = '.'
            if kegg_genes != []:
                cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(kegg_genes)))
                rows = cur.fetchall()
                kegg_genes_names = list(set([i[0] for i in rows]))
                kegg_int = [i for i in kegg_genes_names if i in dp_interaction_partners]
                if kegg_int != []:
                    kegg_int_ann = ';'.join(kegg_int)
            line += '\t' + '\t'.join([kegg_ann,kegg_int_ann])
        else:
            line += '\t' + kegg_ann

    # BioCarta pathway
    biocarta = []
    if c['pw_biocarta'] != ['None']:
        biocarta_genes=[]
        biocarta_ann = '.'
        for pw in c['pw_biocarta']:
            cur.execute('select EntrezID from BIOCARTA_EntrezID where Pathway = %s;', pw)
            rows = cur.fetchall()
            biocarta_genes += rows[0][0].split(',')
            if var[:-1].split('\t')[1] in rows[0][0].split(','):
              biocarta.append(pw.replace('_', ' ').title())
        if biocarta != []:
            biocarta_ann = ', '.join(biocarta)
            pathway_count += 1
        if c["pathway_interactors"] == ['y']:
            biocarta_int_ann = '.'
            if biocarta_genes != []:
                cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(biocarta_genes)))
                rows = cur.fetchall()
                biocarta_genes_names = list(set([i[0] for i in rows]))
                biocarta_int = [i for i in biocarta_genes_names if i in dp_interaction_partners]
                if biocarta_int != []:
                    biocarta_int_ann = ';'.join(biocarta_int)
            line += '\t' + '\t'.join([biocarta_ann,biocarta_int_ann])
        else:
            line += '\t' + biocarta_ann
    
    # Reactome pathway    
    reactome = []
    if c['pw_reactome'] != ['None']:
        reactome_genes=[] 
        reactome_ann = '.'
        for pw in c['pw_reactome']:
            cur.execute('select EntrezID from REACTOME_EntrezID where Pathway = %s;', pw)
            rows = cur.fetchall()
            reactome_genes += rows[0][0].split(',')
            if var[:-1].split('\t')[1] in rows[0][0].split(','):
              reactome.append(pw.replace('_', ' ').title())
        if reactome != []:
            reactome_ann = ', '.join(reactome)
            pathway_count += 1
        if c["pathway_interactors"] == ['y']:
            reactome_int_ann = '.'
            if reactome_genes != []:
                cur.execute("select Gene_Name from gID_gNames where Gene_ID in ( {0} );".format(','.join(reactome_genes))) 
                rows = cur.fetchall()
                reactome_genes_names = list(set([i[0] for i in rows]))
                reactome_int = [i for i in reactome_genes_names if i in dp_interaction_partners]
                if reactome_int != []:
                    reactome_int_ann = ';'.join(reactome_int)
            line += '\t' + '\t'.join([reactome_ann,reactome_int_ann])
        else:
            line += '\t' + reactome_ann
    if c["pathway_filter"] == ['y'] and pathway_count == 0:
        continue
   
    fh_dp.write(line.strip() + '\n')
fh_dp.close()


