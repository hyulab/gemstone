# Title: ddG_comphets.py
# Description: This script is to calculate ddG score for variants (particularly on output from recessive compound heterozygous inheritance model).
# Author: Siwei Chen
# Date created: 18/03/2016
# Date last modified: 29/04/2016
# Python version: 2.7.5

import os
import imp
import argparse
import MySQLdb

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="input file")
args = parser.parse_args()

mydb = MySQLdb.connect(host="", user="", passwd="", db="")
cur = mydb.cursor()

# Create a list of existing preminimized 3D structures
os.system('''ls ../../GeMSTONE-data/ddG/preminimized/ > ../../GeMSTONE-data/ddG/preminimized_all.txt''')
fh_preminimized= open('../../GeMSTONE-data/ddG/preminimized_all.txt').read().split('\n')[:-1]

fh_2ddg = open(args.f).readlines()
fh_ddg = open(args.f[:-3]+'ddg.txt', 'w')
fh_ddg.write(fh_2ddg[0][:-1] + '\tRosetta_ddG\n')
for var in fh_2ddg[1:]:
    consequence = var[:-1].split('\t')[13]
    uniprot = var[:-1].split('\t')[22]
    if consequence != 'missense_variant' or uniprot == '.':
        fh_ddg.write(var[:-1] + '\t.\n')
        continue
    else:
        aa_pos = int(var.split('\t')[21][5:-3])
        cur.execute('''select modelID,uniprot_begin from m3D_pdb_structures_index where uniprot= '%s' and uniprot_begin <= '%d' and uniprot_end >= '%d';''' % (uniprot, aa_pos, aa_pos))
        rows_pdb = cur.fetchall()
        if rows_pdb == ():
            cur.execute('''select modbase_modelID,modpipe_quality_score,target_begin from m3D_modbase_models where uniprot= '%s' and target_begin <= '%d' and target_end >= '%d' and modpipe_quality_score >= 1.1;''' % (
                uniprot, aa_pos, aa_pos))
            rows_modbase = cur.fetchall()
            if rows_modbase == ():
                fh_ddg.write(var[:-1] + '\t.\n')
                continue
            else:
                modelID = rows_modbase[0][0]
                aa_begin = rows_modbase[0][2]
        else:
            modelID = rows_pdb[0][0]
            aa_begin = rows_pdb[0][1]
        # Check if there is an existing preminimized 3D structure for this modelID
        if 'min_cst_0.5.' + modelID + '_0001.pdb' not in fh_preminimized:
            fh_ddg.write(var[:-1] + '\t.\n')
            continue
        else:
            # Create the mut file (using relative position of aa_change)
            # Create a dict of {3_letter_aa:1_letter_aa}
            aa_code = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Glx': 'Z', 'Gly': 'G', 'His': 'H',
                       'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
            aa_change = aa_code[var.split('\t')[21][2:5]] + '\t' + str(aa_pos - aa_begin + 1) + '\t' + aa_code[var.split('\t')[21][-3:]] + '\n'
            open(os.path.dirname(args.f)+'/ddG/' + modelID + '_mutlist.txt', 'w').write('total 1\n1\n' + aa_change + '\n')
            # Run ddg_monomer to calculate ddg
            # (two input files: 1. min_cst_0.5.$modelID_0001.pdb (from preminimize step of $modelID.pdb) 2. $modelID_mutlist.txt)
            # (one output file: modelID_ddg_predictions.out)
            os.system('''/opt/rosetta_bin_linux_2015.39.58186_bundle/main/source/bin/ddg_monomer.linuxgccrelease \
                  -in:file:s ../../GeMSTONE-data/ddG/preminimized/min_cst_0.5.''' + modelID + '_0001.pdb ' +  '''\
                  -ddg:out {0}/ddG/'''.format(os.path.dirname(args.f)) + modelID + '_ddg_predictions.out ' + '''\
                  -ddg::mut_file {0}/ddG/'''.format(os.path.dirname(args.f)) + modelID + '_mutlist.txt ' + '''\
                  -ddg:weight_file soft_rep_design \
                  -database /opt/rosetta_bin_linux_2015.39.58186_bundle/main/database/ \
                  -fa_max_dis 9.0 \
                  -ddg::iterations 20 \
                  -ddg::dump_pdbs false \
                  -ignore_unrecognized_res \
                  -ignore_zero_occupancy false \
                  -ddg::local_opt_only true \
                  -ddg::suppress_checkpointing true \
                  -in::file::fullatom -ddg::mean false \
                  -ddg::min true -ddg::min_cst false \
                  -mute all''')

            # Annotate ddg
            if os.path.exists('{0}/ddG/'.format(os.path.dirname(args.f)) + modelID + '_ddg_predictions.out'):
                fh_ddg_out = open('{0}/ddG/'''.format(os.path.dirname(args.f)) + modelID + '_ddg_predictions.out').readlines()
                ddg= fh_ddg_out[-2].split()[2]
                fh_ddg.write(var[:-1] + '\t' + ddg + '\n')
            else:
                fh_ddg.write(var[:-1] + '\t.\n')
fh_ddg.close()

