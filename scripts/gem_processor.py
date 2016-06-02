# Title: gem_processor_latest.py
# Description: This script lives constantly on the server ready to send runs!
# Author: Juan Felipe Beltran, Siwei Chen
# Date created: 18/08/2015
# Date last modified: 24/04/2016
# Python version: 2.7.5

import time
import os
import subprocess as sub
import imp
from itertools import groupby
import shutil
import MySQLdb
from glob import glob

vcftools = 'vcftools'


def sql_pipe(sql, params):
    mydb = MySQLdb.connect(host="", user="", passwd="", db="")
    cursor = mydb.cursor()
    try:
        cursor.execute(sql, params)
        mydb.commit()
    except:
        mydb.rollback()
    result = [row for row in cursor]
    cursor.close()
    mydb.close()
    return result


def bash_pipe(command):
    print
    print "=" * 80
    print command
    print "-" * 80
    os.system(command)
    print "=" * 80


def mail(to, subject, text):
    if to == 'test@gmail.com':
        return
    import smtplib
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText
    gmail_usr = ""
    gmail_pwd = ""
    msg = MIMEMultipart()
    msg['From'] = gmail_usr
    msg['To'] = to
    msg['Subject'] = subject
    msg.attach(MIMEText(text))
    mailServer = smtplib.SMTP("smtp.gmail.com", 587)
    mailServer.ehlo()
    mailServer.starttls()
    mailServer.ehlo()
    mailServer.login(gmail_usr, gmail_pwd)
    mailServer.sendmail(gmail_usr, to, msg.as_string())
    mailServer.close()


def processFile(filename):
    print 'processing ' + filename
    c = imp.load_source('cache', filename).cache

    if c['qid'] in os.listdir('../query/output'):
        shutil.rmtree('../query/output/{0}'.format(c['qid']))
        shutil.rmtree('../query/D3_output/{0}'.format(c['qid']))

    sql_pipe('UPDATE USER_MATCHING SET status="running" WHERE qid=%s', (c['qid'],))

    os.mkdir('../query/output/{0}'.format(c['qid']))
    os.mkdir('../query/output/{0}/ddG'.format(c['qid']))
    os.mkdir('../query/output/{0}/X_files'.format(c['qid']))
    bash_pipe("setfacl -m m:rwx ../query/output/{0}/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/output/{0}/".format(c['qid']))
    bash_pipe("setfacl -m m:rwx ../query/output/{0}/ddG/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/output/{0}/ddG/".format(c['qid']))

    os.mkdir('../query/output/{0}/burden'.format(c['qid']))
    bash_pipe("setfacl -m m:rwx ../query/output/{0}/burden/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/output/{0}/burden/".format(c['qid']))
    bash_pipe("setfacl -m m:rwx ../query/output/{0}/X_files/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/output/{0}/X_files/".format(c['qid']))
    os.mkdir('../query/D3_output/{0}'.format(c['qid']))
    bash_pipe("setfacl -m m:rwx ../query/D3_output/{0}/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/".format(c['qid']))
    os.mkdir('../query/D3_output/{0}/unfiltered/'.format(c['qid']))
    bash_pipe(
        "setfacl -m m:rwx ../query/D3_output/{0}/unfiltered/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/unfiltered/".format(c['qid']))
    os.mkdir('../query/D3_output/{0}/filtered/'.format(c['qid']))
    bash_pipe(
        "setfacl -m m:rwx ../query/D3_output/{0}/filtered/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/filtered/".format(c['qid']))
    os.mkdir('../query/D3_output/{0}/unfiltered/general/'.format(c['qid']))
    bash_pipe(
        "setfacl -m m:rwx ../query/D3_output/{0}/unfiltered/general/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/unfiltered/general/".format(c['qid']))
    os.mkdir('../query/D3_output/{0}/filtered/general/'.format(c['qid']))
    bash_pipe(
        "setfacl -m m:rwx ../query/D3_output/{0}/filtered/general/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/filtered/general/".format(c['qid']))
    os.mkdir(
        '../query/D3_output/{0}/unfiltered/individual/'.format(c['qid']))
    bash_pipe(
        "setfacl -m m:rwx ../query/D3_output/{0}/unfiltered/individual/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/unfiltered/individual/".format(c['qid']))
    os.mkdir(
        '../query/D3_output/{0}/filtered/individual/'.format(c['qid']))
    bash_pipe(
        "setfacl -m m:rwx ../query/D3_output/{0}/filtered/individual/".format(c['qid']))
    bash_pipe(
        "setfacl -m u:apache:rwx ../query/D3_output/{0}/filtered/individual/".format(c['qid']))

    ped_contents = [line.split('\t') for line in open(
        '../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n')]

    # Uncompress input vcf if .vcf.gz
    if c['vcf_filename'][-7:] == '.vcf.gz':
        bash_pipe("""gunzip ../query/data/{0}/{1}""".format(c['qid'], c['vcf_filename']))
        c['vcf_filename'] = c['vcf_filename'][:-3]

    # First convert chromosome name from chr1 to 1
    bash_pipe("""sed 's/^chr//g' ../query/data/{0}/{1} > ../query/output/{0}/{0}.vcf""".format(c['qid'], c['vcf_filename']))
    # Select chrs and regions for different inheritance models
    if c['inheritance'] in ['Dominant', 'Recessive', 'Rec_compound', 'No']:
        if c['nonPAR'] == 'Y': 
            bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.vcf \
                      --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 \
                      --recode --stdout | \
                      bcftools convert --haploid2diploid \
                      > ../query/output/{0}/{0}.chr.hap.vcf""".format(c['qid'], c['vcf_filename']))

            bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.vcf \
                      --positions ../../GeMSTONE-data/PAR.pos.{2}.txt \
                      --recode --stdout | \
                      bcftools convert --haploid2diploid | \
                      sed '/^#/d' \
                      >> ../query/output/{0}/{0}.chr.hap.vcf""".format(c['qid'], c['vcf_filename'], c['build']))
        if c['nonPAR'] != 'Y':
            bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.vcf \
                      --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 \
                      --chr X --chr Y \
                      --recode --stdout | \
                      bcftools convert --haploid2diploid \
                      > ../query/output/{0}/{0}.chr.hap.vcf""".format(c['qid'], c['vcf_filename']))

    if c['inheritance'] in ['XR', 'XD']:
        bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.vcf \
                  --chr X --exclude-positions ../../GeMSTONE-data/PAR.pos.{2}.txt \
                  --recode --stdout | \
                  bcftools convert --haploid2diploid \
                  > ../query/output/{0}/{0}.chr.hap.vcf""".format(c['qid'], c['vcf_filename'], c['build']))

    if c['inheritance'] == 'Y':
        bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.vcf \
                  --chr Y --exclude-positions ../../GeMSTONE-data/PAR.pos.{2}.txt \
                  --recode --stdout | \
                  bcftools convert --haploid2diploid \
                  > ../query/output/{0}/{0}.chr.hap.vcf""".format(c['qid'], c['vcf_filename'], c['build']))

    # Normalization
    if c['build'] == 'GRCh38.79':
        bash_pipe("""sed 's/ID=AD,Number=./ID=AD,Number=R/' ../query/output/{0}/{0}.chr.hap.vcf | \
                  sed 's/chr//g' | \
                  vt decompose -s - | \
                  vt decompose_blocksub -a - | \
                  vt normalize -n -r ../../GeMSTONE-data/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa - | \
                  vt uniq - > ../query/output/{0}/{0}.alts.vcf """.format(c['qid'], c['vcf_filename']))
    if c['build'] == 'GRCh37.75':
        bash_pipe("""sed 's/ID=AD,Number=./ID=AD,Number=R/' ../query/output/{0}/{0}.chr.hap.vcf | \
                  sed 's/chr//g' | \
                  vt decompose -s - | \
                  vt decompose_blocksub -a - | \
                  vt normalize -n -r ../../GeMSTONE-data/human_g1k_v37.fasta - | \
                  vt uniq - > ../query/output/{0}/{0}.alts.vcf """.format(c['qid'], c['vcf_filename']))

    bash_pipe("python helper_scripts/convert_hap2dip.py \
              -i ../query/output/{0}/{0}.alts.vcf \
              -o ../query/output/{0}/{0}.alts.vcf".format(c['qid']))

    os.remove("../query/output/{0}/{0}.vcf".format(c['qid']))
    os.remove("../query/output/{0}/{0}.chr.hap.vcf".format(c['qid']))

    # Per-site visualization stats before all filters
    bash_pipe(
        """python helper_scripts/stats_persite_BF.py \
        -v ../query/output/{0}/{0}.alts.vcf -c {1}""".format(c['qid'], filename))

    # QC
    bash_pipe("""{0} --vcf ../query/output/{1}/{1}.alts.vcf --minQ {2} {3} --minGQ {4} --minDP {5} """
              .format(vcftools, c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
              """--recode --out ../query/output/{0}/{0}.alts.filter""".format(c['qid'], c['vcf_filename']))

    # Taking individuals from the pedigree file
    bash_pipe("""awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $2}' """ +
              """../query/data/{0}/{1} """ .format(c['qid'], c['ped_filename']) +
              """> ../query/output/{0}/{0}.indv.txt""".format(c['qid'], c['ped_filename']))

    ped_contents = [line.split('\t') for line in open(
        '../query/data/{0}/{1}'.format(c['qid'], c['ped_filename'])).read().strip().split('\n')]
    individual_info = dict([(vector[1], vector[2:])
                            for vector in ped_contents])
    family_members = dict([(key, list([g[1] for g in group])) for key, group in groupby(
        sorted([(p[0], p[1]) for p in ped_contents]), lambda x: x[0])])

    # Recessive compound hete inheritance model
    # Generate input for GEMINI
    if c['inheritance'] == 'Rec_compound':
        joint_control = 'N'
        # Predict impact of the mutation and add it to the vcf under INFO > ANN
        bash_pipe("""java -Xmx4g -jar /opt/snpEff/snpEff.jar eff \
                  -v {2} -canon \
                  -q ../query/output/{0}/{0}.alts.filter.recode.vcf \
                  > ../query/output/{0}/{0}.alts.filter.snpeff.vcf""".format(c['qid'], c['vcf_filename'], c['build']))
        bash_pipe("""grep ^# ../query/output/{0}/{0}.alts.filter.snpeff.vcf """.format(c['qid'], c['vcf_filename']) +
                  """> ../query/output/{0}/{0}.alts.filter.eff.vcf""".format(c['qid'], c['vcf_filename']))
        # Resctrict consequences to default selection as comphets only interpretable for mutations changing aa sequences
        # (synonymous and other non-coding mutations are restricted)
        # consequences_restrict = ['frameshift_variant', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion', 'inframe_deletion', 'inframe_insertion', 'start_lost', 'stop_gained', 'stop_lost', 'missense_variant', 'rare_amino_acid_variant', 'exon_loss_variant', 'intron_gain_variant', 'splice_acceptor_variant', 'splice_donor_variant']
        # Further filter for the variant consequences chosen by the user
        bash_pipe("""egrep "{2}" ../query/output/{0}/{0}.alts.filter.snpeff.vcf """.format(c['qid'], c['vcf_filename'], '|'.join(c['consequences'])) +
                  """>> ../query/output/{0}/{0}.alts.filter.eff.vcf""".format(c['qid']))

        bash_pipe("python helper_scripts/Rec_compound.py -para {0}".format(filename))
        bash_pipe("python helper_scripts/af_sc_ANN_comphets0417.py -c {0}".format(filename))

    else:

        # Controls in joint called vcf according to inheritance model
        # no controls filtering at this step for recessive compound hete -- controled separately within GEMINI
        # -- just change file name with .spc. for convenience
        samples = sub.Popen(["""bcftools query -l ../query/output/{0}/{0}.alts.filter.recode.vcf""".format(c['qid'])],
                            shell=True, stdout=sub.PIPE, stderr=sub.PIPE).communicate()[0].strip().split("\n")
        ped_samples = [line.strip() for line in open("../query/output/{0}/{0}.indv.txt".format(c['qid'])).readlines()]
        if [i for i in samples if i not in ped_samples] == []:
            joint_control = 'N'
        else:
            joint_control = 'Y'
            if c['inheritance'] in ['Dominant', 'XD', 'No']:
                bash_pipe("""{0} """.format(vcftools) +
                          """--vcf ../query/output/{0}/{0}.alts.filter.recode.vcf """.format(c['qid'], c['vcf_filename']) +
                          """--remove ../query/output/{0}/{0}.indv.txt --recode --stdout | """.format(c['qid'], c['vcf_filename']) +
                          """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | """ +
                          """sed '/1[\|\/][\.01]/d' | sed '/[\.0][\|\/]1/d' | """ +
                          """awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $1,$2,$3,$4}' | sort """ +
                          """> ../query/output/{0}/{0}.spc.pos.txt""".format(c['qid'], c['vcf_filename']))
            if c['inheritance'] in ['Recessive', 'XR', 'Y']:
                bash_pipe("""{0} """.format(vcftools) +
                          """--vcf ../query/output/{0}/{0}.alts.filter.recode.vcf """.format(c['qid'], c['vcf_filename']) +
                          """--remove ../query/output/{0}/{0}.indv.txt --recode --stdout | """.format(c['qid'], c['vcf_filename']) +
                          """bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | """ +
                          """sed '/1[\|\/]1/d' | """ +
                          """awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $1,$2,$3,$4}' | sort """ +
                          """> ../query/output/{0}/{0}.spc.pos.txt""".format(c['qid'], c['vcf_filename']))

        ctrl_report = open("../query/output/{0}/{0}.ctrl.report.txt".format(c['qid']), 'w')
        ctrl_report.write("Samples: " + ','.join(samples) + '\nJoint-controls: ' + joint_control + '; ' + ','.join([i for i in samples if i not in ped_samples]) + '\n')
        ctrl_report.write("Separate controls: " + c['control_filename'] + '\n')
        ctrl_report.close()
        # Remove variants in additional controls
        if c['control_filename'] == '':
            # Use CHROM,POS,REF,ALT pairs to create a new VCF file containing only desired
            if joint_control == 'Y':
                bash_pipe("""grep ^# ../query/output/{0}/{0}.alts.filter.recode.vcf > ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))
                bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
                          ../query/output/{0}/{0}.spc.pos.txt ../query/output/{0}/{0}.alts.filter.recode.vcf \
                          >> ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))
            else:
                bash_pipe("""cp ../query/output/{0}/{0}.alts.filter.recode.vcf ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))

        elif c['control_filename'][-4:] == '.txt':
            if joint_control == 'Y':
                bash_pipe("""comm -13 ../query/data/{0}/{1} ../query/output/{0}/{0}.spc.pos.txt \
                          > ../query/output/{0}/{0}.spc.pos.tmp.txt""".format(c['qid'], c['control_filename']))
                bash_pipe("""mv ../query/output/{0}/{0}.spc.pos.tmp.txt ../query/output/{0}/{0}.spc.pos.txt""".format(c['qid'], c['vcf_filename']))

                bash_pipe("""grep ^# ../query/output/{0}/{0}.alts.filter.recode.vcf > ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))
                bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
                          ../query/output/{0}/{0}.spc.pos.txt ../query/output/{0}/{0}.alts.filter.recode.vcf \
                          >> ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))

            else:
                bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} !(($1$2$4$5) in vals) {{print}}' \
                          ../query/data/{0}/{1} ../query/output/{0}/{0}.alts.filter.recode.vcf \
                          > ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['control_filename']))

        else:
            if c['control_filename'][-7:] == '.vcf.gz':
                bash_pipe("""gunzip ../query/data/{0}/{1}""".format(c['qid'], c['control_filename']))
                c['control_filename'] = c['control_filename'][:-3]
            if c['control_filename'][-4:] == '.vcf':
                if c['build'] == 'GRCh38.79':
                    bash_pipe("""sed 's/ID=AD,Number=./ID=AD,Number=R/' ../query/data/{0}/{1} | \
                              sed 's/chr//g' | \
                              vt decompose -s - | \
                              vt normalize -n -r ../../GeMSTONE-data/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa - | \
                              vt uniq - > ../query/output/{0}/{0}.ctrl.alts.vcf""".format(c['qid'], c['control_filename']))
                if c['build'] == 'GRCh37.75':
                    bash_pipe("""sed 's/ID=AD,Number=./ID=AD,Number=R/' ../query/data/{0}/{1} | \
                              sed 's/chr//g' | \
                              vt decompose -s - | \
                              vt normalize -n -r ../../GeMSTONE-data/human_g1k_v37.fasta - | \
                              vt uniq - > ../query/output/{0}/{0}.ctrl.alts.vcf""".format(c['qid'], c['control_filename']))

                bash_pipe("python helper_scripts/convert_hap2dip.py \
                  -i ../query/output/{0}/{0}.ctrl.alts.vcf \
                  -o ../query/output/{0}/{0}.ctrl.alts.vcf".format(c['qid']))

                # Screen genotypes for of at least one hete/homo for Dom/Rec
                if c['inheritance'] == 'Dominant':
                    bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.ctrl.alts.vcf  --minQ {1} {2} --minGQ {3} --minDP {4}  \
                              --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 \
                              --chr X --chr Y """.format(c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
                              """--recode --stdout | \
                              bcftools convert --haploid2diploid |\
                              bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | \
                              egrep '1[\|\/][\.01]|[\.0][\|\/]1' | \
                              sed 's/chr//g' |\
                              awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort \
                              > ../query/output/{0}/{0}.ctrl.pos.txt""".format(c['qid'], c['control_filename']))
                if c['inheritance'] == 'XD':
                    bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.ctrl.alts.vcf  --minQ {1} {2} --minGQ {3} --minDP {4} \
                              --chr X """.format(c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
                              """--recode --stdout | \
                              bcftools convert --haploid2diploid | \
                              bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | \
                              egrep '1[\|\/][\.01]|[\.0][\|\/]1' | \
                              sed 's/chr//g' | \
                              awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort \
                              > ../query/output/{0}/{0}.ctrl.pos.txt""".format(c['qid'], c['control_filename']))
                if c['inheritance'] == 'Recessive':
                    bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.ctrl.alts.vcf  --minQ {1} {2} --minGQ {3} --minDP {4} \
                              --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 \
                              --chr X --chr Y """.format(c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
                              """--recode --stdout | \
                              bcftools convert --haploid2diploid | \
                              bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | \
                              egrep '1[\|\/]1' |\
                              sed 's/chr//g' |\
                              awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort \
                              > ../query/output/{0}/{0}.ctrl.pos.txt""".format(c['qid'], c['control_filename']))
                if c['inheritance'] == 'XR':
                    bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.ctrl.alts.vcf  --minQ {1} {2} --minGQ {3} --minDP {4} \
                              --chr X """.format(c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
                              """--recode --stdout |\
                              bcftools convert --haploid2diploid |\
                              bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |\
                              egrep '1[\|\/]1' |\
                              sed 's/chr//g' |\
                              awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort \
                              > ../query/output/{0}/{0}.ctrl.pos.txt""".format(c['qid'], c['control_filename']))
                if c['inheritance'] == 'Y':
                    bash_pipe("""vcftools --vcf ../query/output/{0}/{0}.ctrl.alts.vcf  --minQ {1} {2} --minGQ {3} --minDP {4} \
                              --chr Y """.format(c['qid'], c['qual'], c['flag'], c['gq'], c['dp']) +
                              """--recode --stdout |\
                              bcftools convert --haploid2diploid |\
                              bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |\
                              egrep '1[\|\/]1' |\
                              sed 's/chr//g' |\
                              awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; {{print $1,$2,$3,$4}}' | sort \
                              > ../query/output/{0}/{0}.ctrl.pos.txt""".format(c['qid'], c['control_filename']))

                if joint_control == 'Y':
                    bash_pipe("""comm -13 ../query/output/{0}/{0}.ctrl.pos.txt ../query/output/{0}/{0}.spc.pos.txt \
                              > ../query/output/{0}/{0}.spc.pos.tmp.txt""".format(c['qid'], c['vcf_filename']))
                    bash_pipe("""mv ../query/output/{0}/{0}.spc.pos.tmp.txt ../query/output/{0}/{0}.spc.pos.txt""".format(c['qid'], c['vcf_filename']))
                    bash_pipe("""grep ^# ../query/output/{0}/{0}.alts.filter.recode.vcf > ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))
                    bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} ($1$2$4$5) in vals' \
                              ../query/output/{0}/{0}.spc.pos.txt ../query/output/{0}/{0}.alts.filter.recode.vcf \
                              >> ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))
                else:
                    bash_pipe("""awk 'NR==FNR {{vals[$1$2$3$4];next}} !(($1$2$4$5) in vals) {{print}}' \
                              ../query/output/{0}/{0}.ctrl.pos.txt ../query/output/{0}/{0}.alts.filter.recode.vcf \
                              > ../query/output/{0}/{0}.alts.filter.spc.vcf""".format(c['qid'], c['vcf_filename']))


        # Predict impact of the mutation and add it to the vcf under INFO > ANN
        bash_pipe("""java -Xmx4g -jar /opt/snpEff/snpEff.jar eff \
                  -v {2} -canon \
                  -q ../query/output/{0}/{0}.alts.filter.spc.vcf \
                  > ../query/output/{0}/{0}.alts.filter.spc.snpeff.vcf""".format(c['qid'], c['vcf_filename'], c['build']))

        bash_pipe("""grep ^# ../query/output/{0}/{0}.alts.filter.spc.snpeff.vcf """.format(c['qid'], c['vcf_filename']) +
                  """> ../query/output/{0}/{0}.alts.filter.spc.eff.vcf""".format(c['qid'], c['vcf_filename']))

        # Filter for the variant consequences chosen by the user
        bash_pipe("""egrep "{2}" ../query/output/{0}/{0}.alts.filter.spc.snpeff.vcf """.format(c['qid'], c['vcf_filename'], '|'.join(c['consequences'])) +
                  """>> ../query/output/{0}/{0}.alts.filter.spc.eff.vcf""".format(c['qid']))

        # Population code for ethnicity-specific AF filter
        pop_code = ['esp_AA', 'esp_ALL', 'esp_EA', 'exac_AFR', 'exac_ALL', 'exac_AMR',
                    'exac_Adj', 'exac_EAS', 'exac_FIN', 'exac_NFE', 'exac_OTH', 'exac_SAS',
                    'kg_ALL', 'kg_AFR', 'kg_AMR', 'kg_EAS', 'kg_EUR', 'kg_SAS', 'tagc_AJ',
                    'exac', 'esp', 'kg', 'tagc']

        for f_id in family_members:
            pop_ids = [individual_info[member][-1]
                       for member in family_members[f_id]]
            pop_ids = ','.join(pop_ids).split(",")
            pop = list(set([p if p in pop_code else 'exac_ALL' for p in pop_ids]))

            # Co-segregation analysis
            bash_pipe(
                "python helper_scripts/inheritance.py -c {0} -f {1}".format(filename, f_id))
            # AF filter and variant annotations (w / w/o filtering according to user's choice)
            if c['build'] == 'GRCh37.75':
                bash_pipe("python helper_scripts/af_sc_ANN.py -c {0} -p {1}".format(filename, ','.join(pop)))
            else:
                bash_pipe("python helper_scripts/af_sc_ANN.hg38.py -c {0} -p {1}".format(filename, ','.join(pop)))
            # Individual output file of each family and sporadic sample
            if len(family_members[f_id]) > 1:
                bash_pipe('mv ../query/output/{0}/{0}.alts.filter.spc.eff.cs.ANN.txt ../query/output/{0}/Family_{2}_co-segregating_variants.txt'.format(
                    c['qid'], c['vcf_filename'], f_id))
            else:
                bash_pipe('mv ../query/output/{0}/{0}.alts.filter.spc.eff.cs.ANN.txt ../query/output/{0}/Sporadic_sample_{2}_variants.txt'.format(
                    c['qid'], c['vcf_filename'], family_members[f_id][0]))

        bash_pipe("python helper_scripts/cross_family.py -c {0}".format(filename))

    # ddG calculation
    if 'ddG' in c['ann_fip_Rosetta']:
        if c['inheritance'] != 'Rec_compound':
            bash_pipe(
                "python helper_scripts/ddG.py -f ../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid']))
        else:
            bash_pipe(
                "python helper_scripts/ddG_comphets.py -f ../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid']))
            
        bash_pipe('mv ../query/output/{0}/Cross-sample_all_variants.ddg.txt \
                  ../query/output/{0}/Cross-sample_all_variants.txt'.format(c['qid'], c['vcf_filename']))
    # Functional prediction annotation (w / w/o filtering according to user's choice)
    if c['ann_fip'] != ['None']:
        bash_pipe(
            "python helper_scripts/functional_scores_dbNSFP.py -c {1} -f ../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid'], filename))
        bash_pipe('mv ../query/output/{0}/Cross-sample_all_variants.fip.txt \
                  ../query/output/{0}/Cross-sample_all_variants.txt'.format(c['qid'], c['vcf_filename']))

    # Generate gene file from variant file
    fh2gene = open("../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid'])).readlines()
    fh_gene = open("../query/output/{0}/Cross-sample_all_genes.txt".format(c['qid']), 'w')
    genes = (fh2gene[0].strip().split("\t")[9:12] + [fh2gene[0].strip().split("\t")[13], fh2gene[0].strip().split("\t")[17]]) if c['inheritance'] != 'Rec_compound' else \
            (fh2gene[0].strip().split("\t")[15:19] + [fh2gene[0].strip().split("\t")[22]])
    fh_gene.write("\t".join(genes) + '\n')
    identifiers5 = set()
    for line in fh2gene[1:]:
        genes = (line.strip().split("\t")[9:12] + [line.strip().split("\t")[13], line.strip().split("\t")[17]]) if c['inheritance'] != 'Rec_compound' else \
                (line.strip().split("\t")[15:19] + [line.strip().split("\t")[22]])
        identifiers5.add(tuple(genes))
    for gene in sorted(list(identifiers5)):
        fh_gene.write("\t".join(gene) + '\n')
    fh_gene.close()

    fh2gene = open("../query/output/{0}/Cross-sample_recurrent_variants.txt".format(c['qid'])).readlines()
    fh_gene = open("../query/output/{0}/Cross-sample_recurrent_genes.txt".format(c['qid']), 'w')
    genes = (fh2gene[0].strip().split("\t")[9:12] + [fh2gene[0].strip().split("\t")[13], fh2gene[0].strip().split("\t")[17]]) if c['inheritance'] != 'Rec_compound' else \
            (fh2gene[0].strip().split("\t")[15:19] + [fh2gene[0].strip().split("\t")[22]])
    fh_gene.write("\t".join(genes) + '\n')
    identifiers5 = set()
    for line in fh2gene[1:]:
        genes = (line.strip().split("\t")[9:12] + [line.strip().split("\t")[13], line.strip().split("\t")[17]]) if c['inheritance'] != 'Rec_compound' else \
                (line.strip().split("\t")[15:19] + [line.strip().split("\t")[22]])
        identifiers5.add(tuple(genes))
    for gene in sorted(list(identifiers5)):
        fh_gene.write("\t".join(gene) + '\n')
    fh_gene.close()

    # Gene level annotations(w / w/o filtering according to user's choice)
    # user upload gene list of interest
    if c['gl_user'] != 'None':
        bash_pipe("""python helper_scripts/genes_short.py \
                  -c {0} \
                  -f ../query/output/{1}/Cross-sample_all_genes.txt \
                  -l {2} \
                  -o ../query/output/{1}/Cross-sample_all_genes.ann.txt""".format(filename, c['qid'], c['gl_user']))
    else:
        bash_pipe("""python helper_scripts/genes_short.py \
                  -c {0} \
                  -f ../query/output/{1}/Cross-sample_all_genes.txt \
                  -o ../query/output/{1}/Cross-sample_all_genes.ann.txt""".format(filename, c['qid']))

    bash_pipe("""mv ../query/output/{1}/Cross-sample_all_genes.ann.txt ../query/output/{1}/Cross-sample_all_genes.txt""".format(filename, c['qid']))

    # Pathway enrichment using the gene list
    if c['pea'] != ['None']:
        open("../query/output/{0}/Cross-sample_all_genes2enrichment.txt".format(c['qid']), 'w').write(
            "\n".join(line.strip().split("\t")[1] for
                      line in open("../query/output/{0}/Cross-sample_all_genes.txt".format(c['qid'])).readlines()[1:]))
        open("../query/output/{0}/Cross-sample_recurrent_genes2enrichment.txt".format(c['qid']), 'w').write(
            "\n".join(line.strip().split("\t")[1] for
                      line in open("../query/output/{0}/Cross-sample_recurrent_genes.txt".format(c['qid'])).readlines()[1:]))

        # Convert to entrez gene list for input (good practice when using the script as a outside the pipeline)
        converted_background_list = {'background_gene_list_homo_sapiens': '../../GeMSTONE-data/Homo_sapiens_Biomart_entrez_coding.txt',
                                     'background_gene_list_vcf': '../query/data/{0}/snpeff_entrez_genes.txt'.format(c['qid']),
                                     'background_gene_list_custom': ''}
        if c['background_gene_list'] == 'background_gene_list_vcf':
            HGNC_entrez = dict(line.strip().split("\t") for line in open("../../GeMSTONE-data/Homo_sapiens_Biomart_entrezHGNC.txt").readlines() if len(line.strip().split("\t")) > 1)
            open("../query/data/{0}/snpeff_entrez_genes.txt".format(c['qid']), 'w').write("\n".join(
                [HGNC_entrez[sgene] for sgene in [line.strip().split("\t")[0] for line in open('snpEff_genes.txt').readlines()]
                 if sgene in HGNC_entrez]))
        if c['background_gene_list'] == 'background_gene_list_custom':
            converted_background_list['background_gene_list_custom'] = "../query/data/{0}/.pathway_background.txt".format(c['qid'])
            if not open(c['background_user']).readlines()[0].strip().isdigit():
                HGNC_entrez = dict(line.strip().split("\t") for line in open("../../GeMSTONE-data/Homo_sapiens_Biomart_entrezHGNC.txt").readlines() if len(line.strip().split("\t")) > 1)
                open("../query/data/{0}/.pathway_background.converted.txt".format(c['qid']), 'w').write("\n".join(
                    [HGNC_entrez[sgene] for sgene in [line.strip().split("\t")[0] for line in open("../query/data/{0}/.pathway_background.converted.txt".format(c['qid'])).readlines()] if sgene in HGNC_entrez]))
                converted_background_list['background_gene_list_custom'] = "../query/data/{0}/.pathway_background.converted.txt".format(c['qid'])

        for gene_list in ['../query/output/{0}/Cross-sample_all_genes2enrichment.txt'.format(c['qid']),
                          '../query/output/{0}/Cross-sample_recurrent_genes2enrichment.txt'.format(c['qid'])]:
            bash_pipe("""python helper_scripts/pathway_enrichment.py \
                      -c {0} \
                      -out ../query/output/{1}/Pathway_enrichment_{2}_genes.txt \
                      -target {3} \
                      -background {4}""".format(filename,
                                                c['qid'],
                                                gene_list.split("/")[-1].split("_")[1],
                                                gene_list,
                                                converted_background_list[c['background_gene_list']]))
    # Gene expression annotation
    if c['expression'] != {}:
        bash_pipe("""python helper_scripts/gene_expression.py -c {0}""".format(filename))
    # Gene damaging index annotation
    if c['gdi'] != ['None']:
        bash_pipe("""python helper_scripts/gdi.py -c {0} -f ../query/output/{1}/Cross-sample_all_genes.txt""".format(filename, c['qid']))
        bash_pipe("""mv ../query/output/{0}/Cross-sample_all_genes.gdi.txt ../query/output/{0}/Cross-sample_all_genes.txt""".format(c['qid']))
    # Gene rvis score annotation
    if c['rvis'] == 'Y':
        bash_pipe("""python helper_scripts/rvis.py -f ../query/output/{0}/Cross-sample_all_genes.txt""".format(c['qid']))
        bash_pipe("""mv ../query/output/{0}/Cross-sample_all_genes.rvis.txt ../query/output/{0}/Cross-sample_all_genes.txt""".format(c['qid']))
    # Gene burden test
    if c['burden'] != ['None']:
        bash_pipe("""python helper_scripts/burden_test.py -c {0}""".format(filename))

    # Gpdate all_ and recurrent_variants from all_genes
    # Gubset recurrent_gene from all_gene after ALL gene annotations
    # Update Family/Sporadic single files based on the all_variant file
    if c['inheritance'] != 'Rec_compound':
        bash_pipe("""awk 'NR==FNR {{vals[$1];next}} ($10) in vals' \
                  ../query/output/{0}/Cross-sample_all_genes.txt ../query/output/{0}/Cross-sample_all_variants.txt \
                  > ../query/output/{0}/Cross-sample_all_variants.functional.txt""".format(c['qid']))
        bash_pipe("""awk 'NR==FNR {{vals[$1$2$4$5];next}} ($1$2$4$5) in vals' \
                  ../query/output/{0}/Cross-sample_recurrent_variants.txt ../query/output/{0}/Cross-sample_all_variants.functional.txt \
                  > ../query/output/{0}/Cross-sample_recurrent_variants.functional.txt""".format(c['qid']))
        bash_pipe("""awk 'NR==FNR {{vals[$10];next}} ($1) in vals' \
                  ../query/output/{0}/Cross-sample_recurrent_variants.functional.txt ../query/output/{0}/Cross-sample_all_genes.txt \
                  > ../query/output/{0}/Cross-sample_recurrent_genes.txt""".format(c['qid']))
        sub_results = glob('../query/output/{0}/Family*'.format(c['qid'])) + glob('../query/output/{0}/Sporadic*'.format(c['qid']))
        # update Family/Sporadic single files and move to sub folder
        all_variants = dict([tuple(line.strip().split("\t")[:2] + line.strip().split("\t")[3:5]), line.strip().split("\t")[7:]]
                            for line in open("../query/output/{0}/Cross-sample_all_variants.functional.txt".format(c['qid'])).readlines())
        for fn in sub_results:
            fh = open(fn).readlines()
            of = open(fn, 'w')
            for line in fh:
                if tuple(line.strip().split("\t")[:2] + line.strip().split("\t")[3:5]) in all_variants:
                    of.write("\t".join(line.strip().split("\t")[:7] +
                                       all_variants[tuple(line.strip().split("\t")[:2] + line.strip().split("\t")[3:5])]) + '\n')
            of.close()

    else:
        bash_pipe("""awk 'NR==FNR {{vals[$1];next}} ($16) in vals' \
                  ../query/output/{0}/Cross-sample_all_genes.txt ../query/output/{0}/Cross-sample_all_variants.txt \
                  > ../query/output/{0}/Cross-sample_all_variants.functional.txt""".format(c['qid']))
        bash_pipe("""awk 'NR==FNR {{vals[$1$2$4$5];next}} ($1$2$4$5) in vals' \
                  ../query/output/{0}/Cross-sample_recurrent_variants.txt ../query/output/{0}/Cross-sample_all_variants.functional.txt \
                  > ../query/output/{0}/Cross-sample_recurrent_variants.functional.txt""".format(c['qid']))

        bash_pipe("""awk 'NR==FNR {{vals[$16];next}} ($1) in vals' \
                  ../query/output/{0}/Cross-sample_recurrent_variants.functional.txt ../query/output/{0}/Cross-sample_all_genes.txt \
                  > ../query/output/{0}/Cross-sample_recurrent_genes.txt""".format(c['qid']))

        for f_id in family_members:
            if len(family_members[f_id]) > 1:
                bash_pipe("""awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; ($8 == "{1}" || $8 == "FAMILY_ID") {{print}}' \
                          ../query/output/{0}/Cross-sample_all_variants.functional.txt \
                          > ../query/output/{0}/Family_{1}_co-segregating_variants.txt""".format(c['qid'], f_id))
            else:
                bash_pipe("""awk 'BEGIN {{FS= "\t"; OFS= "\t"}}; ($8 == "{1}" || $8 == "FAMILY_ID") {{print}}' \
                          ../query/output/{0}/Cross-sample_all_variants.functional.txt \
                          > ../query/output/{0}/Sporadic_sample_{1}_co-segregating_variants.txt""".format(c['qid'], f_id))
    os.rename("../query/output/{0}/Cross-sample_all_variants.functional.txt".format(c['qid']), "../query/output/{0}/Cross-sample_all_variants.txt".format(c['qid']))
    os.rename("../query/output/{0}/Cross-sample_recurrent_variants.functional.txt".format(c['qid']), "../query/output/{0}/Cross-sample_recurrent_variants.txt".format(c['qid']))

    # generate filtered per-site visualization stats
    bash_pipe("""grep '^#' ../query/output/{0}/{0}.alts.vcf > ../query/output/{0}/{0}.alts.stats_AF.vcf""".format(c['qid']))
    bash_pipe("""awk 'NR==FNR {{vals[$1$2$4$5];next}} ($1$2$4$5) in vals' \
              ../query/output/{0}/Cross-sample_all_variants.txt ../query/output/{0}/{0}.alts.vcf\
              >> ../query/output/{0}/{0}.alts.stats_AF.vcf""".format(c['qid']))
    bash_pipe(
        """python helper_scripts/stats_persite_AF.py \
        -v ../query/output/{0}/{0}.alts.stats_AF.vcf""".format(c['qid'], filename))

    shutil.move(filename, '../query/history/')
    os.system('mv ../query/data/{0}/00_* ../query/output/{0}'.format(c['qid']))
    # chdir because zip will remember the path of files zipped -- save all subdirs as paths at the end
    os.chdir('../query/data/{0}'.format(c['qid']))
    os.system('''zip ../../output/{0}/Inputs.zip *'''.format(c['qid']))

    os.chdir('../../output/{0}'.format(c['qid']))
    os.system('''zip sub.zip Family_* Sporadic_sample_*''')
    os.system('''zip Results.zip Cross-sample* Pathway* sub.zip''')
    #os.system('''cp 00_* ../../data/{0}'''.format(c['qid']))

    os.rename("{0}.alts.vcf".format(c['qid']), "01_{0}.normalized.vcf".format(c['qid']))
    os.rename("{0}.alts.filter.recode.vcf".format(c['qid']), "02_{0}.normalized.qc.vcf".format(c['qid']))
    if (c['control_filename'] == '' and joint_control == 'N') and c['inheritance'] != 'Rec_compound':
        bash_pipe("mv {0}.alts.filter.spc.snpeff.vcf 03_{0}.normalized.qc.snpeff.vcf".format(c['qid']))
        bash_pipe("mv {0}.alts.filter.spc.eff.vcf 04_{0}.normalized.qc.snpeff.variant_type.vcf".format(c['qid']))

    elif c['inheritance'] == 'Rec_compound':
        bash_pipe("mv {0}.alts.filter.snpeff.vcf 03_{0}.normalized.qc.snpeff.vcf".format(c['qid']))
        bash_pipe("mv {0}.alts.filter.eff.vcf 03_{0}.normalized.qc.snpeff.vcf".format(c['qid']))

    else:
        bash_pipe("mv {0}.alts.filter.spc.vcf 03_{0}.normalized.qc.controlled.vcf".format(c['qid']))
        bash_pipe("mv {0}.alts.filter.spc.snpeff.vcf 04_{0}.normalized.qc.controlled.snpeff.vcf".format(c['qid']))
        bash_pipe("mv {0}.alts.filter.spc.eff.vcf 05_{0}.normalized.qc.controlled.snpeff.variant_type.vcf").format(c['qid'])

    # Clean output files
    exclude_dirs = ['ddG', 'burden', 'X_files']
    selected_files = glob('''0*''') + \
        glob('''Cross-sample*variants.txt''') + \
        glob('''Cross-sample*genes.txt''') + \
        glob('''Family*''') + \
        glob('''Sporadic*''') + \
        glob('''Pathway*''') + \
        glob('''*zip''')
    X_files = [f for f in glob('''*''') if f not in exclude_dirs + selected_files]
    if X_files:
        for fn in X_files:
            os.system("mv {1} X_files".format(c['qid'], fn))
    os.system('''zip All_{1}.zip {2}'''.format(c['qid'], c['proj'], " ".join(selected_files)))
    os.remove('''sub.zip''')
    os.chdir("/data/germline-web/scripts")

    sql_pipe('UPDATE USER_MATCHING SET status="done" WHERE qid=%s', (c['qid'],))
    mail(c['email'], c['proj'] + ' has finished!', 'Your project {0} has now finished running on GeMSTONE. You can access your results by logging in at http://gemstone.yulab.org'.format(c['proj']))


while True:
    jobs = sql_pipe('SELECT qid, project, email FROM USER_MATCHING WHERE status="queued" ORDER BY Time_Created;', None)
    for qid, project, email in jobs:
        print "Working on", qid, project, email
        try:
            processFile('../query/queue/{0}_para.py'.format(qid))
        except Exception, e:
            sql_pipe('UPDATE USER_MATCHING SET status="failed" WHERE qid=%s', (qid,))
    else:
        print "Sleeping..."
        time.sleep(10)
