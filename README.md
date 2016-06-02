# gemstone


1. Directory structure

GeMSTONE-data/: non-MySQL data sources 
GeMSTONE-main/
	scripts/: GeMSTONE source code
	query/
		data/: Inputs (VCF, PED etc. files as uploaded)
		output/: Result files and intermediate files
		D3_output/: Output data for visualization



2. Tools required

|Tool	|Version	|URL |
|-------|-----------|----|
|VT				|2016-01-12	|	http://genome.sph.umich.edu/wiki/Vt|
|VCFtools		|0.1.14		|	https://vcftools.github.io/man_latest.html|
|BCFtools		|1.2		|	https://samtools.github.io/bcftools/bcftools.html|
|SnpEff			|4.1l		|	http://snpeff.sourceforge.net|
|GEMINI			|0.18.0		|	https://gemini.readthedocs.io/en/latest/|
|CrossMap		|v0.2.3		|	http://crossmap.sourceforge.net|
|dbNSFP			|v3.1a		|	https://sites.google.com/site/jpopgen/dbNSFP|
|Rosetta ddG	|3.5		|	https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer|
|PLINK/SEQ		|v0.10		|	https://atgu.mgh.harvard.edu/plinkseq/index.shtml|
|MySQL			|5.7.9		|	https://www.mysql.com|

3. Databases required

+--------------------------+
| Tables_in_germline       |
+--------------------------+
| 1000G_hg19               |
| 1000G_hg38               |
| BIOCARTA_EntrezID        |
| Cancer_Gene_Census       |
| ClinVar                  |
| DOMAINS                  |
| ESP6500                  |
| ESP6500_hg19             |
| ESP6500_hg38             |
| ExAC_hg19                |
| ExAC_hg38                |
| ExAC_r03_coverage        |
| GO_biological_process    |
| GO_cellular_component    |
| GO_molecular_function    |
| Gene_Damage_Index        |
| Gene_InteractionDB       |
| HGMD_disease_gene        |
| HINT_interactors         |
| HSIN_INTERFACE           |
| KEGG_EntrezID            |
| Kegg_human               |
| LIBRARY                  |
| OMIM_std_all             |
| ONTOLOGIES               |
| PHENOTYPES               |
| PROBABILITIES            |
| PROTEIN                  |
| PROVEAN_SIFT             |
| REACTOME_EntrezID        |
| SEQUENCE                 |
| TAGC128_hg19             |
| TAGC128_hg38             |
| UNIPROT_PFAM_DOMAINS     |
| VARIANTS                 |
| WEIGHTS                  |
| dbNSFP                   |
| enst_ensg                |
| gID_Uniprot              |
| gID_gNames               |
| gID_gSymbol              |
| gSymbol_gID_biomart      |
| gerp_elements            |
| hg19_hg38_mapping        |
| hg38_hg19_mapping        |
| m3D_modbase_models       |
| m3D_pdb_structures_index |
| pdbid_chain_uniprot      |
+--------------------------+


*1000G_hg19*
+-------+--------------+------+-----+---------+-------+
| Field | Type         | Null | Key | Default | Extra |
+-------+--------------+------+-----+---------+-------+
| CHROM | varchar(3)   | NO   | PRI |         |       |
| POS   | varchar(15)  | NO   | PRI |         |       |
| REF   | varchar(250) | NO   | PRI |         |       |
| ALT   | varchar(700) | NO   | PRI |         |       |
| AC    | int(10)      | YES  |     | NULL    |       |
| AF    | float(8,5)   | YES  |     | NULL    |       |
| EAS   | float(8,5)   | YES  |     | NULL    |       |
| AMR   | float(8,5)   | YES  |     | NULL    |       |
| AFR   | float(8,5)   | YES  |     | NULL    |       |
| EUR   | float(8,5)   | YES  |     | NULL    |       |
| SAS   | float(8,5)   | YES  |     | NULL    |       |
+-------+--------------+------+-----+---------+-------+


*1000G_hg38*
+-------+--------------+------+-----+---------+-------+
| Field | Type         | Null | Key | Default | Extra |
+-------+--------------+------+-----+---------+-------+
| CHROM | varchar(3)   | NO   | PRI | NULL    |       |
| POS   | varchar(15)  | NO   | PRI | NULL    |       |
| REF   | varchar(250) | NO   | PRI | NULL    |       |
| ALT   | varchar(700) | NO   | PRI | NULL    |       |
| AC    | int(10)      | YES  |     | NULL    |       |
| AF    | float(9,6)   | YES  |     | NULL    |       |
| EAS   | float(9,6)   | YES  |     | NULL    |       |
| AMR   | float(9,6)   | YES  |     | NULL    |       |
| AFR   | float(9,6)   | YES  |     | NULL    |       |
| EUR   | float(9,6)   | YES  |     | NULL    |       |
| SAS   | float(9,6)   | YES  |     | NULL    |       |
+-------+--------------+------+-----+---------+-------+


*BIOCARTA_EntrezID*
+-----------+----------------+------+-----+---------+-------+
| Field     | Type           | Null | Key | Default | Extra |
+-----------+----------------+------+-----+---------+-------+
| Pathway   | varchar(500)   | YES  |     | NULL    |       |
| Link2gsea | varchar(500)   | YES  |     | NULL    |       |
| EntrezID  | varchar(25000) | YES  |     | NULL    |       |
+-----------+----------------+------+-----+---------+-------+


*Cancer_Gene_Census*
+--------------------+--------------+------+-----+---------+-------+
| Field              | Type         | Null | Key | Default | Extra |
+--------------------+--------------+------+-----+---------+-------+
| Tumor_type         | varchar(255) | YES  |     | NULL    |       |
| Gene_Names_Somatic | varchar(500) | YES  |     | NULL    |       |
| Gene_Name_Germline | varchar(255) | YES  |     | NULL    |       |
+--------------------+--------------+------+-----+---------+-------+


*ClinVar*
+----------------+--------------+------+-----+---------+-------+
| Field          | Type         | Null | Key | Default | Extra |
+----------------+--------------+------+-----+---------+-------+
| UniProt        | varchar(20)  | YES  |     | NULL    |       |
| AASub          | varchar(20)  | YES  |     | NULL    |       |
| GRCh37_Variant | varchar(20)  | YES  |     | NULL    |       |
| dbSNP          | varchar(20)  | YES  |     | NULL    |       |
| Disease        | varchar(200) | YES  |     | NULL    |       |
+----------------+--------------+------+-----+---------+-------+


*DOMAINS*
+-----------+----------+------+-----+---------+-------+
| Field     | Type     | Null | Key | Default | Extra |
+-----------+----------+------+-----+---------+-------+
| id        | int(11)  | NO   | MUL | NULL    |       |
| hmm       | char(15) | NO   | MUL | NULL    |       |
| score     | double   | NO   |     | NULL    |       |
| seq_begin | int(11)  | NO   |     | NULL    |       |
| seq_end   | int(11)  | NO   |     | NULL    |       |
| hmm_begin | int(11)  | NO   |     | NULL    |       |
| align     | text     | NO   |     | NULL    |       |
+-----------+----------+------+-----+---------+-------+


*ESP6500*
+--------+--------------+------+-----+---------+-------+
| Field  | Type         | Null | Key | Default | Extra |
+--------+--------------+------+-----+---------+-------+
| CHROM  | varchar(3)   | NO   | PRI |         |       |
| POS    | varchar(10)  | NO   | PRI |         |       |
| REF    | varchar(250) | NO   | PRI |         |       |
| ALT    | varchar(700) | NO   | PRI |         |       |
| AF_EA  | float(7,5)   | YES  |     | NULL    |       |
| AF_AA  | float(7,5)   | YES  |     | NULL    |       |
| AF_All | float(7,5)   | YES  |     | NULL    |       |
+--------+--------------+------+-----+---------+-------+


*ESP6500_hg19*
+--------+--------------+------+-----+---------+-------+
| Field  | Type         | Null | Key | Default | Extra |
+--------+--------------+------+-----+---------+-------+
| CHROM  | varchar(3)   | NO   | PRI | NULL    |       |
| POS    | varchar(15)  | NO   | PRI | NULL    |       |
| REF    | varchar(250) | NO   | PRI | NULL    |       |
| ALT    | varchar(700) | NO   | PRI | NULL    |       |
| AF_EA  | float(9,6)   | YES  |     | NULL    |       |
| AF_AA  | float(9,6)   | YES  |     | NULL    |       |
| AF_All | float(9,6)   | YES  |     | NULL    |       |
+--------+--------------+------+-----+---------+-------+


*ESP6500_hg38*
+--------+--------------+------+-----+---------+-------+
| Field  | Type         | Null | Key | Default | Extra |
+--------+--------------+------+-----+---------+-------+
| CHROM  | varchar(3)   | NO   | PRI | NULL    |       |
| POS    | varchar(15)  | NO   | PRI | NULL    |       |
| REF    | varchar(250) | NO   | PRI | NULL    |       |
| ALT    | varchar(700) | NO   | PRI | NULL    |       |
| AF_EA  | float(9,6)   | YES  |     | NULL    |       |
| AF_AA  | float(9,6)   | YES  |     | NULL    |       |
| AF_All | float(9,6)   | YES  |     | NULL    |       |
+--------+--------------+------+-----+---------+-------+


*ExAC_hg19*
+--------+--------------+------+-----+---------+-------+
| Field  | Type         | Null | Key | Default | Extra |
+--------+--------------+------+-----+---------+-------+
| CHROM  | varchar(3)   | NO   | PRI |         |       |
| POS    | varchar(10)  | NO   | PRI |         |       |
| REF    | varchar(400) | NO   | PRI |         |       |
| ALT    | varchar(700) | NO   | PRI |         |       |
| AC     | int(11)      | YES  |     | NULL    |       |
| AF     | float(9,6)   | YES  |     | NULL    |       |
| AC_AFR | int(11)      | YES  |     | NULL    |       |
| AF_AFR | float(9,6)   | YES  |     | NULL    |       |
| AC_AMR | int(11)      | YES  |     | NULL    |       |
| AF_AMR | float(9,6)   | YES  |     | NULL    |       |
| AC_Adj | int(11)      | YES  |     | NULL    |       |
| AF_Adj | float(9,6)   | YES  |     | NULL    |       |
| AC_EAS | int(11)      | YES  |     | NULL    |       |
| AF_EAS | float(9,6)   | YES  |     | NULL    |       |
| AC_FIN | int(11)      | YES  |     | NULL    |       |
| AF_FIN | float(9,6)   | YES  |     | NULL    |       |
| AC_NFE | int(11)      | YES  |     | NULL    |       |
| AF_NFE | float(9,6)   | YES  |     | NULL    |       |
| AC_OTH | int(11)      | YES  |     | NULL    |       |
| AF_OTH | float(9,6)   | YES  |     | NULL    |       |
| AC_SAS | int(11)      | YES  |     | NULL    |       |
| AF_SAS | float(9,6)   | YES  |     | NULL    |       |
+--------+--------------+------+-----+---------+-------+


*ExAC_hg38*
+--------+--------------+------+-----+---------+-------+
| Field  | Type         | Null | Key | Default | Extra |
+--------+--------------+------+-----+---------+-------+
| CHROM  | varchar(3)   | NO   | PRI | NULL    |       |
| POS    | varchar(10)  | NO   | PRI | NULL    |       |
| REF    | varchar(400) | NO   | PRI | NULL    |       |
| ALT    | varchar(700) | NO   | PRI | NULL    |       |
| AC     | int(11)      | YES  |     | NULL    |       |
| AF     | float(9,6)   | YES  |     | NULL    |       |
| AC_AFR | int(11)      | YES  |     | NULL    |       |
| AF_AFR | float(9,6)   | YES  |     | NULL    |       |
| AC_AMR | int(11)      | YES  |     | NULL    |       |
| AF_AMR | float(9,6)   | YES  |     | NULL    |       |
| AC_Adj | int(11)      | YES  |     | NULL    |       |
| AF_Adj | float(9,6)   | YES  |     | NULL    |       |
| AC_EAS | int(11)      | YES  |     | NULL    |       |
| AF_EAS | float(9,6)   | YES  |     | NULL    |       |
| AC_FIN | int(11)      | YES  |     | NULL    |       |
| AF_FIN | float(9,6)   | YES  |     | NULL    |       |
| AC_NFE | int(11)      | YES  |     | NULL    |       |
| AF_NFE | float(9,6)   | YES  |     | NULL    |       |
| AC_OTH | int(11)      | YES  |     | NULL    |       |
| AF_OTH | float(9,6)   | YES  |     | NULL    |       |
| AC_SAS | int(11)      | YES  |     | NULL    |       |
| AF_SAS | float(9,6)   | YES  |     | NULL    |       |
+--------+--------------+------+-----+---------+-------+


*ExAC_r03_coverage*
+-------------+-------------+------+-----+---------+-------+
| Field       | Type        | Null | Key | Default | Extra |
+-------------+-------------+------+-----+---------+-------+
| CHROM       | varchar(3)  | NO   | PRI |         |       |
| POS         | varchar(10) | NO   | PRI |         |       |
| MEAN        | float(7,5)  | YES  |     | NULL    |       |
| MEDIAN      | float(7,5)  | YES  |     | NULL    |       |
| one         | float(7,5)  | YES  |     | NULL    |       |
| five        | float(7,5)  | YES  |     | NULL    |       |
| ten         | float(7,5)  | YES  |     | NULL    |       |
| fifteen     | float(7,5)  | YES  |     | NULL    |       |
| twenty      | float(7,5)  | YES  |     | NULL    |       |
| twenty_five | float(7,5)  | YES  |     | NULL    |       |
| thirty      | float(7,5)  | YES  |     | NULL    |       |
| fifty       | float(7,5)  | YES  |     | NULL    |       |
| hundred     | float(7,5)  | YES  |     | NULL    |       |
+-------------+-------------+------+-----+---------+-------+


*GO_biological_process*
+----------+----------------+------+-----+---------+-------+
| Field    | Type           | Null | Key | Default | Extra |
+----------+----------------+------+-----+---------+-------+
| BP       | varchar(100)   | YES  |     | NULL    |       |
| LINK     | varchar(500)   | YES  |     | NULL    |       |
| EntrezID | varchar(50000) | YES  |     | NULL    |       |
+----------+----------------+------+-----+---------+-------+


*GO_cellular_component*
+----------+----------------+------+-----+---------+-------+
| Field    | Type           | Null | Key | Default | Extra |
+----------+----------------+------+-----+---------+-------+
| CC       | varchar(100)   | YES  |     | NULL    |       |
| LINK     | varchar(500)   | YES  |     | NULL    |       |
| EntrezID | varchar(50000) | YES  |     | NULL    |       |
+----------+----------------+------+-----+---------+-------+


*GO_molecular_function*
+----------+----------------+------+-----+---------+-------+
| Field    | Type           | Null | Key | Default | Extra |
+----------+----------------+------+-----+---------+-------+
| MF       | varchar(100)   | YES  |     | NULL    |       |
| LINK     | varchar(500)   | YES  |     | NULL    |       |
| EntrezID | varchar(50000) | YES  |     | NULL    |       |
+----------+----------------+------+-----+---------+-------+


*Gene_Damage_Index*
+--------------+-------------+------+-----+---------+-------+
| Field        | Type        | Null | Key | Default | Extra |
+--------------+-------------+------+-----+---------+-------+
| Gene         | varchar(20) | NO   | PRI | NULL    |       |
| GDI          | varchar(20) | YES  |     | NULL    |       |
| GDI_Phred    | varchar(20) | YES  |     | NULL    |       |
| All_diseases | varchar(20) | YES  |     | NULL    |       |
| Mendelian    | varchar(20) | YES  |     | NULL    |       |
| Mendelian_AD | varchar(20) | YES  |     | NULL    |       |
| Mendelian_AR | varchar(20) | YES  |     | NULL    |       |
| PID          | varchar(20) | YES  |     | NULL    |       |
| PID_AD       | varchar(20) | YES  |     | NULL    |       |
| PID_AR       | varchar(20) | YES  |     | NULL    |       |
| Cancer       | varchar(20) | YES  |     | NULL    |       |
| Cancer_AR    | varchar(20) | YES  |     | NULL    |       |
| Cancer_AD    | varchar(20) | YES  |     | NULL    |       |
+--------------+-------------+------+-----+---------+-------+


*Gene_InteractionDB*
+-----------------+--------------+------+-----+---------+-------+
| Field           | Type         | Null | Key | Default | Extra |
+-----------------+--------------+------+-----+---------+-------+
| Gene_name       | varchar(100) | YES  |     | NULL    |       |
| Ensembl_gene    | varchar(100) | YES  |     | NULL    |       |
| Entrez_gene_id  | varchar(100) | YES  |     | NULL    |       |
| IntAct          | longtext     | YES  |     | NULL    |       |
| BioGRID         | longtext     | YES  |     | NULL    |       |
| ConsensusPathDB | longtext     | YES  |     | NULL    |       |
+-----------------+--------------+------+-----+---------+-------+


*HGMD_disease_gene*
+-----------+---------------+------+-----+---------+-------+
| Field     | Type          | Null | Key | Default | Extra |
+-----------+---------------+------+-----+---------+-------+
| Disease   | varchar(1000) | YES  |     | NULL    |       |
| Gene_ID   | varchar(20)   | YES  |     | NULL    |       |
| Gene_Name | varchar(100)  | YES  |     | NULL    |       |
+-----------+---------------+------+-----+---------+-------+


*HINT_interactors*
+--------+-------------+------+-----+---------+-------+
| Field  | Type        | Null | Key | Default | Extra |
+--------+-------------+------+-----+---------+-------+
| id_A   | varchar(40) | YES  |     | NULL    |       |
| id_B   | varchar(40) | YES  |     | NULL    |       |
| Gene_A | varchar(40) | YES  |     | NULL    |       |
| Gene_B | varchar(40) | YES  |     | NULL    |       |
+--------+-------------+------+-----+---------+-------+


*HSIN_INTERFACE*
+----------------------+-------------+------+-----+---------+-------+
| Field                | Type        | Null | Key | Default | Extra |
+----------------------+-------------+------+-----+---------+-------+
| ENTREZ_ID            | varchar(20) | YES  |     | NULL    |       |
| UNIPROT              | varchar(20) | NO   | PRI | NULL    |       |
| PFAM                 | varchar(20) | YES  |     | NULL    |       |
| START                | int(11)     | NO   | PRI | NULL    |       |
| END                  | int(11)     | NO   | PRI | NULL    |       |
| INTERFACE_PREDICTION | varchar(50) | YES  |     | NULL    |       |
+----------------------+-------------+------+-----+---------+-------+


*KEGG_EntrezID*
+-----------+----------------+------+-----+---------+-------+
| Field     | Type           | Null | Key | Default | Extra |
+-----------+----------------+------+-----+---------+-------+
| Pathway   | varchar(500)   | YES  |     | NULL    |       |
| Link2gsea | varchar(500)   | YES  |     | NULL    |       |
| EntrezID  | varchar(25000) | YES  |     | NULL    |       |
+-----------+----------------+------+-----+---------+-------+


*Kegg_human*
+-----------+--------------+------+-----+---------+-------+
| Field     | Type         | Null | Key | Default | Extra |
+-----------+--------------+------+-----+---------+-------+
| Kegg_ID   | varchar(20)  | YES  |     | NULL    |       |
| Pathway   | varchar(100) | YES  |     | NULL    |       |
| Gene_Name | varchar(40)  | YES  |     | NULL    |       |
| Gene_ID   | varchar(20)  | YES  |     | NULL    |       |
+-----------+--------------+------+-----+---------+-------+


*LIBRARY*
+-------------+----------+------+-----+---------+-------+
| Field       | Type     | Null | Key | Default | Extra |
+-------------+----------+------+-----+---------+-------+
| id          | char(15) | NO   | PRI | NULL    |       |
| accession   | char(30) | NO   | MUL | NULL    |       |
| description | text     | YES  |     | NULL    |       |
+-------------+----------+------+-----+---------+-------+


*OMIM_std_all*
+-----------+---------------+------+-----+---------+-------+
| Field     | Type          | Null | Key | Default | Extra |
+-----------+---------------+------+-----+---------+-------+
| Disease   | varchar(1000) | YES  |     | NULL    |       |
| Gene_Name | varchar(100)  | YES  |     | NULL    |       |
+-----------+---------------+------+-----+---------+-------+


*ONTOLOGIES*
+-------------+---------+------+-----+---------+-------+
| Field       | Type    | Null | Key | Default | Extra |
+-------------+---------+------+-----+---------+-------+
| id          | char(2) | NO   | PRI | NULL    |       |
| description | text    | NO   |     | NULL    |       |
+-------------+---------+------+-----+---------+-------+


*PHENOTYPES*
+-------------+---------------+------+-----+---------+-------+
| Field       | Type          | Null | Key | Default | Extra |
+-------------+---------------+------+-----+---------+-------+
| id          | char(30)      | NO   |     | NULL    |       |
| type        | char(2)       | NO   |     | NULL    |       |
| accession   | char(30)      | NO   | MUL | NULL    |       |
| description | char(150)     | NO   |     | NULL    |       |
| score       | double        | NO   |     | NULL    |       |
| origin      | enum('0','1') | NO   |     | NULL    |       |
+-------------+---------------+------+-----+---------+-------+


*PROBABILITIES*
+-------------+----------+------+-----+---------+-------+
| Field       | Type     | Null | Key | Default | Extra |
+-------------+----------+------+-----+---------+-------+
| id          | char(15) | NO   | PRI | NULL    |       |
| position    | int(11)  | NO   | PRI | NULL    |       |
| A           | double   | NO   |     | NULL    |       |
| C           | double   | NO   |     | NULL    |       |
| D           | double   | NO   |     | NULL    |       |
| E           | double   | NO   |     | NULL    |       |
| F           | double   | NO   |     | NULL    |       |
| G           | double   | NO   |     | NULL    |       |
| H           | double   | NO   |     | NULL    |       |
| I           | double   | NO   |     | NULL    |       |
| K           | double   | NO   |     | NULL    |       |
| L           | double   | NO   |     | NULL    |       |
| M           | double   | NO   |     | NULL    |       |
| N           | double   | NO   |     | NULL    |       |
| P           | double   | NO   |     | NULL    |       |
| Q           | double   | NO   |     | NULL    |       |
| R           | double   | NO   |     | NULL    |       |
| S           | double   | NO   |     | NULL    |       |
| T           | double   | NO   |     | NULL    |       |
| V           | double   | NO   |     | NULL    |       |
| W           | double   | NO   |     | NULL    |       |
| Y           | double   | NO   |     | NULL    |       |
| information | double   | NO   |     | NULL    |       |
+-------------+----------+------+-----+---------+-------+


*PROTEIN*
+-------+-----------+------+-----+---------+-------+
| Field | Type      | Null | Key | Default | Extra |
+-------+-----------+------+-----+---------+-------+
| id    | int(11)   | NO   | MUL | NULL    |       |
| name  | char(100) | NO   | PRI | NULL    |       |
+-------+-----------+------+-----+---------+-------+


*PROVEAN_SIFT*
+---------------------+--------------+------+-----+---------+-------+
| Field               | Type         | Null | Key | Default | Extra |
+---------------------+--------------+------+-----+---------+-------+
| DBSNP_ID            | varchar(100) | YES  |     | NULL    |       |
| VARIANT             | varchar(500) | YES  |     | NULL    |       |
| PROTEIN_ID          | varchar(100) | YES  |     | NULL    |       |
| LENGTH              | varchar(100) | YES  |     | NULL    |       |
| STRAND              | varchar(100) | YES  |     | NULL    |       |
| CODON_CHANGE        | varchar(500) | YES  |     | NULL    |       |
| POS                 | varchar(100) | YES  |     | NULL    |       |
| RESIDUE_REF         | varchar(100) | YES  |     | NULL    |       |
| RESIDUE_ALT         | varchar(100) | YES  |     | NULL    |       |
| TYPE                | varchar(100) | YES  |     | NULL    |       |
| PROVEAN_SCORE       | varchar(100) | YES  |     | NULL    |       |
| PROVEAN_PREDICTION  | varchar(100) | YES  |     | NULL    |       |
| PROVEAN_NUM_SEQ     | varchar(100) | YES  |     | NULL    |       |
| PROVEAN_NUM_CLUSTER | varchar(100) | YES  |     | NULL    |       |
| SIFT_SCORE          | varchar(100) | YES  |     | NULL    |       |
| SIFT_PREDICTION     | varchar(100) | YES  |     | NULL    |       |
| SIFT_MEDIAN_INFO    | varchar(100) | YES  |     | NULL    |       |
| SIFT_NUM_SEQ        | varchar(100) | YES  |     | NULL    |       |
+---------------------+--------------+------+-----+---------+-------+


*REACTOME_EntrezID*
+-----------+----------------+------+-----+---------+-------+
| Field     | Type           | Null | Key | Default | Extra |
+-----------+----------------+------+-----+---------+-------+
| Pathway   | varchar(500)   | YES  |     | NULL    |       |
| Link2gsea | varchar(500)   | YES  |     | NULL    |       |
| EntrezID  | varchar(25000) | YES  |     | NULL    |       |
+-----------+----------------+------+-----+---------+-------+


*SEQUENCE*
+----------+---------+------+-----+---------+-------+
| Field    | Type    | Null | Key | Default | Extra |
+----------+---------+------+-----+---------+-------+
| id       | int(11) | NO   | PRI | NULL    |       |
| sequence | text    | NO   | MUL | NULL    |       |
+----------+---------+------+-----+---------+-------+


*TAGC128_hg19*
+-------+--------------+------+-----+---------+-------+
| Field | Type         | Null | Key | Default | Extra |
+-------+--------------+------+-----+---------+-------+
| CHROM | varchar(3)   | NO   | PRI | NULL    |       |
| POS   | varchar(15)  | NO   | PRI | NULL    |       |
| REF   | varchar(250) | NO   | PRI | NULL    |       |
| ALT   | varchar(700) | NO   | PRI | NULL    |       |
| AC    | int(11)      | YES  |     | NULL    |       |
| AN    | int(11)      | YES  |     | NULL    |       |
| AF    | float(9,6)   | YES  |     | NULL    |       |
+-------+--------------+------+-----+---------+-------+


*TAGC128_hg38*
+-------+--------------+------+-----+---------+-------+
| Field | Type         | Null | Key | Default | Extra |
+-------+--------------+------+-----+---------+-------+
| CHROM | varchar(3)   | NO   | PRI | NULL    |       |
| POS   | varchar(15)  | NO   | PRI | NULL    |       |
| REF   | varchar(250) | NO   | PRI | NULL    |       |
| ALT   | varchar(700) | NO   | PRI | NULL    |       |
| AC    | int(11)      | YES  |     | NULL    |       |
| AN    | int(11)      | YES  |     | NULL    |       |
| AF    | float(9,6)   | YES  |     | NULL    |       |
+-------+--------------+------+-----+---------+-------+


*UNIPROT_PFAM_DOMAINS*
+-------------+--------------+------+-----+---------+-------+
| Field       | Type         | Null | Key | Default | Extra |
+-------------+--------------+------+-----+---------+-------+
| UNIPROT     | varchar(20)  | NO   | PRI | NULL    |       |
| PFAM        | varchar(20)  | YES  |     | NULL    |       |
| START       | int(11)      | NO   | PRI | NULL    |       |
| END         | int(11)      | NO   | PRI | NULL    |       |
| DOMAIN_NAME | varchar(500) | YES  |     | NULL    |       |
+-------------+--------------+------+-----+---------+-------+


*VARIANTS*
+--------------+-----------+------+-----+---------+-------+
| Field        | Type      | Null | Key | Default | Extra |
+--------------+-----------+------+-----+---------+-------+
| id           | char(25)  | NO   | MUL | NULL    |       |
| protein      | char(100) | NO   |     | NULL    |       |
| substitution | char(10)  | NO   |     | NULL    |       |
+--------------+-----------+------+-----+---------+-------+


*WEIGHTS*
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------+-----+---------+-------+
| Field   | Type                                                                                                                                                                                                                                                     | Null | Key | Default | Extra |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------+-----+---------+-------+
| id      | char(15)                                                                                                                                                                                                                                                 | NO   | PRI | NULL    |       |
| type    | enum('BLOOD','BLOOD_COAGULATION','CANCER','DEVELOPMENTAL','DIGESTIVE','EAR_NOSE_THROAT','ENDOCRINE','EYE','GENITOURINARY','HEART','IMMUNE','INHERITED','METABOLIC','MUSCULOSKELETAL','NERVOUS_SYSTEM','PSYCHIATRIC','REPRODUCTIVE','RESPIRATORY','SKIN') | NO   | PRI | NULL    |       |
| disease | double                                                                                                                                                                                                                                                   | NO   |     | NULL    |       |
| other   | double                                                                                                                                                                                                                                                   | NO   |     | NULL    |       |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------+-----+---------+-------+


*dbNSFP*
+--------------------------+--------------+------+-----+---------+-------+
| Field                    | Type         | Null | Key | Default | Extra |
+--------------------------+--------------+------+-----+---------+-------+
| CHROM                    | varchar(5)   | YES  |     | NULL    |       |
| POS                      | varchar(50)  | YES  |     | NULL    |       |
| REF                      | varchar(100) | YES  |     | NULL    |       |
| ALT                      | varchar(100) | YES  |     | NULL    |       |
| AA_REF                   | varchar(100) | YES  |     | NULL    |       |
| AA_ALT                   | varchar(100) | YES  |     | NULL    |       |
| Ensembl_geneid           | varchar(200) | YES  |     | NULL    |       |
| Ensembl_transcriptid     | varchar(500) | YES  |     | NULL    |       |
| Ensembl_proteinid        | varchar(500) | YES  |     | NULL    |       |
| AA_POS                   | varchar(100) | YES  |     | NULL    |       |
| SIFT_score               | varchar(100) | YES  |     | NULL    |       |
| SIFT_pred                | varchar(100) | YES  |     | NULL    |       |
| Uniprot_acc_Polyphen2    | varchar(500) | YES  |     | NULL    |       |
| Polyphen2_HDIV_score     | varchar(200) | YES  |     | NULL    |       |
| Polyphen2_HDIV_pred      | varchar(100) | YES  |     | NULL    |       |
| Polyphen2_HVAR_score     | varchar(200) | YES  |     | NULL    |       |
| Polyphen2_HVAR_pred      | varchar(100) | YES  |     | NULL    |       |
| LRT_score                | varchar(200) | YES  |     | NULL    |       |
| LRT_pred                 | varchar(100) | YES  |     | NULL    |       |
| MutationTaster_score     | varchar(200) | YES  |     | NULL    |       |
| MutationTaster_pred      | varchar(100) | YES  |     | NULL    |       |
| MutationAssessor_score   | varchar(200) | YES  |     | NULL    |       |
| MutationAssessor_pred    | varchar(100) | YES  |     | NULL    |       |
| FATHMM_score             | varchar(200) | YES  |     | NULL    |       |
| FATHMM_pred              | varchar(100) | YES  |     | NULL    |       |
| PROVEAN_score            | varchar(200) | YES  |     | NULL    |       |
| PROVEAN_pred             | varchar(100) | YES  |     | NULL    |       |
| CADD_phred               | varchar(200) | YES  |     | NULL    |       |
| DANN_score               | varchar(200) | YES  |     | NULL    |       |
| fathmm_MKL_coding_score  | varchar(200) | YES  |     | NULL    |       |
| fathmm_MKL_coding_pred   | varchar(100) | YES  |     | NULL    |       |
| MetaSVM_score            | varchar(200) | YES  |     | NULL    |       |
| MetaSVM_pred             | varchar(100) | YES  |     | NULL    |       |
| MetaLR_score             | varchar(200) | YES  |     | NULL    |       |
| MetaLR_pred              | varchar(100) | YES  |     | NULL    |       |
| GERP_RS                  | varchar(200) | YES  |     | NULL    |       |
| phyloP7way_vertebrate    | varchar(200) | YES  |     | NULL    |       |
| phyloP20way_mammalian    | varchar(200) | YES  |     | NULL    |       |
| phastCons7way_vertebrate | varchar(200) | YES  |     | NULL    |       |
| phastCons20way_mammalian | varchar(200) | YES  |     | NULL    |       |
| SiPhy_29way_logOadds     | varchar(200) | YES  |     | NULL    |       |
+--------------------------+--------------+------+-----+---------+-------+


*enst_ensg*
+-----------------------+-------------+------+-----+---------+-------+
| Field                 | Type        | Null | Key | Default | Extra |
+-----------------------+-------------+------+-----+---------+-------+
| ENSEMBL_TRANSCRIPT_ID | varchar(20) | YES  |     | NULL    |       |
| ENSEMBL_GENE_ID       | varchar(20) | YES  |     | NULL    |       |
+-----------------------+-------------+------+-----+---------+-------+


*gID_Uniprot*
+---------+-------------+------+-----+---------+-------+
| Field   | Type        | Null | Key | Default | Extra |
+---------+-------------+------+-----+---------+-------+
| Gene_ID | varchar(20) | YES  |     | NULL    |       |
| Uniprot | varchar(10) | YES  |     | NULL    |       |
+---------+-------------+------+-----+---------+-------+


*gID_gNames*
+-----------+-------------+------+-----+---------+-------+
| Field     | Type        | Null | Key | Default | Extra |
+-----------+-------------+------+-----+---------+-------+
| Gene_ID   | varchar(20) | YES  |     | NULL    |       |
| Gene_Name | varchar(40) | YES  |     | NULL    |       |
+-----------+-------------+------+-----+---------+-------+


*gID_gSymbol*
+--------+-------------+------+-----+---------+-------+
| Field  | Type        | Null | Key | Default | Extra |
+--------+-------------+------+-----+---------+-------+
| ID     | varchar(40) | YES  |     | NULL    |       |
| Symbol | varchar(40) | YES  |     | NULL    |       |
+--------+-------------+------+-----+---------+-------+


*gSymbol_gID_biomart*
+--------+-------------+------+-----+---------+-------+
| Field  | Type        | Null | Key | Default | Extra |
+--------+-------------+------+-----+---------+-------+
| Symbol | varchar(40) | YES  |     | NULL    |       |
| ID     | varchar(40) | YES  |     | NULL    |       |
+--------+-------------+------+-----+---------+-------+


*gerp_elements*
+---------+-------------+------+-----+---------+-------+
| Field   | Type        | Null | Key | Default | Extra |
+---------+-------------+------+-----+---------+-------+
| CHROM   | varchar(5)  | YES  |     | NULL    |       |
| START   | varchar(20) | YES  |     | NULL    |       |
| END     | varchar(20) | YES  |     | NULL    |       |
| LENGTH  | varchar(20) | YES  |     | NULL    |       |
| SCORE   | varchar(20) | YES  |     | NULL    |       |
| P_VALUE | varchar(20) | YES  |     | NULL    |       |
+---------+-------------+------+-----+---------+-------+


*hg19_hg38_mapping*
+-------+-------------+------+-----+---------+-------+
| Field | Type        | Null | Key | Default | Extra |
+-------+-------------+------+-----+---------+-------+
| CHROM | varchar(3)  | NO   | PRI | NULL    |       |
| hg_19 | varchar(30) | NO   | PRI | NULL    |       |
| hg_38 | varchar(30) | YES  |     | NULL    |       |
+-------+-------------+------+-----+---------+-------+


*hg38_hg19_mapping*
+-------+-------------+------+-----+---------+-------+
| Field | Type        | Null | Key | Default | Extra |
+-------+-------------+------+-----+---------+-------+
| CHROM | varchar(3)  | NO   | PRI | NULL    |       |
| hg_38 | varchar(30) | NO   | PRI | NULL    |       |
| hg_19 | varchar(30) | YES  |     | NULL    |       |
+-------+-------------+------+-----+---------+-------+


*m3D_modbase_models*
+-----------------------+-------------+------+-----+---------+-------+
| Field                 | Type        | Null | Key | Default | Extra |
+-----------------------+-------------+------+-----+---------+-------+
| uniprot               | varchar(10) | YES  |     | NULL    |       |
| template_length       | int(11)     | YES  |     | NULL    |       |
| target_length         | int(11)     | YES  |     | NULL    |       |
| template_pdb          | varchar(5)  | YES  |     | NULL    |       |
| target_begin          | int(11)     | YES  |     | NULL    |       |
| target_end            | int(11)     | YES  |     | NULL    |       |
| sequence_identity     | float(5,2)  | YES  |     | NULL    |       |
| model_score           | float(3,2)  | YES  |     | NULL    |       |
| modpipe_quality_score | float(8,6)  | YES  |     | NULL    |       |
| zDOPE                 | float(4,2)  | YES  |     | NULL    |       |
| eVALUE                | float(4,2)  | YES  |     | NULL    |       |
| modbase_modelID       | varchar(40) | YES  |     | NULL    |       |
+-----------------------+-------------+------+-----+---------+-------+


*m3D_pdb_structures_index*
+-----------------+-------------+------+-----+---------+-------+
| Field           | Type        | Null | Key | Default | Extra |
+-----------------+-------------+------+-----+---------+-------+
| uniprot         | varchar(10) | YES  |     | NULL    |       |
| template_length | int(11)     | YES  |     | NULL    |       |
| uniprot_length  | int(11)     | YES  |     | NULL    |       |
| template_pdb    | varchar(5)  | YES  |     | NULL    |       |
| uniprot_begin   | int(11)     | YES  |     | NULL    |       |
| uniprot_end     | int(11)     | YES  |     | NULL    |       |
| pdb_begin       | int(11)     | YES  |     | NULL    |       |
| pdb_end         | int(11)     | YES  |     | NULL    |       |
| modelID         | varchar(30) | YES  |     | NULL    |       |
+-----------------+-------------+------+-----+---------+-------+


*pdbid_chain_uniprot*
+---------+-------------+------+-----+---------+-------+
| Field   | Type        | Null | Key | Default | Extra |
+---------+-------------+------+-----+---------+-------+
| pdbid   | varchar(5)  | YES  |     | NULL    |       |
| chain   | varchar(1)  | YES  |     | NULL    |       |
| uniprot | varchar(10) | YES  |     | NULL    |       |
+---------+-------------+------+-----+---------+-------+

