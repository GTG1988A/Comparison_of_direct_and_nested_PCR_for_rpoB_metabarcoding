# Comparison of direct and nested PCR for rpoB metabarcoding
![image](https://github.com/GTG1988A/Comparison_of_direct_and_nested_PCR_for_rpoB_metabarcoding/assets/83345122/7e23c06f-a8a6-4d30-a595-06a7926a00b9)

## 1. Download genomes with genome_update
### Genome_update installation
Documentation:
https://github.com/pirovc/genome_updater
https://ftp.ncbi.nlm.nih.gov/genomes/all/README.txt
https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt


```bash=
wget --quiet --show-progress https://raw.githubusercontent.com/pirovc/genome_updater/master/genome_updater.sh
chmod +x genome_updater.sh
```
### Download 

I download genomes at all assembly levels: complete genome, chromosome, scaffold and contig. So I get all the rpob sequences I can.
FecthMGS needs the protein dna of the genes. If you also want to have the nucleic sequences, you have to give it that too. So I download this:

```bash=
# protein sequence file FAA
./genome_updater.sh -o "arc_refseq_cg" -d "refseq" -g "bacteria" -f "protein.faa.gz" -t 12

# protein nucleic acid sequence file FNA
./genome_updater.sh -o "arc_refseq_cg" -d "refseq" -g "bacteria" -f "cds_from_genomic.fna.gz" -t 12

```

## 2. Recovering rpoBs with fetchMGs

I need the .faa protein files and the .fna sequence files so that fetchMGS can give me the nucleic acid version.

But for this to work, the two files must have the same sequence identifiers 


### 2.a Preparing .faa files 
the files are in the form:
```bash=
$ head GCF_000219415.3_ASM21941v3_protein.faa
>WP_002964358.1 MULTISPECIES: 30S ribosomal protein S19 [Hyphomicrobiales]
MARSVWKGPFVDGYLLKKAEKVREGGRNEVIKMWSRRSTILPQFVGLTFGVYNGNKHVPVSVSEEMVGHKFGEFAPTRTYYGHGADKKAKRK
```
so I made a bash script to keep only the protein identifier (WP_ ...)

/!\please don't forget to change your variable paths in the script 

```bash=
./rename_faa.sh
```

### 2.b Preparing .fna files

Files are in format:

```bash=
$ head GCF_000219415.3_ASM21941v3_cds_from_genomic.fna
>lcl|NZ_CP137016.1_cds_WP_244428070.1_1 [locus_tag=SFGR64A_RS00035] [protein=hypothetical protein] [protein_id=WP_244428070.1] [location=complement(6240..6629)] [gbkey=CDS]
GTGGTTGCCCTTCATGGCGCCGTTCCGCGCCGTTCCGCCGCGTCCTTCACCAATATCAAGGTTCGCGTGCGCGACGACCG(..)ACAGGAAGCGCTCGTCGCCAACATCGCCGTGGCGCTCAGAAAGCTCTCGATAGGCGAATTGTTGCTTTAG
```

I only keep the protein_id section with the script:

/!\please don't forget to change your variable paths in the script 
```bash=
./rename_fna.sh

```

NB: This script does not save sequences if they do not have a protein identifier because if there is no protein equivalent, FetchMGS does not run on nucleic acids. They are therefore useless. 

### 2.c  Launch fetchMGS

You need to run fetch MGS for about 80,000 files. So if I give the command line as it is : 

```bash=
./fetchMGs/fetchMGs.pl -m extraction test/arc_refseq_cg/2024-01-08_11-31-26/files/*.faa -d arc_refseq_cg/2024-01-08_07-45-09/files/*.fna -c COG0085
```
I'm going to get an error because there are too many entries. 

So, I get the file names: 

```
ls protein_seq/arc_refseq_cg/2024-02-05_09-01-16/files/modified/ > tmp.txt
sed 's/_protein_modified\.faa//' > names_files.txt
rm tmp.txt
```
and I create a script that loops through the command line by file 

/!\please don't forget to change your variable paths in the script 
```bash=
./launch_FetchMGS.sh
```

It takes a while to run, but it's fetchMGS and you only have to do it once.


Fetch creates a folder for each genome, containing the files COG0085.faa and COG0085.fna. COG0085 corresponds to rpoB, which is the one we are interested in. 
## 3. Formatting the file

### 3.a Add Taxid
The taxid must be added to launch obiconvert before EcoPCR. I make a table of correspondence between the genome identifier and its species taxids, which are columns 1 and 7 of the assembly summary downloaded by genome_updater:

```
cut -f1,7 genome/arc_refseq_cg/2024-02-05_09-00-49/assembly_summary.txt > correspondance_tab.txt
```

Thanks to this correspondence table, I've made a python script which will go into each folder resulting from FetchMGS; find the ID and the corresponding taxid then it will add it to the name of the sequence next to the protein ID.

/!\please don't forget to change your variable paths in the script 
```python=
python add_taxid.py
```


### 3.b Concatenation of fna files to obtain our rpob regions 

They must then be concatenated to obtain the complete file with the rpob genes in nucleic acid.

```bash=
find output/ -name "COG0085.fna" -exec cat {} + > concatenated_COG0085.fna
```

## 4. Formatting databases for EcoPCR: Obitools
### Obitools for rpod sequences

```bach!
mkdir ecoPCR_db/
obiconvert --fasta concatenated_COG0085.fna --ecopcrdb-output=ecoPCR_db/rpob -t ncbi_tax_dumb/
```

NB: it is possible to have an error because sometimes genome updater has more IDs than the ncbi taxonomy which has not yet been updated. I didn't get this error here. 

## 5. Launch EcoPCR

### 1st PCR for Nested 

| rpoB_F  | rpoB_R |
| -------- | -------- | 
| CAGYTDTCNCARTTYATGGAYCA|AGTTRTARCCDTYCCANGKCAT| 


```bash=
$ ecoPCR -d ecoPCR_db/rpob -e 2 CAGYTDTCNCARTTYATGGAYCA AGTTRTARCCDTYCCANGKCAT  > primers/rpob_l70_L8000_e2.ecopcr
```

I retrieve the names of the genomes that have been amplified, those present in the EcoPCR result. I retrieve their species taxid and then their taxonomy.
The genomes that are not in EcoPCR are also recovered, i.e. those that have not been amplified. I also recover their species ID.
If they are in the EcoPCR result, they have a "True" in the "amplified" column, and if they are not, they have a "False".
I then run the script for the Krona

```bash=
#Recovers amplified protein IDs
awk '!/^#/{print $1}' FS="|" ecoPCR_rpob1/rpob_l70_L8000_e2.ecopcr > id_amplified.txt

# Total IDs in all genomes
grep ">" concatenated_COG0085.fna | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt

#To find out the difference between the two
grep -v -F -f id_amplified.txt id_seq_total.txt > not_amplified.txt
```

I've got the sequence identifiers now I want the taxids:

```bash=
#my file with sequence names only
grep ">" concatenated_COG0085.fna > name_seq.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' not_amplified.txt name_seq.txt > tmp_no_amplified.txt
cut -d'|' -f2 tmp_no_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > no_amplified_speciesTAXID.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' id_amplified.txt name_seq.txt > tmp_amplified.txt
cut -d'|' -f2 tmp_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > amplified_speciesTAXID.txt

```
I can use the TaxonKit tool to access their taxonomies. 

To do this, download the folder containing the ncbi taxonomy from this link: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

```bash=
cd ncbi_tax_dumb/2024-15-02
wget -r https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
gunzip new_taxdump.tar.gz
tar -xvf new_taxdump.tar.gz
```

```bash=
# not amplified
$ cat no_amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > no_amplified_species.taxo
# amplified
$ cat amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > amplified_species.taxo
```

I remove all possible spaces in the names:

```bash=
$ sed 's/ /_/g' amplified_species.taxo > amplified_species_wc.taxo
$ sed 's/ /_/g' no_amplified_species.taxo > no_amplified_species_wc.taxo
```
I reformat these files so that they are in the format accepted by krona:
* I add a "Sequence" column which corresponds to a sequence sum. They are all 1
* I add an "amplified" column which is True or False. So True for the amplified file and False for the non-amplified file.

NB: I'm changing the True and False by hand in the script because I didn't add the argument.

/!\please don't forget to change your variable paths in the script 
```python=
python format_krona.py input_file output_file
```


Then I create a header:
```bash=
$ cat header.tsv
seqID	sequence	taxonomy	amplified
```
I concatenate everything:

```bash=
$ cat header.tsv amplified_species_wc_reformat.tsv > tmp.tsv
$ cat tmp.tsv no_amplified_species_wc_reformat.tsv > final_results.tsv
```

I'm removing the column corresponding to the taxid because it's not accepted by the krona script:

```bash=
cut -f 2-4 final_results.tsv > input.tsv
```

I run the script to get the krona:
TODO: add  make_krona_xml_bis.py
```bash=
$ python make_krona_xml_bis.py input.tsv output.xml  --name rpob -v --color_label amplified  --dataset_labels sequence
# launch krona
$ ktImportXML output.xml

```

### 2nde PCR 
| Univ_rpoB_F_deg  | Univ_rpoB_R_deg |
| -------- | -------- | 
| GGYTWYGAAGTNCGHGACGTDCA|TGACGYTGCATGTTBGMRCCCATMA| 


```bash=
$ ecoPCR -d ecoPCR_db/rpob -e 2 GGYTWYGAAGTNCGHGACGTDCA TGACGYTGCATGTTBGMRCCCATMA  > primers/rpob_l70_L8000_e2.ecopcr

```

=> same treatment as above

### Nested PCR

For this one, I took the EcoPCR results for the first pair of primers called rpoB_F and rpoB_R.
Then I kept only the genomes that appeared in the EcoPCR results, so when the primers worked.

```bash!
#Recovers amplified protein IDs
awk '!/^#/{print $1}' FS="|" ecoPCR_rpob1/rpob_l70_L8000_e2.ecopcr > id_amplified.txt

# Total IDs in all genomes
grep ">" concatenated_COG0085.fna | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt

# To find out the difference between the two
grep -v -F -f liste_id_ecoPCR_sort.txt liste_id_concatened_COG0085_sort.txt > not_amplified.txt
```

I use this script to remove non-amplified genomes 

/!\please don't forget to change your variable paths in the script 
```bash=
python deletion_seq_in_fasta.py
```
Now I'm going to reformat the fasta with Obitools to run EcoPCR on it.

```bash=
mkdir ecoPCR_db/
$ obiconvert --fasta fichier_filtre.fasta --ecopcrdb-output=ecoPCR_db/rpob -t ../ncbi_tax_dumb/2024-15-02/
```

ecoPCR launch with 2nd primers
| Univ_rpoB_F_deg  | Univ_rpoB_R_deg |
| -------- | -------- | 
| GGYTWYGAAGTNCGHGACGTDCA|TGACGYTGCATGTTBGMRCCCATMA| 

```bash=
$ ecoPCR -d ecoPCR_db/rpob -e 2  GGYTWYGAAGTNCGHGACGTDCA TGACGYTGCATGTTBGMRCCCATMA  > primers/rpob_l70_L8000_e2.ecopcr

```
I recover the identifiers of the species amplified by the primers

```bash=
$ awk '!/^#/{print $1}' FS="|" primers/rpob_l70_L8000_e2.ecopcr > id_seq_amplified_species.txt

```

I also retrieve the same list of identifiers from my input fasta file in order to have a list of all those that are not in common. So that's all the species that haven't been amplified.
```bash=
grep ">" fichier_filtre.fasta | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt
```

to obtain the list of non-amplified:
```bash=
grep -v -F -f id_seq_amplified_species.txt id_seq_total.txt > id_no_amplified_species.txt
```

I recover the taxids of those that have been amplified as well as those that have not been amplified:
```bash=
#my file with sequence names only
grep ">" concatenated_COG0085.fna > name_seq.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' id_no_amplified_species.txt name_seq.txt > tmp_no_amplified.txt
cut -d'|' -f2 tmp_no_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > no_amplified_speciesTAXID.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' id_seq_amplified_species.txt name_seq.txt > tmp_amplified.txt
cut -d'|' -f2 tmp_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > amplified_speciesTAXID.txt
```
Taxonkit:

```bash=
# not amplified
$ cat no_amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > no_amplified_species.taxo
# amplified
$ cat amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > amplified_species.taxo
```
I remove all possible spaces in the names:

```bash=
$ sed 's/ /_/g' amplified_species.taxo > amplified_species_wc.taxo
$ sed 's/ /_/g' no_amplified_species.taxo > no_amplified_species_wc.taxo
```

I reformat these files so that they are in the format accepted by krona:
* I add a "Sequence" column which corresponds to a sequence sum. They are all 1
* I add an "amplified" column which is True or False. So True for the amplified file and False for the non-amplified file.
  
NB: I change the True and False by hand in the script because I didn't add the argument.

header
```bash=
$ cat header.tsv
seqID	sequence	taxonomy	amplified
```

concatenate:

```bash=
$ cat header.tsv amplified_species_wc_reformat.tsv > tmp.tsv
$ cat tmp.tsv no_amplified_species_wc_reformat.tsv > final_results.tsv
```

I'm removing the column corresponding to the taxid because it's not accepted by the krona script:

```bash=
cut -f 2-4 final_results.tsv > input.tsv
```

script for krona
```bash=
$ python make_krona_xml_bis.py input.tsv output.xml  --name rpob -v --color_label amplified  --dataset_labels sequence
# launch krona
$ ktImportXML output.xml

```

# Comparison with  16S
I have a database of 16S sequences extract from refseq genomes which use gff, which I can download using genome_updater.

```bash=
./../genome_updater.sh -o "arc_refseq_cg" -d "refseq" -g "bacteria" -f "genomic.gff.gz,genomic.fna.gz" -t 12 -N
```
for this we need to uses the script src/extract_16S23S.py

we need 2 files genome_dirs.txt and genome_files_with_taxid.txt created respectively by the scripts find_dirs.sh and find_name.py by adapting the paths of the download databases. 
the two files resemble a: 
```bash=
gagoutin@genobioinfo2 ~/work/test_16S/base_16S_Geraldine $ head -2 genome_dirs.txt
/home/gagoutin/work/base_données/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415
/home/gagoutin/work/base_données/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395
$ head -2 genome_files_with_taxid.txt
/home/gagoutin/work/base_données//NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415/GCA_018630415.1_ASM1863041v1_genomic.fna.gz GCA_018630415.1 1280
/home/gagoutin/work/base_données//NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395/GCA_018630395.1_ASM1863039v1_genomic.fna.gz GCA_018630395.1 1280
```
and then run the script:

```bash=
python extract_16S23S.py -a genome_dirs.txt -t genome_files_with_taxid.txt
```

The genomes in this database are kept only if they also have an rpoB sequence to make the comparison reliable.
I retrieve the folder names of the fetchMGS results and I retrieve the names of the genomes that have a non-empty fna file.

```bash=
l output_fetchMGS/*/*.fna > fna_files.txt
python find_fna_not_empty.py > genomes_names_rpoB.txt
```

For the 16S genomes, after running the extract script  extract_16S23S.py  we get a tsv file that looks like this:
```bash=
$ head 16S.tsv
gene_id seqid   product gene_name       start   end     strand  circular        seqid_length    partial genome_name     species_taxid
rna-MAG715_00878        JAANXF010000016.1       16S ribosomal RNA       16S     22324   23807   +       False   83241   False   GCA_018263395.1 2687197
rna-MAG471_00900        JAANXB010000052.1       16S ribosomal RNA       16S     6007    7522    -       False   42085   False   GCA_018263555.1 2024894
rna-MAG458_00410        JAANXA010000013.1       16S ribosomal RNA       16S     19451   20920   +       False   59501   False   GCA_018263505.1 2024843
rna-MAG453_01709        JAANWZ010000126.1       16S ribosomal RNA       16S     16263   17787   +       False   24252   False   GCA_018263585.1 2212469
rna-MAG551_02765        JAANXD010000102.1       16S ribosomal RNA       16S     242     1810    +       False   17091   False   GCA_018263435.1 1127984
rna-MAG451_03223        JAANWY010000396.1       16S ribosomal RNA       16S     5367    6870    -       False   14192   False   GCA_018263655.1 2073117
rna-MAG431_00147        JAANWX010000004.1       16S ribosomal RNA       16S     75386   76912   +       False   107298  False   GCA_018263575.1 2026724
rna-MAG581_01280        JAANXE010000022.1       16S ribosomal RNA       16S     141497  143009  -       False   162169  False   GCA_018263475.1 2026735
rna-KGY51_12075 JAGXLS010000180.1       16S ribosomal RNA       16S     1       273     -       False   2851    True    GCA_018334605.1 1921580
```

So we get the 11th column which is the gene_name, I take it out -u and I get this list back

```bash=
zcat 16S.tsv.gz| cut -f11 | sort -u > /PATH/genome_name_16S.txt
#...
cut -d'_' -f1,2 names_files.txt > genome_names_rpob.txt
```

I recover those that are in common, so I look for the rpob genomes in the 16S genomes and if it's found it's good:
grep -F -f genome_names_rpob.txt  genome_name_16S.txt > common_genomes.txt

We therefore need to remove those that are not in common in the 16S fasta.
I'm looking for which 16S genome is not in the genomes that have rpob

```bash
grep -v -F -f genomes_names_rpoB.txt genome_name_16S.txt > not_common.txt
```

I retrieve the rows corresponding to these genes from my tsv:
```bash=
awk -F'\t'  'NR==FNR{genes[$1]; next} $11 in genes' not_common.txt 16S.tsv > output.tsv
gene_id seqid   product gene_name       start   end     strand  circular        seqid_length    partial genome_name     species_taxid
rna-HQ400_RS13095       NZ_CP053879.1   16S ribosomal RNA       16S     2835345 2836889 -       True    4535455 False   GCF_018802265.1 650
rna-HQ400_RS14170       NZ_CP053879.1   16S ribosomal RNA       16S     3049597 3051141 -       True    4535455 False   GCF_018802265.1 650
rna-HQ400_RS14960       NZ_CP053879.1   16S ribosomal RNA       16S     3225998 3227542 -       True    4535455 False   GCF_018802265.1 650
rna-HQ400_RS16380       NZ_CP053879.1   16S ribosomal RNA       16S     3513214 3514758 -       True    4535455 False   GCF_018802265.1 650
```

This will allow me to retrieve the seqid that will allow me to delete this sequence in my fasta using the delete_sequence_in_a_fasta.py script.

```bash=
cut -f2 output.tsv > seqid_not_com.txt
python delete_sequence_in_a_fasta.py
```

Formatting my new fasta so that it is compatible with EcoPCR with obiconvert 

```bash=
mkdir ecoPCR_db/
obiconvert --fasta new_fasta_file.fasta --ecopcrdb-output=ecoPCR_db/16S -t ncbi_tax_dumb/2024-15-04/
```
if you ever get obiconvert errors due to an unknown ID, delete the sequences concerned. They are due to a lack of coherence between the different taxonomies.

Now going to use this fasta to launch my EcoPCR 16S.


| 16S_F  | 16S_R |
| -------- | -------- | 
| GAAGAGTTTGATCATGGCTC|AAGGAGGTGATCCAGCCGCA|

```bash=
ecoPCR -d ecoPCR_db/16S -e 2  GAAGAGTTTGATCATGGCTC AAGGAGGTGATCCAGCCGCA  > 16S.ecopcr
```

We retrieve the identifiers of the species amplified by the primers
```bash=
awk '!/^#/{print $1}' FS="|" 16S.ecopcr > id_amplified.txt
```

We also retrieve the same list of identifiers from my input fasta file in order to have a list of all those that are not in common. So that's all the species that haven't been amplified.
```bash=
grep ">" fasta_file.fasta | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt
```
To find out the difference between the two
```bash=
grep -v -F -f id_amplified.txt id_seq_total.txt > not_amplified.txt
```

Then turn the steps as above to make the krona graph
