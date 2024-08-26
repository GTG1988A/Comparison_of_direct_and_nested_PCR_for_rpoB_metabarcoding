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

I download genomes at assembly levels complete genome and chromosomeL. I don't take scaffold and contig, to get the best possible quality rpoB sequences.
FecthMGS needs the protein dna of the genes. If you also want to have the nucleic sequences, you have to give it that too. So I download this:

```bash=
# protein sequence file FAA
./genome_updater.sh -o "arc_refseq_cg" -d "refseq" -g "bacteria" -f "protein.faa.gz" -l "complete genome,chromosome" -t 12

# protein nucleic acid sequence file FNA
./genome_updater.sh -o "arc_refseq_cg" -d "refseq" -g "bacteria" -f "cds_from_genomic.fna.gz" -l "complete genome,chromosome" -t 12

```

## 2. Recovering rpoBs with fetchMGs

We require the .faa protein files and the .fna sequence files to allow fetchMGS to generate the nucleic acid version.
To ensure this functionality, the .faa and .fna files must have matching sequence identifiers.

### 2.a Preparing .faa files 
the files are in the form:
```bash=
$ head GCF_000219415.3_ASM21941v3_protein.faa
>WP_002964358.1 MULTISPECIES: 30S ribosomal protein S19 [Hyphomicrobiales]
MARSVWKGPFVDGYLLKKAEKVREGGRNEVIKMWSRRSTILPQFVGLTFGVYNGNKHVPVSVSEEMVGHKFGEFAPTRTYYGHGADKKAKRK
```
You need to use the script rename_faa.sh which only retains the protein identifier 

```bash=
./script/rename_faa.sh protein_seq/ protein_seq/modified_faa
```

### 2.b Preparing .fna files

Files are in format:

```bash=
$ head GCF_000219415.3_ASM21941v3_cds_from_genomic.fna
>lcl|NZ_CP137016.1_cds_WP_244428070.1_1 [locus_tag=SFGR64A_RS00035] [protein=hypothetical protein] [protein_id=WP_244428070.1] [location=complement(6240..6629)] [gbkey=CDS]
GTGGTTGCCCTTCATGGCGCCGTTCCGCGCCGTTCCGCCGCGTCCTTCACCAATATCAAGGTTCGCGTGCGCGACGACCG(..)ACAGGAAGCGCTCGTCGCCAACATCGCCGTGGCGCTCAGAAAGCTCTCGATAGGCGAATTGTTGCTTTAG
```

We only keep the protein_id section with the script rename_fna.sh 

/!\please don't forget to change your variable paths in the script 
```bash=
./script/rename_fna.sh cds_genomic_seq/ cds_genomic_seq/modified_fna
```
NB: This script does not save sequences if they do not have a protein identifier because if there is no protein equivalent, FetchMGS does not run on nucleic acids. They are therefore useless. 

### 2.c  Launch fetchMGS
You need a file with all the file names in the protein folder. You can use the command:

```
protein_seq/modified_faa/ | sed 's/_protein_modified\.faa//' > names_files.txt
```

and then run the script_fetch.sh.

NB: this can be a lengthy process.

The perl version of FetchMGS is used here. There is now a python version, not used here. (https://github.com/motu-tool/fetchMGs)

```bash=
./script/script_fetch.sh names_files.txt output_fetchMGS cds_genomic_seq/modified_fna/ protein_seq/modified_faa/
```
Fetch creates a folder for each genome, containing the files COG0085.faa and COG0085.fna. COG0085 corresponds to rpoB, which is the one we are interested in.

## 3. Formatting the file

### 3.a Add Taxid
The taxid must be added to launch obiconvert before EcoPCR. You need to make  a table of correspondence between the genome identifier and its species taxids, which are columns 1 and 7 of the assembly summary downloaded by genome_updater:

```
cut -f1,7 protein_seq/assembly_summary_light.txt > correspondance_tab.txt
```

We now need to use a python script which will go into each FetchMGS result folder and add the taxid corresponding to the genome in the output files. This taxid is added to the sequence name in the form "| taxid=XXX;", which is very important for what follows.

```python=
module load devel/python/Python-3.11.1
python script/add_taxid.py correspondance_tab.txt output_fetchMGS/
```

### 3.b Concatenation of fna files to obtain our rpob regions 

They must then be concatenated to obtain the complete file with the rpob genes in nucleic acid.

```bash=
find output_fetchMGS/ -name "COG0085.fna" -exec cat {} + > concatenated_COG0085.fna
```

## 4. Formatting databases for EcoPCR: Obitools
### Obitools for rpod sequences
To do this, download the folder containing the ncbi taxonomy from this link: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

```bash=
cd ncbi_tax_dumb/2024-15-02
wget -r https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
gunzip new_taxdump.tar.gz
tar -xvf new_taxdump.tar.gz
```


Run a script that adds a number tag to indicate the sequence number so that there is no ambiguity if there are several rpoB or 16S sequences in the same genome. I also cut the identifier so that it is identical with the EcoPCR output because ecoPCR is limited to 19 characters. 

```bash
python format_fasta_for_ecopcr.py concatenated_COG0085.fna concatenated_COG0085_format.fna
```
then, turn obiconvert:

```bach!
mkdir ecoPCR_db/
module load devel/Miniconda/Miniconda3
module load bioinfo/OBITools/1.2.11
obiconvert --fasta concatenated_COG0085_format.fna --ecopcrdb-output=ecoPCR_db/rpob -t ncbi_tax_dumb/
```
Warning: errors may occur if a taxid is out of date in the ncbi taxonomy or if it has been deleted. Delete the problematic sequences.

## 5. Launch EcoPCR

### 1st PCR for Nested 

```bash
mkdir nested_1st_PCR
cd nested_1st_PCR
```

| rpoB_F  | rpoB_R |
| -------- | -------- | 
| CAGYTDTCNCARTTYATGGAYCA|AGTTRTARCCDTYCCANGKCAT| 


```bash=
ecoPCR -d ecoPCR_db/rpob -e 2 CAGYTDTCNCARTTYATGGAYCA AGTTRTARCCDTYCCANGKCAT  > nested_1st_PCRrpob.ecopcr
```

You must retrieve the names of the genomes that have been amplified, i.e. those present in the EcoPCR result. 
Genomes that are not in the EcoPCR are also retrieved, i.e. those that have not been amplified.
Their taxids are also recovered for taxonomy.
 
```bash=
awk -F"|" '!/^#/{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $4); print ">" $1 "| taxid=" $4 ";"}' nested_1st_PCRrpob.ecopcr > id_amplified.txt

# Total seq name in all genomes
grep ">" ../concatenated_COG0085_format.fna  > name_seq.txt

#To find out the difference between the two
grep -v -F -f id_amplified.txt name_seq.txt > not_amplified.txt
```

For the taxid:

```bash=
cut -d'|' -f2 not_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > no_amplified_speciesTAXID.txt
cut -d'|' -f2 id_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > amplified_speciesTAXID.txt

```
I can use the TaxonKit tool to access their taxonomies. 
```bash=
module load bioinfo/TaxonKit/0.15.0
# not amplified
cat no_amplified_speciesTAXID.txt | taxonkit lineage --data-dir ncbi_tax_dumb/ | taxonkit reformat --data-dir ncbi_tax_dumb/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > no_amplified_species.taxo
# amplified
cat amplified_speciesTAXID.txt | taxonkit lineage --data-dir ncbi_tax_dumb/ | taxonkit reformat --data-dir ncbi_tax_dumb/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > amplified_species.taxo
```

Remove all possible spaces in the names:

```bash=
sed 's/ /_/g' amplified_species.taxo > amplified_species_wc.taxo
sed 's/ /_/g' no_amplified_species.taxo > no_amplified_species_wc.taxo
```

Reformat these files so that they are in the format accepted by krona:
* Add a "Sequence" column which corresponds to a sequence sum. They are all 1
* Aadd an "amplified" column which is True or False. So True for the amplified file and False for the non-amplified file.

```python=
python ../script/format_krona.py amplified_species_wc.taxo amplified_species_wc_reformat.tsv True
python ../script/format_krona.py no_amplified_species_wc.taxo no_amplified_species_wc_reformat.tsv False
```


Then I create a header:
```bash=
$ cat header.tsv
seqID	sequence	taxonomy	amplified
```
I concatenate everything:

```bash=
cat header.tsv amplified_species_wc_reformat.tsv > tmp.tsv
cat tmp.tsv no_amplified_species_wc_reformat.tsv > final_results.tsv
```

I'm removing the column corresponding to the taxid because it's not accepted by the krona script:

```bash=
cut -f 2-4 final_results.tsv > input.tsv
```

I run the script to get the krona:


```bash=
$ python make_krona_xml_bis.py input.tsv output.xml  --name rpob -v --color_label amplified  --dataset_labels sequence
# launch krona
$ ktImportXML output.xml

```

### direct_stategy

```bash
mkdir direct_stategy
cd direct_stategy
```
| Univ_rpoB_F_deg  | Univ_rpoB_R_deg |
| -------- | -------- | 
| GGYTWYGAAGTNCGHGACGTDCA|TGACGYTGCATGTTBGMRCCCATMA| 


```bash=
ecoPCR -d ecoPCR_db/rpob -e 2 GGYTWYGAAGTNCGHGACGTDCA TGACGYTGCATGTTBGMRCCCATMA  > direct_stategy.ecopcr
```

=> same treatment as above

### Nested_1st_PCR_then_2nd_PCR
```bash
mkdir nested_1st_PCR_then_2nd_PCR
cd nested_1st_PCR_then_2nd_PCR
```
For this one, We need to took the EcoPCR results for the first pair of primers called rpoB_F and rpoB_R.
Then I kept only the genomes that appeared in the EcoPCR results, so when the primers worked.

Use the deletion_seq_in_fasta.py script to remove non-amplified genomes 

/!\please don't forget to change your variable paths in the script 
```bash=
python deletion_seq_in_fasta.py ../concatenated_COG0085.fna nested_1st_PCR/id_amplified.txt new_fasta.fasta
```
No we need to reformat the fasta with Obitools to run EcoPCR on it.

```bash=
mkdir ecoPCR_db/
obiconvert --fasta new_fasta.fasta --ecopcrdb-output=ecoPCR_db/rpob -t ncbi_tax_dumb/
```

ecoPCR launch with 2nd primers
| Univ_rpoB_F_deg  | Univ_rpoB_R_deg |
| -------- | -------- | 
| GGYTWYGAAGTNCGHGACGTDCA|TGACGYTGCATGTTBGMRCCCATMA| 

```bash=
ecoPCR -d ecoPCR_db/rpob -e 2  GGYTWYGAAGTNCGHGACGTDCA TGACGYTGCATGTTBGMRCCCATMA  > rpob_nested_1st_PCR_then_2nd_PCR.ecopcr

```
Then, same treatment as above.

# Comparison with  16S
I have a complete database of 16S sequences extract from refseq genomes which use gff, which I can download using genome_updater.

For more information on this database: gabryelle.agoutin@inrae.fr

We need to keep the genomes in common between my rpob base and my 16S base so that they are comparable. 
I retrieve the folder names of the fetchMGS results and I retrieve the names of the genomes that have a non-empty fna file.

```bash=
l output_fetchMGS/*/*.fna > fna_files.txt
python find_fna_not_empty.py fna_files.txt  > genomes_names_rpoB.txt
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
zcat 16S.tsv.gz| cut -f11 | sort -u > genome_name_16S.txt
```

Look for the rpob genomes in the 16S genomes and if it's found it's good:
grep -F -f genome_names_rpob.txt  genome_name_16S.txt > common_genomes.txt

We therefore need to remove those that are not in common in the 16S fasta.
I'm looking for common between 16S and rpoB

```bash
grep -F -f genome_name_16S.txt genomes_names_rpoB.txt > common.txt
```

I retrieve the rows corresponding to these genes from my tsv:
```bash=
awk -F'\t'  'NR==FNR{genes[$1]; next} $11 in genes' common.txt 16S.tsv > output.tsv
gene_id seqid   product gene_name       start   end     strand  circular        seqid_length    partial genome_name     species_taxid
rna-HQ400_RS13095       NZ_CP053879.1   16S ribosomal RNA       16S     2835345 2836889 -       True    4535455 False   GCF_018802265.1 650
rna-HQ400_RS14170       NZ_CP053879.1   16S ribosomal RNA       16S     3049597 3051141 -       True    4535455 False   GCF_018802265.1 650
rna-HQ400_RS14960       NZ_CP053879.1   16S ribosomal RNA       16S     3225998 3227542 -       True    4535455 False   GCF_018802265.1 650
rna-HQ400_RS16380       NZ_CP053879.1   16S ribosomal RNA       16S     3513214 3514758 -       True    4535455 False   GCF_018802265.1 650
```

This will allow me to retrieve the seqid that will allow me to delete this sequence in my fasta using the delete_sequence_in_a_fasta.py script.

```bash=
cut -f2 output.tsv > seqid_com.txt
python keep_seq_in_fasta.py 16S.fna seqid_com.txt new_fasta_file.fasta
```

Formatting the new fasta so that it is compatible with EcoPCR with obiconvert 

```bash=
python format_fasta_for_ecopcr.py new_fasta_file.fasta new_fasta_file_format.fasta
mkdir ecoPCR_db/
obiconvert --fasta new_fasta_file_format.fasta --ecopcrdb-output=ecoPCR_db/16S -t ncbi_tax_dumb/2024-15-04/
```
if you ever get obiconvert errors due to an unknown ID, delete the sequences concerned. They are due to a lack of coherence between the different taxonomies.

Now going to use this fasta to launch my EcoPCR 16S.

```

| 16S_F  | 16S_R |
| -------- | -------- | 
| ACGGRAGGCAGCAG|TACCAGGGTATCTAATCCT|

```bash=
ecoPCR -d ecoPCR_db/16S -e 2  ACGGRAGGCAGCAG TACCAGGGTATCTAATCCT  > 16S.ecopcr
```

We retrieve the identifiers of the species amplified by the primers and the taxid in good format
```bash=
awk -F"|" '!/^#/{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $4); print ">" $1 "| taxid=" $4 ";"}' 16S.ecopcr > id_amplified.txt
```

We also retrieve the same list of identifiers from my input fasta file in order to have a list of all those that are not in common. So that's all the species that haven't been amplified.
```bash=
grep ">" new_fasta_file_format.fasta > name_seq.txt
```
To find out the difference between the two
```bash=
grep -v -F -f id_amplified.txt name_seq.txt > not_amplified.txt
```

Then turn the steps as above to make the krona graph
