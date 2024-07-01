# Comparison of direct and nested PCR for rpoB metabarcoding


[TOC]

# rpoB
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
### Telechargement 

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

Mais pour que ça fonctionne, il faut que les deux fichier ai les mêmes identifiants de séquence 


### 2.a Preparing .faa files 
les fichiers sont sous la forme:

```bash=
$ head GCF_000219415.3_ASM21941v3_protein.faa
>WP_002964358.1 MULTISPECIES: 30S ribosomal protein S19 [Hyphomicrobiales]
MARSVWKGPFVDGYLLKKAEKVREGGRNEVIKMWSRRSTILPQFVGLTFGVYNGNKHVPVSVSEEMVGHKFGEFAPTRTYYGHGADKKAKRK
```

donc J'ai fais un script bash pour ne garder que l'identifiant protéique (WP_ ...)

```bash=
$ cat rename_faa.sh

#!/bin/bash

for file in /home/gagoutin/work/rpob_propre/protein_seq/arc_refseq_cg/2024-02-05_09-01-02/files/*.faa.gz; do

    filename=$(basename "$file" .faa.gz)

    zcat "$file" | awk '/^>/ {print $1; next} {print}' > "/home/gagoutin/work/rpob_propre/protein_seq/arc_refseq_cg/2024-02-05_09-01-02/files/modified/${filename}_modified.faa"
done

```

### 2.b Préparation fichier .fna

Ils sont en format:

```bash=
$ head GCF_000219415.3_ASM21941v3_cds_from_genomic.fna
>lcl|NZ_CP137016.1_cds_WP_244428070.1_1 [locus_tag=SFGR64A_RS00035] [protein=hypothetical protein] [protein_id=WP_244428070.1] [location=complement(6240..6629)] [gbkey=CDS]
GTGGTTGCCCTTCATGGCGCCGTTCCGCGCCGTTCCGCCGCGTCCTTCACCAATATCAAGGTTCGCGTGCGCGACGACCG(..)ACAGGAAGCGCTCGTCGCCAACATCGCCGTGGCGCTCAGAAAGCTCTCGATAGGCGAATTGTTGCTTTAG
```

Je ne garde que la section du protein_id avec le script:

```python=
#!/bin/bash

for file in PATH/fna/arc_refseq_cg/2024-02-29_13-59-50/files/*.fna.gz; do
    filename=$(basename "$file" .fna.gz)

    zcat "$file" | awk '/^>/ {
        if (match($0, /\[protein_id=([^]]*)\]/, arr)) {
            print ">" arr[1];
            skip_seq = 0;
        } else {
            skip_seq = 1;
        }
        next
    } !skip_seq {print}' > "PATH/fna/arc_refseq_cg/2024-02-29_13-59-50/files/modified/${filename}_modified.fna"
done

```

NB: Ce script ne garde pas les séquences si elle n'ont pas d'identifiant protéique car si il n'y à pas l'équivalent en protéine, FetchMGS ne tourne pas sur les acides nucléiques. Elles ne servent donc à rien. 

### 2.c Lancement fetch MGS

Il faut lancer fetch MGS pour environ 80 000 fichier. Donc si je donne la ligne de commande tel quel : 

```bash=
./fetchMGs/fetchMGs.pl -m extraction test/arc_refseq_cg/2024-01-08_11-31-26/files/*.faa -d arc_refseq_cg/2024-01-08_07-45-09/files/*.fna -c COG0085
```
Je vais avoir une erreur car il y a trop d'entrée. 

J'ai tout de même essayé en lançant sur seulement 2 fichiers mais en fait fetch n'en gère qu'un à la fois.

Je récupère les noms des fichiers: 

```
ls cds_genomic_seq/arc_refseq_cg/2024-02-05_09-01-16/files/modified/ > tmp.txt
sed 's/_cds_from_genomic\.fna\.gz_modified\.fna//' > noms_des_fichiers.txt
rm tmp.txt
```

et je fais un script qui fait une boucle pour faire ligne de commande par fichier 


```bash=
#!/bin/bash

# Chemin vers le fichier contenant la liste des noms de génomes
liste_noms="nom_genomes.txt"

# Vérification de l'existence du fichier contenant la liste des noms de génomes
if [ ! -f "$liste_noms" ]; then
    echo "Le fichier $liste_noms n'existe pas."
    exit 1
fi


mkdir -p .output_fetchMGS

# Chemin vers les répertoires
repertoire_faa="/PATH/rpob/protein_seq/arc_refseq_cg/2024-02-05_09-01-02/files/modified/"
repertoire_fna="/PATH/rpob/cds_genomic_seq/arc_refseq_cg/2024-02-05_09-01-16/files/modified/"
chemin_sortie_global="output_fetchMGS/"

# Lecture de la liste de noms de génomes
while IFS= read -r genome; do
    # Construction des chemins des fichiers FAAs et FNAs
    fichier_faa="${repertoire_faa}${genome}_protein_modified.faa"
    fichier_fna="${repertoire_fna}${genome}_cds_from_genomic.fna.gz_modified.fna"

    # Vérification de l'existence du fichier FNA
    if [ -e "$fichier_fna" ]; then
        # Création du dossier de sortie
        dossier_sortie="${chemin_sortie_global}${genome}_output"

        # Appel de fetchMGs.pl avec les chemins des fichiers FAAs et FNAs
        ./../fetchMGs/fetchMGs.pl -m extraction "$fichier_faa" -d "$fichier_fna" -c COG0085 -o "$dossier_sortie"
    else
        echo "Fichier .fna manquant pour le génome $genome"
    fi

```

C'est assez long, mais c'est fetch et c'est à faire qu'une fois


Fetch crée un dossier par génome avec dedans les fichiers COG0085.faa et COG0085.fna. Le COG0085 correspond à rpoB c'est celui qui nous intéresse. 

## 3. Formatage du fichier

### 3.a Ajout du Taxid
Je me suis souvenu qu'il fallait ajouter le taxid pour lancer obiconvert avant EcoPCR. Je fais un tableau de correspondance entre l'identifiant du génome et son species taxids qui sont les colonnes 1 et 7 de l'assembly summary téléchargé par genome_updater:

```
cut -f1,7 genome/arc_refseq_cg/2024-02-05_09-00-49/assembly_summary.txt > correspondance_tab.txt
```



Grâce à ce tableau de correspondance, j'ai fais un script python qui va aller dans chaque dossier résultant de FetchMGS; trouvé l'ID et le taxid correspondant puis il l'ajoutera au nom de la sequence a coter de l'ID protéique

```python=
import os

# Chemin vers le fichier des correspondances refseqid et taxid
chemin_correspondances = "~/work/rpob_propre/result_correspondance_tab/tableau_resultat.tsv"

# Charger les correspondances depuis le fichier dans un dictionnaire
correspondances = {}
with open(chemin_correspondances, 'r') as f:
    for ligne in f:
        elements = ligne.strip().split()
        if len(elements) == 2:
            refseq_id, taxid = elements
            correspondances[refseq_id] = int(taxid)

# Dossier contenant les fichiers .faa et .fna
dossier_parent = "/home/gagoutin/work/rpob_propre/output_fetchMGS"

# Parcourir les dossiers dans le dossier parent
for dossier in os.listdir(dossier_parent):
    chemin_dossier = os.path.join(dossier_parent, dossier)

    # S'assurer que l'élément est un dossier
    if os.path.isdir(chemin_dossier):
        # Extraire l'identifiant RefSeq du nom du dossier
        refseq_id = dossier.split("_")[0] + "_" + dossier.split("_")[1]

        # Récupérer le taxid associé au dossier
        taxid = correspondances.get(refseq_id, None)

        if taxid is not None:
            # Parcourir les fichiers .faa et .fna
            for fichier in ["COG0085.faa", "COG0085.fna"]:
                chemin_fichier = os.path.join(chemin_dossier, fichier)

                # Lire le contenu du fichier
                with open(chemin_fichier, 'r') as f:
                    lignes = f.readlines()

                # Mettre à jour le nom de la séquence dans chaque ligne
                for i in range(len(lignes)):
                    if lignes[i].startswith('>'):
                        lignes[i] = f">{lignes[i][1:].strip()} {taxid}\n"

                # Écrire le contenu mis à jour dans le fichier
                with open(chemin_fichier, 'w') as f:
                    f.writelines(lignes)

            print(f"Les noms des séquences dans {dossier} ont été mis à jour avec le taxid {taxid}.")
        else:
            print(f"Taxid non trouvé pour le dossier {dossier}.")

```


### 3.b Concaténation des fichiers fna pour avoir nos régions rpob 

Il faut en suite les concaténer pour avoir le fichier  complet avec les gènes rpob en acide nucléique

```bash=
find output/ -name "COG0085.fna" -exec cat {} + > concatenated_COG0085.fna
```


### Nettoyage
Je me suis rendu compte que parfois, on se retrouve avec des séquences qui n'ont pas de nom de protéine. Dans ce genre de cas, fetch remplace le nom par ">1" mais il n'y à pas de résultat faa. Donc si on tombe sur ce cas, il faut enlever les séquences avec cette commande:

```bash=
awk '/^>1/{getline;next} 1' concatenated_COG0085.fna > output_file.fna
``` 

Mais pour les génomes à tout niveau d'assemblage, il n'y avait pas ce cas au vu du nettoyage des fasta avant (ceux qui n'ont pas de protein_id)

=> Donc maintenant j'ai un fichier fasta de séquence nucléique de tout les gènes rpoB des génomes de refseq. 

## 4. Formatage des bases pour EcoPCR: Obitools
### Obitools pour les séquences rpod

```bach!
$ ml devel/Miniconda/Miniconda3
$ ml bioinfo/OBITools/1.2.11
(obitools-v1.2.11_env)  $ obiconvert --fasta concatenated_COG0085.fna --ecopcrdb-output=ecoPCR_db/rpob -t ncbi_tax_dumb/

```

NB : il est possible d'avoir une erreur car parfois genome updater à des ID en plus que la taxonomy du ncbi qui n'a pas encore été mise à jour. Ici, je n'ai pas eu cette erreur. 


## 5. Lancement EcoPCR

### 1ère PCR pour la Nested 

| rpoB_F  | rpoB_R |
| -------- | -------- | 
| CAGYTDTCNCARTTYATGGAYCA|AGTTRTARCCDTYCCANGKCAT| 


```bash=
$ module load bioinfo/ecoPCR/1.0.1;
$ ecoPCR -d ecoPCR_db/rpob -e 2 CAGYTDTCNCARTTYATGGAYCA AGTTRTARCCDTYCCANGKCAT  > primers/rpob_l70_L8000_e2.ecopcr

```

Je récupère les noms des génomes qui ont été amplifiés, ceux qui sont présents dans le résultat de EcoPCR. J'en récupère le species taxid, puis je récupère leur taxonomie.
Ceux qui ne sont pas dans EcoPCR sont aussi récupérés, ce sont ceux qui n'ont pas été amplifiés. Je récupère aussi leur species ID.
S'ils sont dans le résultat de EcoPCR, ils ont un "True" dans la colonne "amplified", et s'ils n'y sont pas, ils ont un "False".
Ensuite, je lance le script pour le Krona

```bash=
#Recupère les ID protéiques des amplifiés
awk '!/^#/{print $1}' FS="|" ecoPCR_rpob1/rpob_l70_L8000_e2.ecopcr > id_amplified.txt

# Les ID total dans tout les génomes
grep ">" concatenated_COG0085.fna | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt

#Pour avoir la différence entre les deux
grep -v -F -f liste_id_ecoPCR_sort.txt id_amplified.txt > not_amplified.txt
```

J'ai les identifiants de séquence maintenant je veux les taxids:

```bash=
#mon fichier avec les nom de sequences uniquement
grep ">" concatenated_COG0085.fna > name_seq.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' not_amplified.txt name_seq.txt > tmp_no_amplified.txt
cut -d'|' -f2 tmp_no_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > no_amplified_speciesTAXID.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' id_amplified.txt name_seq.txt > tmp_amplified.txt
cut -d'|' -f2 tmp_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > amplified_speciesTAXID.txt

#pas très propre tout ces cut mais pas très grave
```
ça me permet de me servir de l'outil TaxonKit pour avoir leurs taxonomies:
```bash=
$ ml bioinfo/TaxonKit/0.15.0
# non amplifiéeead
$ cat no_amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > no_amplified_species.taxo
gagoutin@n001 ~/work/rpob/chromosome/ecoPCR_nested $ head -n1 no_amplified_species.taxo
210     k__Bacteria;p__Campylobacterota;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori

# amplifiée
$ cat amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > amplified_species.taxo

```

J'enlève tout les possibles espace dans les noms:

```bash=
$ sed 's/ /_/g' amplified_species.taxo > amplified_species_wc.taxo
$ sed 's/ /_/g' no_amplified_species.taxo > no_amplified_species_wc.taxo
```
Je reformate ces fichiers pour qu'ils soit au format accepter par krona:
* Je rajoute une colonne "Sequence" qui correspond à une somme de séquence. Elles sont toute à 1
* Je rajoute une colonne "amplified" qui est a True ou False. Donc True pour le fichier des amplifiées et False pour le fichier des non amplifiées

NB: Je change le True et False à la main dans le script car je n'ai pas ajouté l'argument mais je peux le faire au besoin

```python=
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Ouvrir le fichier d'entrée en mode lecture
with open(input_file, 'r') as f_in:
    # Ouvrir un nouveau fichier en mode écriture pour stocker les résultats
    with open(output_file, 'w') as f_out:
        # Lire chaque ligne du fichier d'entrée
        for line in f_in:
            # Séparer l'ID et la taxonomie en utilisant l'onglet comme délimiteur
            id_, taxonomy = line.strip().split('\t')
            # Diviser la taxonomie en différents niveaux en utilisant ';'
            taxonomy_levels = [level.split('__')[1] for level in taxonomy.split(';')]
            # Extraire les niveaux requis dans l'ordre phylum, class, ordre, famille, genre, espèce
            kingdom, phylum, classe, ordre, family, genus, species = taxonomy_levels
            # Écrire les informations dans le nouveau format dans le fichier de sortie
            new_line = f"{id_}\t1\t{species};{genus};{family};{ordre};{classe};{phylum};{kingdom}\tFalse\n"
            f_out.write(new_line)

```



Puis je crée un header:
```bash=
$ cat header.tsv
seqID	sequence	taxonomy	amplified
```

Je concatene tout:

```bash=
$ cat header.tsv amplified_species_wc_reformat.tsv > tmp.tsv
$ cat tmp.tsv no_amplified_species_wc_reformat.tsv > final_results.tsv
```

J'enlève la colonne correspondant au taxid car ce n'est pas accepter par le script de krona:

```bash=
cut -f 2-4 final_results.tsv > input.tsv
```

Je lance le script me permettant d'avoir le krona:
```bash=
$ python make_krona_xml_bis.py input.tsv output.xml  --name rpob -v --color_label amplified  --dataset_labels sequence
$ ml bioinfo/Krona/2.8.1
$ ktImportXML output.xml

```

### 2nde PCR  (amorces de metabarcoding)
| Univ_rpoB_F_deg  | Univ_rpoB_R_deg |
| -------- | -------- | 
| GGYTWYGAAGTNCGHGACGTDCA|TGACGYTGCATGTTBGMRCCCATMA| 


```bash=
$ module load bioinfo/ecoPCR/1.0.1;
$ ecoPCR -d ecoPCR_db/rpob -e 2 GGYTWYGAAGTNCGHGACGTDCA TGACGYTGCATGTTBGMRCCCATMA  > primers/rpob_l70_L8000_e2.ecopcr

```

=> même traitement que ci dessus

![image](https://hackmd.io/_uploads/HyLgl96Tp.png)

### 1ere puis 2nde PCR

Pour celui-ci j'ai donc pris les résultats d'EcoPCR du premier couple d'amorce appelé rpoB_F et rpoB_R
Puis j'ai garder uniquement les génomes qui apparaissent dans le résultat de EcoPCR, donc quand les amorces ont fonctionnée.

```bash!
#Recupère les ID protéiques des amplifiés
awk '!/^#/{print $1}' FS="|" ecoPCR_rpob1/rpob_l70_L8000_e2.ecopcr > id_amplified.txt

# Les ID total dans tout les génomes
grep ">" concatenated_COG0085.fna | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt

#Pour avoir la différence entre les deux
grep -v -F -f liste_id_ecoPCR_sort.txt liste_id_concatened_COG0085_sort.txt > not_amplified.txt
wc -l not_amplified.txt #7 481
```

J'enlève donc 7 481 génomes avec ce script:

```python=
from Bio import SeqIO

def filter_sequences(fasta_file, id_list_file, output_file):
    # Charger les identifiants de séquence depuis le fichier
    with open(id_list_file, 'r') as f:
        id_list = set(line.strip() for line in f)

    # Parcourir le fichier FASTA et écrire les séquences dans un nouveau fichier
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Extraire l'identifiant de séquence complet du titre du FASTA
            seq_id = record.description.strip()

            # Vérifier si l'identifiant est dans la liste à filtrer
            if seq_id.split('|')[0] not in id_list:
                out_f.write('>' + record.description + '\n')
                out_f.write(str(record.seq) + '\n')
if __name__ == "__main__":
    fasta_file = "concatenated_COG0085.fna"
    id_list_file = "Wp_id_not_amplified.txt"
    output_file = "fichier_filtre.fasta"

    filter_sequences(fasta_file, id_list_file, output_file)


```

```bash=
python suppression_sequence_dans_un_fasta.py
```

Maintenant je vais reformater le fasta avec Obitools pour lancer EcoPCR dessus

```bash=
$ ml devel/Miniconda/Miniconda3
$ ml bioinfo/OBITools/1.2.11
$ obiconvert --fasta fichier_filtre.fasta --ecopcrdb-output=ecoPCR_db/rpob -t ../ncbi_tax_dumb/2024-15-02/
```

Lancement ecoPCR avec les 2nd amorces

| Univ_rpoB_F_deg  | Univ_rpoB_R_deg |
| -------- | -------- | 
| GGYTWYGAAGTNCGHGACGTDCA|TGACGYTGCATGTTBGMRCCCATMA| 

```bash=
$ module load bioinfo/ecoPCR/1.0.1;
$ ecoPCR -d ecoPCR_db/rpob -e 2  GGYTWYGAAGTNCGHGACGTDCA TGACGYTGCATGTTBGMRCCCATMA  > primers/rpob_l70_L8000_e2.ecopcr

```

Je récupère les identifiants des espèces amplifiées par les amorces

```bash=
$ awk '!/^#/{print $1}' FS="|" primers/rpob_l70_L8000_e2.ecopcr > id_seq_amplified_species.txt
$ wc -l id_seq_amplified_species.txt

```
Il y a 48823 séquences amplifiées sur les 71 525


Je récupère aussi une même liste d'identifiant de mon fichier fasta d'entrée afin d'avoir la liste de tout ceux qui ne sont pas en commun. C'est donc toute les espèces qui n'ont pas été amplifié

```bash=
grep ">" fichier_filtre.fasta | cut -d'|' -f1 |cut -d '>' -f2 > id_seq_total.txt
# pas propre pour enlever le > mais awk ajoutait des lignes vides alors j'ai fais ça rapidement. ça fonctionne très bien.
```

pour avoir la liste des non amplifiés:
```bash=
grep -v -F -f id_seq_amplified_species.txt id_seq_total.txt > id_no_amplified_species.txt
```


Je récupère les taxid de ceux qui ont été amplifié ainsi que des non amplifiés:

```bash=
#mon fichier avec les nom de sequences uniquement
grep ">" concatenated_COG0085.fna > name_seq.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' id_no_amplified_species.txt name_seq.txt > tmp_no_amplified.txt
cut -d'|' -f2 tmp_no_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > no_amplified_speciesTAXID.txt
awk 'NR==FNR {ids[$1]; next} {for (id in ids) {if (index($0, id)) {print; break}}}' id_seq_amplified_species.txt name_seq.txt > tmp_amplified.txt
cut -d'|' -f2 tmp_amplified.txt | cut -d'=' -f2 |cut -d ';' -f1  > amplified_speciesTAXID.txt

# encore une fois pas très propre tout ces cut mais pas très grave
```

ça me permet de me servir de l'outil TaxonKit pour avoir leurs taxonomies:

```bash=
$ ml bioinfo/TaxonKit/0.15.0
# non amplifiéeead
$ cat no_amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > no_amplified_species.taxo
gagoutin@n001 ~/work/rpob/chromosome/ecoPCR_nested $ head -n1 no_amplified_species.taxo
210     k__Bacteria;p__Campylobacterota;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori

# amplifiée
$ cat amplified_speciesTAXID.txt | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-02/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-02/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > amplified_species.taxo

```

J'enlève tout les possibles espace dans les noms:

```bash=
$ sed 's/ /_/g' amplified_species.taxo > amplified_species_wc.taxo
$ sed 's/ /_/g' no_amplified_species.taxo > no_amplified_species_wc.taxo
```

Je reformate ces fichiers pour qu'ils soit au format accepter par krona:
* Je rajoute une colonne "Sequence" qui correspond à une somme de séquence. Elles sont toute à 1
* Je rajoute une colonne "amplified" qui est a True ou False. Donc True pour le fichier des amplifiées et False pour le fichier des non amplifiées

NB: Je change le True et False à la main dans le script car je n'ai pas ajouté l'argument mais je peux le faire au besoin

```python=
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Ouvrir le fichier d'entrée en mode lecture
with open(input_file, 'r') as f_in:
    # Ouvrir un nouveau fichier en mode écriture pour stocker les résultats
    with open(output_file, 'w') as f_out:
        # Lire chaque ligne du fichier d'entrée
        for line in f_in:
            # Séparer l'ID et la taxonomie en utilisant l'onglet comme délimiteur
            id_, taxonomy = line.strip().split('\t')
            # Diviser la taxonomie en différents niveaux en utilisant ';'
            taxonomy_levels = [level.split('__')[1] for level in taxonomy.split(';')]
            # Extraire les niveaux requis dans l'ordre phylum, class, ordre, famille, genre, espèce
            kingdom, phylum, classe, ordre, family, genus, species = taxonomy_levels
            # Écrire les informations dans le nouveau format dans le fichier de sortie
            new_line = f"{id_}\t1\t{species};{genus};{family};{ordre};{classe};{phylum};{kingdom}\tFalse\n"
            f_out.write(new_line)

```



Puis je crée un header:
```bash=
$ cat header.tsv
seqID	sequence	taxonomy	amplified
```

Je concatene tout:

```bash=
$ cat header.tsv amplified_species_wc_reformat.tsv > tmp.tsv
$ cat tmp.tsv no_amplified_species_wc_reformat.tsv > final_results.tsv
```

J'enlève la colonne correspondant au taxid car ce n'est pas accepter par le script de krona:

```bash=
cut -f 2-4 final_results.tsv > input.tsv
```

Je lance le script me permettant d'avoir le krona:
```bash=
$ python make_krona_xml_bis.py input.tsv output.xml  --name rpob -v --color_label amplified  --dataset_labels sequence
$ ml bioinfo/Krona/2.8.1
$ ktImportXML output.xml

```

![image](https://hackmd.io/_uploads/H1uXl5pa6.png)


## 6. EcoPCR avec 10 représentant par espèces
Nous avions peur que certaine espèces soit sur representé ce qui pourrait pousser à un mauvaise interpréation. J'ai donc réduit notre fichier fasta avec seulement 10 répresentant par espèce. 
Chaque génome à un species id, donc au bout du 10 eme species ID identique, si il y avait d'autre séquence avec ce speciesID là elle ne seront pas conservé.

script:
```python=
#!/bin/bash

# Nom du fichier d'entrée
input_file="concatenated_COG0085.fna"

# Créer un tableau associatif pour stocker le nombre de représentants par taxid
declare -A taxid_count

# Variable pour stocker l'identifiant de séquence
seq_id=""

# Lire le fichier ligne par ligne
while IFS= read -r line; do
    # Vérifier si la ligne contient un identifiant de séquence
    if [[ $line == ">"* ]]; then
        # Extraire le taxid de l'identifiant de séquence
        taxid=$(echo "$line" | grep -oP 'taxid=\K\d+')
        # Vérifier si le nombre de représentants pour ce taxid est inférieur à 10
        if (( taxid_count[$taxid] < 10 )); then
            # Afficher l'identifiant de la séquence
            echo "$line"
            # Stocker l'identifiant de la séquence actuelle
            seq_id="$line"
            # Incrémenter le compteur pour ce taxid
            (( taxid_count[$taxid]++ ))
        else
            # Réinitialiser l'identifiant de la séquence pour ne pas afficher la séquence suivante
            seq_id=""
        fi
    else
        # Afficher la séquence si l'identifiant de la séquence actuelle a été stocké
        if [[ -n $seq_id ]]; then
            echo "$line"
        fi
    fi
done < "$input_file"
```

Une fois ceci réalisé je lance exactement le même traitement que la partie ci-dessus. 
