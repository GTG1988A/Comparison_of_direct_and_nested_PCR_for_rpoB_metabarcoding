import os

# Path to the refseqid and taxid correspondence file
chemin_correspondances = "/path/result_correspondance_tab/tableau_resultat.tsv"

# Load matches from the file into a dictionary
correspondances = {}
with open(chemin_correspondances, 'r') as f:
    for ligne in f:
        elements = ligne.strip().split()
        if len(elements) == 2:
            refseq_id, taxid = elements
            correspondances[refseq_id] = int(taxid)

# Folder containing .faa and .fna files
dossier_parent = "/path/output_fetchMGS"

# Browse folders in the parent folder
for dossier in os.listdir(dossier_parent):
    chemin_dossier = os.path.join(dossier_parent, dossier)

    # Ensure that the item is a folder
    if os.path.isdir(chemin_dossier):
        # Extraire l'identifiant RefSeq du nom du dossier
        refseq_id = dossier.split("_")[0] + "_" + dossier.split("_")[1]

        # Retrieve the taxid associated with the file
        taxid = correspondances.get(refseq_id, None)

        if taxid is not None:
            for fichier in ["COG0085.faa", "COG0085.fna"]:
                chemin_fichier = os.path.join(chemin_dossier, fichier)

                with open(chemin_fichier, 'r') as f:
                    lignes = f.readlines()

                # Update the sequence name in each line
                for i in range(len(lignes)):
                    if lignes[i].startswith('>'):
                        lignes[i] = f">{lignes[i][1:].strip()} {taxid}\n"

                # Write the updated content to the file
                with open(chemin_fichier, 'w') as f:
                    f.writelines(lignes)

            print(f"The names of the sequences in {folder} have been updated with the taxid {taxid}.")
        else:
            print(f"Taxid not found for folder {folder}.")
