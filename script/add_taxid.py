import os
import argparse

def update_sequence_names(correspondence_path, parent_folder):
    # Load correspondences from the file into a dictionary
    correspondences = {}
    with open(correspondence_path, 'r') as f:
        for line in f:
            elements = line.strip().split()
            if len(elements) == 2:
                refseq_id, taxid = elements
                correspondences[refseq_id] = int(taxid)

    # Iterate through folders in the parent folder
    for folder in os.listdir(parent_folder):
        folder_path = os.path.join(parent_folder, folder)

        # Ensure the element is a folder
        if os.path.isdir(folder_path):
            # Extract the RefSeq identifier from the folder name
            refseq_id = folder.split("_")[0] + "_" + folder.split("_")[1]

            # Retrieve the associated taxid for the folder
            taxid = correspondences.get(refseq_id, None)

            if taxid is not None:
                # Iterate through .faa and .fna files
                for file in ["COG0085.faa", "COG0085.fna"]:
                    file_path = os.path.join(folder_path, file)

                    # Read the file content
                    with open(file_path, 'r') as f:
                        lines = f.readlines()

                    # Update the sequence name in each line
                    for i in range(len(lines)):
                        if lines[i].startswith('>'):
                            lines[i] = f">{lines[i][1:].strip()}| taxid={taxid};\n"
                    with open(file_path, 'w') as f:
                        f.writelines(lines)

                print(f"Sequence names in {folder} have been updated with taxid {taxid}.")
            else:
                print(f"Taxid not found for folder {folder}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update sequence names in .faa and .fna files with taxid.")
    parser.add_argument("correspondence_path", type=str, help="Path to the correspondence file")
    parser.add_argument("parent_folder", type=str, help="Path to the parent folder containing subfolders with .faa and .fna files")

    args = parser.parse_args()

    update_sequence_names(args.correspondence_path, args.parent_folder)
