from Bio import SeqIO

def filter_sequences(fasta_file, id_list_file, output_file):
    with open(id_list_file, 'r') as f:
        id_list = set(line.strip() for line in f)

    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Extract the full sequence identifier from the FASTA title
            seq_id = record.description.strip()

            # Check whether the ID is in the list to be filtered
            if seq_id.split('|')[0] not in id_list:
                out_f.write('>' + record.description + '\n')
                out_f.write(str(record.seq) + '\n')
if __name__ == "__main__":
    fasta_file = "16S.fna"
    id_list_file = "seqid_not_com.txt"
    output_file = "new_fasta_file.fasta"

    filter_sequences(fasta_file, id_list_file, output_file)
