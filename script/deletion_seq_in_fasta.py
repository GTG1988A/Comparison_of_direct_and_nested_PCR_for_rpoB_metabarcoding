import argparse
from Bio import SeqIO

def filter_sequences(fasta_file, id_list_file, output_file):
    # Load sequence names from the file
    with open(id_list_file, 'r') as f:
        id_list = set(line.strip() for line in f)

    # Process the FASTA file and write filtered sequences to a new file
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Extract the full sequence ID from the FASTA header
            seq_id = record.description.strip()

            # Check if the full sequence ID is in the filter list
            if '>' + seq_id in id_list:
                out_f.write('>' + record.description + '\n')
                out_f.write(str(record.seq) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences from a FASTA file based on a list of sequence names.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file")
    parser.add_argument("id_list_file", type=str, help="Path to the file containing sequence names to filter")
    parser.add_argument("output_file", type=str, help="Path to the output file for filtered sequences")

    args = parser.parse_args()

    filter_sequences(args.fasta_file, args.id_list_file, args.output_file)
