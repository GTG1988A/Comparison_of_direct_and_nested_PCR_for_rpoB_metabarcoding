import re
import argparse

def find_non_empty_genomes(input_file):
    capture_assembly_acc = re.compile(r'([A-Z]{3}_\d+\.\d+)')
    with open(input_file, 'r') as file:
        for line in file:
            columns = line.split()
            if columns[4] != '0':
                match = capture_assembly_acc.search(line)
                if match:
                    genome_id = match.group(1)
                    print(genome_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find non-empty genomes based on a formatted input file.")
    parser.add_argument("input_file", type=str)

    args = parser.parse_args()

    find_non_empty_genomes(args.input_file)

