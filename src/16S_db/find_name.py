import re
import argparse
import subprocess

# Regular expression to capture the genome name
capture_assembly_acc = re.compile(r'([A-Z]{3}_\d+\.\d+)')

def extract_species_taxid(assembly_summary_file):
    genome_species_map = {}
    with open(assembly_summary_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 7:
                genome_name = fields[0]
                species_taxid = fields[6]
                genome_species_map[genome_name] = species_taxid
    return genome_species_map

def main(assembly_summary_file, genome_files_output):
    genome_species_map = extract_species_taxid(assembly_summary_file)
    # Run the find command to obtain the file paths
    try:
        command = ['find', '-L', '/PATH/downloads/2024-05-15_11-08-32/files/', '-mindepth', '5', '-maxdepth', '5', '-type', 'f']
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        output_lines=[]
        for line in process.stdout:
            line=line.decode().strip()
            output_lines.append(line)
        process.wait()
        output, _ = process.communicate()
    except subprocess.CalledProcessError as e:
        print("Error:", e)
    # Open output file
    with open(genome_files_output, 'w') as output_file:
        for line in output_lines:
            genome_name_match = capture_assembly_acc.search(line)
            if genome_name_match:
                genome_name = genome_name_match.group(1)
                if genome_name in genome_species_map:
                    species_taxid = genome_species_map[genome_name]
                    output_file.write(f"{line.strip()} {genome_name} {species_taxid}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract genome files and species taxid from assembly summary.')
    parser.add_argument('assembly_summary_file', type=str, help='Path to the assembly summary file')
    parser.add_argument('genome_files_output', type=str, help='Path to the output file containing genome files and species taxid')
    args = parser.parse_args()

    main(args.assembly_summary_file, args.genome_files_output)
