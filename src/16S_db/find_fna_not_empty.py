import re
capture_assembly_acc = re.compile(r'([A-Z]{3}_\d+\.\d+)')
with open('fna_files.txt', 'r') as file:
    for line in file:
        columns = line.split()
        if columns[4] != '0':
            match = capture_assembly_acc.search(line)
            if match:
                genome_id = match.group(1)
                print(genome_id)
