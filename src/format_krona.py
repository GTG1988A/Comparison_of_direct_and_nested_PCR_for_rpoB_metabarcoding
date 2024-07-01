import sys

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as f_in:
    with open(output_file, 'w') as f_out:
        for line in f_in:
            id_, taxonomy = line.strip().split('\t')
            taxonomy_levels = [level.split('__')[1] for level in taxonomy.split(';')]
            kingdom, phylum, classe, ordre, family, genus, species = taxonomy_levels
            new_line = f"{id_}\t1\t{species};{genus};{family};{ordre};{classe};{phylum};{kingdom}\tFalse\n"
            f_out.write(new_line)
