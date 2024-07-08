import argparse

def process_taxonomy(input_file, output_file, use_true):
    # Open the input file for reading
    with open(input_file, 'r') as f_in:
        # Open a new file for writing results
        with open(output_file, 'w') as f_out:
            # Read each line from the input file
            for line in f_in:
                # Split ID and taxonomy using tab delimiter
                id_, taxonomy = line.strip().split('\t')
                # Split taxonomy into different levels using ';'
                taxonomy_levels = [level.split('__')[1] for level in taxonomy.split(';')]
                # Extract required levels in the order: phylum, class, order, family, genus, species
                kingdom, phylum, classe, ordre, family, genus, species = taxonomy_levels
                # Write information in the new format to the output file
                if use_true:
                    new_line = f"{id_}\t1\t{species};{genus};{family};{ordre};{classe};{phylum};{kingdom}\tTrue\n"
                else:
                    new_line = f"{id_}\t1\t{species};{genus};{family};{ordre};{classe};{phylum};{kingdom}\tFalse\n"
                f_out.write(new_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process taxonomy data and write to output file.")
    parser.add_argument("input_file", type=str, help="Path to the input file")
    parser.add_argument("output_file", type=str, help="Path to the output file")
    parser.add_argument("use_true", type=str, choices=["True", "False"], help="Use 'True' or 'False' for the last field")

    args = parser.parse_args()

    use_true = args.use_true == "True"  # Convert string "True" or "False" to boolean

    process_taxonomy(args.input_file, args.output_file, use_true)
