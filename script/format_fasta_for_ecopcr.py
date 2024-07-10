import argparse

def format_fasta(input_file, output_file):
    """
    Format FASTA file identifiers to limit the total length of line_number_prefix to 19 characters and add number of sequence.

    Args:
    - input_file (str): Path to input FASTA file.
    - output_file (str): Path to output formatted FASTA file.
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        line_number = 1
        
        for line in f_in:
            if line.startswith('>'):
                # Split identifier into parts before and after "|"
                parts = line.strip().split('|', 1)
                
                # Get the first part before "|"
                identifier = parts[0][1:]  # Remove ">" at the beginning
                
                # Construct line_number_prefix
                line_number_prefix = f"{line_number}_{identifier}"
                
                # Limit total length to 19 characters
                max_line_number_prefix_length = 19
                if len(line_number_prefix) > max_line_number_prefix_length:
                    line_number_prefix = line_number_prefix[:max_line_number_prefix_length]
                
                # Rebuild the line with the formatted identifier and the rest of the sequence description
                if len(parts) > 1:
                    formatted_line = f">{line_number_prefix}|{parts[1]}\n"
                else:
                    formatted_line = f">{line_number_prefix}\n"
                
                # Write to the output file
                f_out.write(formatted_line)
                
                # Increment line number
                line_number += 1
            else:
                # Write sequence lines without modification
                f_out.write(line)

def main():
    parser = argparse.ArgumentParser(description='Format FASTA file identifiers.')
    parser.add_argument('input_file', help='Input FASTA file path')
    parser.add_argument('output_file', help='Output formatted FASTA file path')
    
    example_command = 'python format_fasta.py input.fasta output_formatted.fasta'
    parser.epilog = f'Example command: {example_command}'
    
    args = parser.parse_args()
    
    format_fasta(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

