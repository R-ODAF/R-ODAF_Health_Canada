
# Example use: python biospyder_manifest_to_reffiles.py -i your_manifest.csv -p Probe_Name -s Probe_Sequence
# By default, output file names will match the manifest file name
# You may use --fasta-output and --gtf_output flags to override the default output filenames
	
import csv
import sys
import os 
from argparse import ArgumentParser

def main():
    parser = ArgumentParser(
        description="Takes a BioSpyder manifest file (CSV) and outputs FASTA and GTF reference files for the BioSpyder probes in that kit."
    )
    parser.add_argument('-i', '--manifest', required=True,
                        help='Path to the input manifest CSV file.')
    parser.add_argument('-p', '--probe_ID_column', required=True,
                        help='Name of the column containing the identifier to be used for FASTA headers and GTF entries (e.g., Probe_Name)
    parser.add_argument('-s', '--sequence_column', required=True,
                        help='Name of the column containing the probe sequences (e.g., Probe_Sequence).')

    # Get manifest file name, use for default output file names
    known_args, _ = parser.parse_known_args()
    manifest_basename = os.path.basename(known_args.manifest)
    manifest_stem, _ = os.path.splitext(manifest_basename)

    default_fasta_output = f"{manifest_stem}.fa"
    default_gtf_output = f"{manifest_stem}.gtf"
	
	# Allow output filename overrides from arguments
    parser.add_argument('-f', '--fasta_output', default=default_fasta_output,
                        help=f'Name for the output FASTA file (default: {default_fasta_output}).')
    parser.add_argument('-g', '--gtf_output', default=default_gtf_output,
                        help=f'Name for the output GTF file (default: {default_gtf_output}).')

    args = parser.parse_args() # Re-parse all arguments, now with updated defaults

    try:
        # Open input CSV and output FASTA/GTF files
        with open(args.manifest, 'r', newline='', encoding='utf-8') as infile, \
             open(args.fasta_output, 'w', encoding='utf-8') as fasta_outfile, \
             open(args.gtf_output, 'w', encoding='utf-8') as gtf_outfile:

            reader = csv.DictReader(infile)

            # Check if the specified columns exist in the CSV header
            if args.probe_ID_column not in reader.fieldnames:
                sys.exit(f"Error: Column '{args.probe_ID_column}' not found in the manifest file header. "
                         f"Please ensure you provide the correct column name for the probe identifier (e.g., 'Probe_Name' from your example).")
            if args.sequence_column not in reader.fieldnames:
                sys.exit(f"Error: Column '{args.sequence_column}' not found in the manifest file header.")

            # Process each row in the CSV
            for row in reader:
                probe_identifier = row[args.probe_ID_column]
                probe_sequence = row[args.sequence_column]

                # Write to FASTA file
                fasta_outfile.write(f">{probe_identifier}\n")
                fasta_outfile.write(f"{probe_sequence}\n")

                # Write to GTF file
                sequence_length = len(probe_sequence)
                gtf_line = (
                    f"{probe_identifier}\tprotein_coding\texon\t1\t{sequence_length}\t.\t.\t.\t"
                    f'gene_id "{probe_identifier}";transcript_id "{probe_identifier}";\n'
                )
                gtf_outfile.write(gtf_line)

        print(f"Successfully generated {args.fasta_output} and {args.gtf_output} from '{args.manifest}'.")

    except FileNotFoundError:
        sys.exit(f"Error: The manifest file '{args.manifest}' was not found. Please check the path and filename.")
    except Exception as e:
        sys.exit(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()