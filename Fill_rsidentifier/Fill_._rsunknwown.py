#!/usr/bin/env python3

import gzip
import argparse

def process_vcf(input_vcf, output_vcf):
    """
    Process a gzipped vcf file to replace '.' in the ID column with a unique identifier.

    Args:
        input_vcf (str): Path to the gzipped input vcf file.
        output_vcf (str): Path to the gzipped output vcf file.
    """
    counter = 1  # Counter to generate unique IDs

    with gzip.open(input_vcf, 'rt') as infile, gzip.open(output_vcf, 'wt') as outfile:
        for line in infile:
            # Write header lines directly to the output
            if line.startswith("#"):
                outfile.write(line)
            else:
                columns = line.strip().split("\t")
                # Replace '.' in the ID column (column index 2) with a unique ID
                if columns[2] == ".":
                    columns[2] = f"rsunknown{counter}"
                    counter += 1
                # Write the modified line to the output file
                outfile.write("\t".join(columns) + "\n")
    
    print(f"Processing completed. Output saved to '{output_vcf}'.")

def main():
    parser = argparse.ArgumentParser(description="Replace '.' in the ID column of a vcf file with unique identifiers.")
    parser.add_argument("--input", "-i", required=True, help="Path to the gzipped input vcf file.")
    parser.add_argument("--output", "-o", required=True, help="Path to the gzipped output vcf file.")
    
    args = parser.parse_args()
    process_vcf(args.input, args.output)

if __name__ == "__main__":
    main()
