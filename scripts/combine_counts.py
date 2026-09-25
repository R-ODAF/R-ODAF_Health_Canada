#!/usr/bin/env python3

import sys
import pandas as pd

# Get input files and output file from command line
md_file = sys.argv[1]
input_files = sys.argv[2:-1]  # All but last argument
output_file = sys.argv[-1] # Last argument only

# Read metadata
md = pd.read_csv(md_file, sep="\t")
SAMPLES = md["sample_ID"].tolist()

# Read all count tables
tables = []
for f in input_files:
    print(f"  Reading: {f}")
    df = pd.read_csv(f, sep="\t", index_col=0)
    tables.append(df)

# Combine with outer join (preserves all genes and samples)
# Fill missing values with 0
combined = pd.concat(tables, axis=1).fillna(0).astype(int)

# Keep only columns with names matching the samples in metadata
# Needed in cases where, ex. one library uses 96 barcodes and another only 48
# All barcodes will be searched for in both libraries by STARsolo, 
# and even barcodes that weren't used will get a few hits, creating new (not real!) columns
combined_filt = combined[SAMPLES]

# Sort by index (gene names) for consistency
combined_filt = combined_filt.sort_index()

print(f"Combined table: {combined_filt.shape[0]} genes x {combined_filt.shape[1]} samples")
print(f"Writing to: {output_file}")

# Write output
combined.to_csv(output_file, sep="\t")