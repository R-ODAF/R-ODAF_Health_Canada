#!/usr/bin/env python3

import sys
import pandas as pd

# Get input files and output file from command line
input_files = sys.argv[1:-1]  # All but last argument
output_file = sys.argv[-1] # Last argument only

# Read all tables
tables = []
for f in input_files:
    print(f"  Reading: {f}")
    df = pd.read_csv(f, sep="\t", index_col=0)
    tables.append(df)

# Combine with outer join (preserves all genes and samples)
# Fill missing values with 0
combined = pd.concat(tables, axis=1).fillna(0).astype(int)

# Sort by index (gene names) for consistency
combined = combined.sort_index()

print(f"Combined table: {combined.shape[0]} genes x {combined.shape[1]} samples")
print(f"Writing to: {output_file}")

# Write output
combined.to_csv(output_file, sep="\t")