# This script produces a summary table from benchmark files output by Snakemake
## when running the R-ODAF_Health_Canada pipeline
#Run it in the directory with the benchmark files

import pandas as pd
import glob
import os
from argparse import ArgumentParser

# Arguments
parser = ArgumentParser(
    description = 
    '''
    Produces a summary table from all the benchmark files in the directory where it's run.
    Use the -c flag to enter columns of interest, separated by spaces.
    example usage: python3 scripts/summarize_sm_benchmarks.py -i output/logs/ -c s mean_load -o output/logs/summary.txt -t test
    ''')
parser.add_argument('-i', '--input', nargs = 1, required = True,
                    help='Path to the directory containing benchmark files.')
parser.add_argument('-c', '--columns', nargs = '+', required = True,
    choices={"s", "max_rss", "max_vms", "max_uss", "max_pss", "io_in", "io_out", "mean_load", "cpu_time"},
    help = 'Columns for which you want summary data.')
parser.add_argument('-o', '--output_file', help='Path to the output TSV file', required=True)
parser.add_argument('-t', '--tag', required = True, type = str,
    help='Unique identifier for the dataset/infrastructure/environment. To facilitate comparison across different runs.')
args = parser.parse_args()

# Add trailing slash to input path if needed
if args.input[0].endswith('/') == False:
	args.input[0] = args.input[0] + '/'

# Get file names
files_pattern = os.path.join(args.input[0], 'benchmark.*.txt')

# Use glob to get a list of file paths matching the pattern
file_paths = glob.glob(files_pattern)

# Initialize an empty list to store DataFrames
dfs = []

# Loop through each file and read it into a DataFrame
for file_path in file_paths:
    # Extract the part of the file name between "benchmark." and ".txt"
    info_from_file_name = os.path.basename(file_path).split('benchmark.')[1].replace('.txt', '')
    
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')
    
    # Add a new column with information from the file name
    df.insert(0, 'source', info_from_file_name)
    
    # Split the 'source' column and create a new 'rule' column
    split_source = df['source'].str.split('.', n=1, expand=True)
    df['rule'] = split_source[1] if len(split_source.columns) > 1 else df['source']
        
    # Append the DataFrame to the list
    dfs.append(df)

# Concatenate all DataFrames in the list
combined_df = pd.concat(dfs, axis=0, ignore_index=True)

# Reset index
combined_df = combined_df.reset_index(drop=True)

grouped_df = combined_df.groupby('rule')

# Create an empty DataFrame to store the results
result_df = pd.DataFrame()

# Calculate count for the "rule" column
result_df['rule_count'] = grouped_df['rule'].count()

# Iterate through the specified columns and calculate mean and total
for column in args.columns:
    if column in df.columns:
        result_df[f'{column}_mean'] = grouped_df[column].mean()
        result_df[f'{column}_total'] = grouped_df[column].sum()

# Reset the index to make "rule" a regular column
result_df.reset_index(inplace=True)

# Add tag column
result_df['Tag'] = args.tag


# Save the result DataFrame to a TSV file
result_df.to_csv(args.output_file, sep='\t', index=False)

print(f"Result saved to {args.output_file}")