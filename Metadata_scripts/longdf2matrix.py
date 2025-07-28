import pandas as pd
import argparse

def convert_to_count_matrix(input_file, output_file):
    # Read the two-column file
    df = pd.read_csv(input_file, sep='\t', header=None, names=['Gene', 'Sample'])

    # Create a pivot table to count occurrences
    count_matrix = df.pivot_table(index='Gene', columns='Sample', aggfunc='size', fill_value=0)

    # Write the count matrix to the output file
    count_matrix.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert two-column file to gene presence count matrix.")
    parser.add_argument('-i', '--input', required=True, help="Input two-column file")
    parser.add_argument('-o', '--output', required=True, help="Output count matrix file")
    
    args = parser.parse_args()
    
    convert_to_count_matrix(args.input, args.output)
