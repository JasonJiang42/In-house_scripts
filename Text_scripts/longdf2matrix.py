import pandas as pd
import argparse

def convert_to_matrix(input_file, output_file):
    # Read the two-column file
    df = pd.read_csv(input_file, sep='\t', header=None, names=['Gene', 'Sample'])

    # Create a pivot table to form a binary matrix
    binary_matrix = df.pivot_table(index='Gene', columns='Sample', aggfunc='size', fill_value=0)

    # Convert counts to binary presence/absence
    binary_matrix = (binary_matrix > 0).astype(int)

    # Write the binary matrix to the output file
    binary_matrix.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert two-column file to binary presence/absence matrix.")
    parser.add_argument('-i', '--input', required=True, help="Input two-column file")
    parser.add_argument('-o', '--output', required=True, help="Output binary matrix file")
    
    args = parser.parse_args()
    
    convert_to_matrix(args.input, args.output)
