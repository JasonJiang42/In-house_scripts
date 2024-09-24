## python script to convert the two list dataframe into group matrix for figure plot

import pandas as pd
import argparse

def group_numbers(input_file, output_file):
    # Load the data
    df = pd.read_csv(input_file, sep='\t', header=None, names=['Group', 'Number'])

    # Group by 'Group' and aggregate 'Number' into lists
    grouped = df.groupby('Group')['Number'].apply(list).reset_index()

    # Find the maximum number of numbers in any group
    max_length = max(grouped['Number'].apply(len))

    # Create a DataFrame for transposed data
    transposed_data = {row['Group']: row['Number'] + [''] * (max_length - len(row['Number'])) for _, row in grouped.iterrows()}

    # Convert to DataFrame and transpose it
    transposed_df = pd.DataFrame(transposed_data)

    # Write the transposed output to a TSV file
    transposed_df.to_csv(output_file, sep='\t', index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Group numbers by their group and transpose output.")
    parser.add_argument('-i', '--input', required=True, help="Input TSV file")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")
    
    args = parser.parse_args()
    
    group_numbers(args.input, args.output)
