import argparse
import pandas as pd

def generate_sankey_links(input_file, output_file, cols):
    # Load the TSV file with header
    data = pd.read_csv(input_file, sep="\t", header=0)
    # Select only the specified columns
    selected = data[cols]
    # Generate all adjacent source-target pairs
    links_list = []
    for i in range(len(cols) - 1):
        pair = selected[[cols[i], cols[i+1]]].rename(columns={cols[i]: 'source', cols[i+1]: 'target'})
        links_list.append(pair)
    links = pd.concat(links_list)
    # Group by source and target to calculate the value (counts)
    links = links.groupby(['source', 'target']).size().reset_index(name='value')
    # Save the result as a CSV file
    links.to_csv(output_file, index=False)
    print(f"Sankey link file saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Sankey link file from a TSV file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input TSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file.")
    parser.add_argument("--cols", required=True, help="Comma-separated list of column names to use (e.g., A,B,C,D,E)")
    args = parser.parse_args()
    cols = [col.strip() for col in args.cols.split(",")]
    if len(cols) < 2:
        raise ValueError("You must specify at least two column names with --cols")
    generate_sankey_links(args.input, args.output, cols)