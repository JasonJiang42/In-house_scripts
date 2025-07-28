import argparse
import pandas as pd
import numpy as np
import os

def extract_and_create_heatmap(input_file, reference_length, output_file):
    segments = np.zeros(reference_length // 100)
    
    with open(input_file) as infile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            try:
                start_pos = int(parts[0])
                end_pos = int(parts[1])
                identity = float(parts[6])
                
                start_segment = (start_pos - 1) // 100
                end_segment = (end_pos - 1) // 100
                
                for segment in range(start_segment, end_segment + 1):
                    if segment < len(segments):
                        segments[segment] = max(segments[segment], identity)
            except ValueError:
                continue
    
    row_name = os.path.basename(input_file)
    heatmap_df = pd.DataFrame([segments], index=[row_name])
    heatmap_df.to_csv(output_file, sep='\t', header=True, index=True)
    print(heatmap_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract columns and convert to a heatmap matrix.")
    parser.add_argument('-i', '--input', required=True, help="Input file with nucmer coords data")
    parser.add_argument('-r', '--reference_length', type=int, required=True, help="Length of the reference sequence")
    parser.add_argument('-o', '--output', required=True, help="Output file for heatmap matrix")
    
    args = parser.parse_args()
    
    extract_and_create_heatmap(args.input, args.reference_length, args.output)
