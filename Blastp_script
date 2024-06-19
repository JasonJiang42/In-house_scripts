import pandas as pd
import subprocess
import argparse
import os

def run_blastp(input_file, db, num_threads):
    blastp_output = 'blastp_results.tsv'
    blastp_command = [
        'blastp',
        '-query', input_file,
        '-db', db,
        '-out', blastp_output,
        '-outfmt', '6 qseqid sseqid length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore pident',
        '-max_target_seqs', '1',
        '-evalue', '1e-10',
        '-num_threads', str(num_threads)
    ]
    
    try:
        subprocess.run(blastp_command, check=True)
        return blastp_output
    except subprocess.CalledProcessError as e:
        print(f"Error running BLASTP: {e}")
        return None

def add_headers_calculate_coverage_and_filter(blastp_output, output_file, min_identity, min_coverage):
    # Define the headers
    headers = ["qseqid", "sseqid", "length", "mismatch", "gapopen", 
               "qlen", "qstart", "qend", "slen", "sstart", "send", 
               "evalue", "bitscore", "pident"]
    
    try:
        # Read the BLASTP output TSV file without headers
        df = pd.read_csv(blastp_output, sep='\t', header=None)
        
        # Assign the headers to the DataFrame
        df.columns = headers
        
        # Calculate the coverage as a percentage
        df['coverage'] = (df['length'] / df['slen']) * 100
        
        # Filter the DataFrame based on min_identity and min_coverage
        filtered_df = df[(df['pident'] >= min_identity) & (df['coverage'] >= min_coverage)]

        # Write the filtered DataFrame to the output TSV file
        filtered_df.to_csv(output_file, sep='\t', index=False)
        
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run BLASTP, assign headers to the output, calculate coverage, and filter based on identity and coverage.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input FASTA file.")
    parser.add_argument('-db', '--database', required=True, help="Path to the BLASTP database.")
    parser.add_argument('-o', '--output', required=True, help="Path to the output TSV file.")
    parser.add_argument('--minid', type=float, required=True, help="Minimum identity percentage to filter the results.")
    parser.add_argument('--mincov', type=float, required=True, help="Minimum coverage percentage to filter the results.")
    parser.add_argument('-t', '--num_threads', type=int, default=1, help="Number of threads to use for BLASTP.")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Run BLASTP
    blastp_output = run_blastp(args.input, args.database, args.num_threads)
    
    # If BLASTP ran successfully, process the output
    if blastp_output:
        add_headers_calculate_coverage_and_filter(blastp_output, args.output, args.minid, args.mincov)
        
        # Optionally, remove the temporary BLASTP output file
        os.remove(blastp_output)
