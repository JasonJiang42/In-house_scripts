import pandas as pd
import argparse

def process_blastn_results(input_file, output_file):
    """
    Process a BLASTn result file, summarize the output per subject and per query,
    and aggregate statistics like total score, max score, and average percent identity.
    """
    # Define the column names of the input BLASTn file
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp"
    ]

    # Read the BLASTn result file into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t", header=None, names=columns)

    # Group by both query (`qseqid`) and subject (`sseqid`) and calculate summary statistics
    summary = df.groupby(["qseqid", "sseqid"]).agg(
        Total_score=("bitscore", "sum"),         # Sum of all scores for query-subject pair
        Max_score=("bitscore", "max"),          # Maximum score for query-subject pair
        Per_identity=("pident", "mean"),        # Average percent identity
        E_value=("evalue", "min"),              # Minimum e-value
        Avg_qcov=("qcovs", "mean")              # Average query coverage
    ).reset_index()

    # Rename columns for output consistency
    summary = summary.rename(columns={
        "qseqid": "query",
        "sseqid": "subject",
        "Per_identity": "Per.identity",
        "E_value": "E.value",
        "Avg_qcov": "qcov"
    })

    # Save the summarized output to the specified file
    summary.to_csv(output_file, sep="\t", index=False)
    print(f"Summarized output written to {output_file}")

# Main function to handle command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize BLASTn results per subject and per query.")
    parser.add_argument("-i", "--input", required=True, metavar="INPUT_FILE",
                        help="Input BLASTn result file.")
    parser.add_argument("-o", "--output", required=True, metavar="OUTPUT_FILE",
                        help="Output file for summarized results.")
    args = parser.parse_args()

    # Process the BLASTn results
    process_blastn_results(args.input, args.output)
