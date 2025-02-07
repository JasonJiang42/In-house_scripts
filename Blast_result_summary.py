import pandas as pd
import argparse

def process_blastn_results(input_file, output_file):
    """
    Process a BLASTn result file, reformat it, and summarize the output with one subject per query.
    """
    # Define the column names of the input BLASTn file
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp"
    ]

    # Read the BLASTn result file into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t", header=None, names=columns)

    # Calculate Total_score as the sum of bitscores for each query-subject pair
    df["Total_score"] = df.groupby(["qseqid", "sseqid"])["bitscore"].transform("sum")

    # Group by query (`qseqid`) and select the best hit for each query based on the highest Total_score
    best_hits = df.sort_values("Total_score", ascending=False).groupby("qseqid").first().reset_index()

    # Rename and select the required columns for output
    best_hits = best_hits.rename(columns={
        "qseqid": "query",
        "sseqid": "subject",
        "bitscore": "Max_score",
        "pident": "Per.identity",
        "evalue": "E.value",
        "qcovs": "qcov"
    })

    # Rearrange columns in the desired order, excluding mapped_length
    best_hits = best_hits[[
        "subject", "query", "Max_score", "Total_score", "qcov",
        "Per.identity", "E.value"
    ]]

    # Save the reformatted and summarized output to the specified file
    best_hits.to_csv(output_file, sep="\t", index=False)
    print(f"Output written to {output_file}")

# Main function to handle command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat and summarize BLASTn results.")
    parser.add_argument("-i", "--input", required=True, metavar="INPUT_FILE",
                        help="Input BLASTn result file.")
    parser.add_argument("-o", "--output", required=True, metavar="OUTPUT_FILE",
                        help="Output file for reformatted and summarized results.")
    args = parser.parse_args()

    # Process the BLASTn results
    process_blastn_results(args.input, args.output)
