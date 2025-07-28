import subprocess
import argparse
import os
from Bio import SeqIO
import pandas as pd


def get_contig_lengths(query_fasta):
    """
    Extract the lengths of all contigs from the query FASTA file.

    Parameters:
        query_fasta (str): Path to the query FASTA file.

    Returns:
        dict: A dictionary where keys are contig headers (query names) and values are their lengths.
    """
    contig_lengths = {}
    for record in SeqIO.parse(query_fasta, "fasta"):
        contig_lengths[record.id] = len(record.seq)
    return contig_lengths


def add_contig_length_and_calculate_coverage(coords_file, query_fasta, coverage_threshold, output_tab_file, output_fasta_file):
    """
    Add query contig lengths, calculate coverage, filter based on coverage, and output results.

    Parameters:
        coords_file (str): Path to the coords.tab file.
        query_fasta (str): Path to the query FASTA file.
        coverage_threshold (float): Minimum coverage threshold for filtering contigs.
        output_tab_file (str): Path to the output tab file with contig lengths and coverage.
        output_fasta_file (str): Path to the output FASTA file for filtered contigs.
    """
    # Load contig lengths from the query FASTA file
    contig_lengths = get_contig_lengths(query_fasta)

    # Read the coords.tab file into a DataFrame
    coords_df = pd.read_csv(coords_file, sep="\t", header=None)

    # Assuming coords.tab has 9 columns, add column headers
    coords_df.columns = [
        "Start1", "End1", "Start2", "End2", "Len1", "Len2", "Identity", "Ref", "Query"
    ]

    # Map the contig lengths to the "Query" column (col9)
    coords_df["Query_Length"] = coords_df["Query"].map(contig_lengths)

    # Calculate coverage for each row: Len2 / Query_Length
    coords_df["Coverage"] = coords_df["Len2"] / coords_df["Query_Length"]

    # Group by Query to calculate total mapped length and cap it at Query_Length
    grouped = coords_df.groupby("Query").agg(
        Total_Mapped_Length=("Len2", "sum"),
        Query_Length=("Query_Length", "first"),  # Query length is the same for all rows in the group
        Ref=("Ref", "first")  # Get the first reference name (assuming one reference per contig)
    )

    # Cap the total mapped length at the query length
    grouped["Capped_Mapped_Length"] = grouped[["Total_Mapped_Length", "Query_Length"]].min(axis=1)

    # Calculate total coverage
    grouped["Total_Coverage"] = grouped["Capped_Mapped_Length"] / grouped["Query_Length"]

    # Save the detailed tab file with lengths, coverage, and reference
    grouped.reset_index().to_csv(output_tab_file, sep="\t", index=False)
    print(f"Detailed tab file saved to: {output_tab_file}")

    # Filter contigs based on the coverage threshold
    filtered_contigs = grouped[grouped["Total_Coverage"] >= coverage_threshold]

    # Extract the filtered contig names
    filtered_contig_names = filtered_contigs.index.tolist()

    # Write the filtered contigs to the output FASTA file
    with open(output_fasta_file, "w") as output_handle:
        for record in SeqIO.parse(query_fasta, "fasta"):
            if record.id in filtered_contig_names:
                SeqIO.write(record, output_handle, "fasta")

    print(f"Filtered contigs saved to: {output_fasta_file}")


def run_nucmer(reference, query, output_prefix, coverage_threshold, output_tab_file, output_fasta_file):
    """
    Run nucmer on the given reference and query, process the output, and filter contigs by coverage.

    Parameters:
        reference (str): File path to the reference input.
        query (str): File path to the query input.
        output_prefix (str): Prefix for nucmer output files.
        coverage_threshold (float): Minimum coverage threshold for filtering contigs.
        output_tab_file (str): Path to the output tab file with contig lengths and coverage.
        output_fasta_file (str): Path to the output FASTA file for filtered contigs.
    """
    try:
        # Step 1: Run nucmer
        nucmer_command = ["nucmer", reference, query, "-p", output_prefix]
        print(f"Running: {' '.join(nucmer_command)}")
        subprocess.run(nucmer_command, check=True)

        # Step 2: Run show-coords to generate coords.tab
        delta_file = f"{output_prefix}.delta"
        coords_file = "coords.tab"
        show_coords_command = ["show-coords", "-T", "-H", delta_file]
        print(f"Running: {' '.join(show_coords_command)} > {coords_file}")
        with open(coords_file, "w") as coords_output:
            subprocess.run(show_coords_command, stdout=coords_output, check=True)

        # Step 3: Add query contig lengths, calculate coverage, and filter contigs
        add_contig_length_and_calculate_coverage(
            coords_file, query, coverage_threshold, output_tab_file, output_fasta_file
        )

        # Step 4: Clean up temporary files
        print(f"Cleaning up temporary files: {delta_file} and {coords_file}")
        os.remove(delta_file)
        os.remove(coords_file)

        print("Process completed successfully.")

    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
    except FileNotFoundError as e:
        print(f"Command not found: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run nucmer, calculate coverage, and filter query contigs by coverage."
    )
    parser.add_argument(
        "-r", "--reference", required=True, help="Path to the reference input file."
    )
    parser.add_argument(
        "-q", "--query", required=True, help="Path to the query input file."
    )
    parser.add_argument(
        "-p", "--prefix", default="nucmer_output", help="Prefix for nucmer output files (default: nucmer_output)."
    )
    parser.add_argument(
        "-c", "--coverage", type=float, required=True, help="Minimum coverage threshold for filtering contigs."
    )
    parser.add_argument(
        "-t", "--tab_output", required=True, help="Path to the output tab file with contig lengths, coverage, and reference."
    )
    parser.add_argument(
        "-o", "--fasta_output", required=True, help="Path to the output FASTA file for filtered contigs."
    )
    args = parser.parse_args()

    # Run nucmer workflow
    run_nucmer(args.reference, args.query, args.prefix, args.coverage, args.tab_output, args.fasta_output)


if __name__ == "__main__":
    main()
