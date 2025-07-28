from Bio import SeqIO
import argparse


def extract_region_from_genbank(input_file, output_file, start, end):
    """
    Extract a region from a GenBank file based on start and end positions.

    :param input_file: Path to the input GenBank file
    :param output_file: Path to the output GenBank file
    :param start: Start position (1-based index)
    :param end: End position (1-based index)
    """
    # Parse the GenBank file using BioPython
    with open(input_file, "r") as infile:
        records = list(SeqIO.parse(infile, "genbank"))

    # Extract the region from each record in the file
    extracted_records = []
    for record in records:
        # Extract the subsequence based on the start and end positions
        sub_record = record[start - 1 : end]  # 1-based to 0-based indexing

        # Add missing molecule_type annotation if not present
        if "molecule_type" not in sub_record.annotations:
            sub_record.annotations["molecule_type"] = "DNA"

        # Update metadata for the extracted region
        sub_record.id = f"{record.id}_region_{start}_{end}"
        sub_record.description = f"Extracted region {start}-{end} from {record.id}"
        extracted_records.append(sub_record)

    # Write the extracted region(s) to the output file
    with open(output_file, "w") as outfile:
        SeqIO.write(extracted_records, outfile, "genbank")


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Extract a specific region from a GenBank file."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input GenBank file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output GenBank file.")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start position (1-based).")
    parser.add_argument("-e", "--end", type=int, required=True, help="End position (1-based).")
    args = parser.parse_args()

    # Extract the region from the GenBank file
    extract_region_from_genbank(args.input, args.output, args.start, args.end)


if __name__ == "__main__":
    main()
