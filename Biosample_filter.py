import csv
import argparse


def filter_tsv(input_file, output_file):
    """
    Filters a TSV file based on the first column and specific headers.
    
    :param input_file: Path to the input TSV file.
    :param output_file: Path to save the filtered TSV output.
    """
    # Define the keywords to filter columns based on headers
    keywords = ["host", "source", "env", "collection_date", "geo_loc_name"]

    with open(input_file, "r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        
        # Identify columns to keep based on the header keywords
        selected_columns = [col for col in reader.fieldnames if any(keyword in col.lower() for keyword in keywords)]
        selected_columns.insert(0, reader.fieldnames[0])  # Always include the first column
        
        # Open the output file and write the filtered results
        with open(output_file, "w", newline="", encoding="utf-8") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=selected_columns, delimiter="\t")
            writer.writeheader()
            
            for row in reader:
                # Write only the selected columns for each row
                filtered_row = {col: row[col] for col in selected_columns}
                writer.writerow(filtered_row)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Filter a TSV file by selecting the first column and specific columns containing keywords."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input TSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file.")
    args = parser.parse_args()

    # Call the filter function
    filter_tsv(args.input, args.output)


if __name__ == "__main__":
    main()
