import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def extract_region_from_genbank(input_file, start, end, output_file):
    """
    Extracts a region from a GenBank file based on the start and end positions,
    while preserving annotations and adjusting feature locations.

    Args:
        input_file (str): Path to the input GenBank file.
        start (int): Start position (1-based).
        end (int): End position (1-based).
        output_file (str): Path to the output GenBank file.
    """
    with open(input_file, "r") as infile:
        # Parse the GenBank file
        record = SeqIO.read(infile, "genbank")
        
        # Extract the sequence region
        extracted_sequence = record.seq[start - 1:end]  # Convert to 0-based indexing
        
        # Update the sequence in the record
        record.seq = extracted_sequence

        # Adjust features to fit within the extracted region
        updated_features = []
        for feature in record.features:
            # Get the start and end positions of the feature
            feature_start = int(feature.location.start)
            feature_end = int(feature.location.end)
            
            # Check if the feature overlaps with the extracted region
            if feature_end >= start and feature_start <= end:
                # Adjust the feature's location relative to the extracted region
                new_start = max(0, feature_start - start)
                new_end = min(end - start, feature_end - start)

                # Create a new feature with the adjusted location
                new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)
                new_feature = SeqFeature(
                    location=new_location,
                    type=feature.type,
                    qualifiers=feature.qualifiers
                )
                updated_features.append(new_feature)

        # Update the record's features
        record.features = updated_features

        # Correct the source feature if it exists
        for feature in record.features:
            if feature.type == "source":
                feature.location = FeatureLocation(0, len(extracted_sequence))  # Entire extracted region

        # Write the updated record to the output file
        with open(output_file, "w") as outfile:
            SeqIO.write(record, outfile, "genbank")
    
    print(f"Extracted region from {start} to {end} saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a region from a GenBank file.")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank file")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start position (1-based)")
    parser.add_argument("-e", "--end", type=int, required=True, help="End position (1-based)")
    parser.add_argument("-o", "--output", required=True, help="Output GenBank file for the extracted region")
    
    args = parser.parse_args()
    
    # Call the function with parsed arguments
    extract_region_from_genbank(args.input, args.start, args.end, args.output)