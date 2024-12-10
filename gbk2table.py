from Bio import SeqIO
import argparse

def genbank_to_table(input_file, output_file):
    # Open the output file
    with open(output_file, 'w') as out:
        # Write the header line
        out.write("CDS_start\tCDS_end\tgene\tproduct\tlocus_tag\n")
        
        # Parse the GenBank file
        for record in SeqIO.parse(input_file, "genbank"):
            for feature in record.features:
                # Process only CDS features
                if feature.type == "CDS":
                    # Extract start and end positions using int(location) to avoid the warning
                    start = int(feature.location.start) + 1  # Convert to 1-based position
                    end = int(feature.location.end)
                    
                    # Extract gene, product, and locus_tag information (use 'NA' if not available)
                    gene = feature.qualifiers.get('gene', ['NA'])[0]
                    product = feature.qualifiers.get('product', ['NA'])[0]
                    locus_tag = feature.qualifiers.get('locus_tag', ['NA'])[0]
                    
                    # Write the information to the output file
                    out.write(f"{start}\t{end}\t{gene}\t{product}\t{locus_tag}\n")

    print(f"Conversion complete! Table written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GenBank file to table format.")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank file.")
    parser.add_argument("-o", "--output", required=True, help="Output table file.")
    
    args = parser.parse_args()
    genbank_to_table(args.input, args.output)
