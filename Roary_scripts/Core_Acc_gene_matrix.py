## Convert the gene pres_abs matrix from roary to specific core or accessory gene matrix

import pandas as pd
import argparse

def filter_gene_matrix(gene_matrix_file, cluster_protein_file, output_file):
    # Read the cluster_protein file and extract genes
    with open(cluster_protein_file, 'r') as f:
        cluster_genes = {line.split(': ')[0] for line in f}

    # Load the gene presence/absence matrix from CSV
    gene_matrix = pd.read_csv(gene_matrix_file, index_col='Gene', low_memory=False)

    # Filter the matrix to include only genes in cluster_genes
    filtered_matrix = gene_matrix.loc[gene_matrix.index.intersection(cluster_genes)]

    # Write the filtered matrix to the output file as TSV
    filtered_matrix.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter gene matrix by cluster_protein file.")
    parser.add_argument('-i', '--input', required=True, help="Input gene matrix CSV file")
    parser.add_argument('-c', '--cluster', required=True, help="Core or accessory cluster protein file")
    parser.add_argument('-o', '--output', required=True, help="Output filtered gene matrix TSV file")
    
    args = parser.parse_args()
    
    filter_gene_matrix(args.input, args.cluster, args.output)
