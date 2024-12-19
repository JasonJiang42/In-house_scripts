from Bio import Entrez
import csv
import argparse
import time


def fetch_project_title_and_description(accession_list, email, delay=0.5):
    """
    Fetch Project_Title and Project_Description for a list of BioProject accessions.

    :param accession_list: List of BioProject accession numbers
    :param email: Email address for NCBI Entrez
    :param delay: Delay (in seconds) between API requests to avoid being blocked
    :return: List of dictionaries with BioProject, Project_Title, and Project_Description
    """
    Entrez.email = email  # Set the email address for NCBI Entrez
    results = []

    for accession in accession_list:
        print(f"Processing BioProject: {accession}")
        try:
            # Use esearch to get the UID for the BioProject accession
            handle = Entrez.esearch(db="bioproject", term=accession)
            record = Entrez.read(handle)
            handle.close()

            # If no UID is found, log a warning and skip
            if not record.get("IdList"):
                print(f"Warning: No UID found for BioProject {accession}")
                results.append({"BioProject": accession, "Project_Title": "Not Found", "Project_Description": "Not Found"})
                continue

            uid = record["IdList"][0]

            # Use esummary to fetch details for the UID
            handle = Entrez.esummary(db="bioproject", id=uid)
            summary = Entrez.read(handle)
            handle.close()

            # Parse the Entrez output and extract Project_Title and Project_Description
            document = summary.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            if document:
                project = document[0]
                title = project.get("Project_Title", "Not Found")
                description = project.get("Project_Description", "Not Found")
            else:
                title = "Not Found"
                description = "Not Found"

            # Save results
            results.append({"BioProject": accession, "Project_Title": title, "Project_Description": description})
            print(f"Accession: {accession}, Title: {title}, Description: {description}")

        except Exception as e:
            print(f"Error processing BioProject {accession}: {e}")
            results.append({"BioProject": accession, "Project_Title": "Error", "Project_Description": "Error"})

        # Add a delay between requests to avoid hitting NCBI rate limits
        time.sleep(delay)

    return results


def read_accession_file(input_file):
    """
    Read a list of BioProject accessions from a file.

    :param input_file: Path to the input file containing BioProject accessions
    :return: List of BioProject accessions
    """
    with open(input_file, "r") as file:
        return [line.strip() for line in file if line.strip()]


def write_output_file(output_file, results):
    """
    Write the results to a TSV file.

    :param output_file: Path to the output file
    :param results: List of dictionaries with BioProject, Project_Title, and Project_Description
    """
    with open(output_file, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=["BioProject", "Project_Title", "Project_Description"], delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Fetch Project_Title and Project_Description for a list of BioProject accessions from NCBI."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input file containing BioProject accessions.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file.")
    parser.add_argument(
        "-e", "--email", required=True, help="Email address to use for NCBI Entrez (required for making queries)."
    )
    args = parser.parse_args()

    # Read the input file
    accession_list = read_accession_file(args.input)

    # Fetch project details
    results = fetch_project_title_and_description(accession_list, args.email)

    # Write the results to the output file
    write_output_file(args.output, results)


if __name__ == "__main__":
    main()
