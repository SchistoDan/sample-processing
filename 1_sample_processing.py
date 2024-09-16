import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime
from Bio import Entrez
import time
import random
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed

def locate_files_in_directory(directory):
    """
    Locates specific TSV files (voucher, collection, specimen, taxonomy, lab) in the given directory.
    
    Args:
        directory (str): Path to the directory containing the TSV files.

    Returns:
        tuple: Paths to the voucher file, collection file, specimen file, taxonomy file, and lab file (if found).
    """
    voucher_file = collection_file = specimen_file = taxonomy_file = lab_file = None
    for filename in os.listdir(directory):
        if 'voucher' in filename.lower():
            voucher_file = os.path.join(directory, filename)
        elif 'collection' in filename.lower():
            collection_file = os.path.join(directory, filename)
        elif 'specimen' in filename.lower():
            specimen_file = os.path.join(directory, filename)
        elif 'taxonomy' in filename.lower():
            taxonomy_file = os.path.join(directory, filename)
        elif 'lab' in filename.lower():
            lab_file = os.path.join(directory, filename)
    
    return voucher_file, collection_file, specimen_file, taxonomy_file, lab_file

def load_tsv_files(tsv_directory_path):
    """
    Loads TSV files using pandas and selects the necessary columns.

    Args:
        tsv_directory_path (str): Directory containing the TSV files.

    Returns:
        tuple: DataFrames for voucher, collection, specimen, taxonomy, and lab data.
    """
    voucher_df = pd.read_csv(os.path.join(tsv_directory_path, 'voucher.tsv'), sep='\t', usecols=['Sample ID', 'Museum ID', 'Institution Storing'])
    collection_df = pd.read_csv(os.path.join(tsv_directory_path, 'collection_data.tsv'), sep='\t', usecols=['Sample ID', 'Collection Date', 'Country/Ocean', 'Collectors', 'Habitat', 'Exact Site', 'Lat', 'Lon'])
    specimen_df = pd.read_csv(os.path.join(tsv_directory_path, 'specimen_details.tsv'), sep='\t', usecols=['Sample ID', 'Tissue Descriptor', 'Sex', 'Life Stage'])
    taxonomy_df = pd.read_csv(os.path.join(tsv_directory_path, 'taxonomy.tsv'), sep='\t', usecols=['Sample ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Identifier'])
    lab_df = pd.read_csv(os.path.join(tsv_directory_path, 'lab.tsv'), sep='\t', usecols=['Sample ID', 'Process ID'])
    
    return voucher_df, collection_df, specimen_df, taxonomy_df, lab_df

def fetch_taxid(genus, species, process_id=None):
    """
    Fetches the taxonomic ID (taxid) for a given genus and species using the NCBI Entrez API.

    Args:
        genus (str): The genus name of the organism.
        species (str): The species name of the organism (optional).
    Returns:
        tuple: A tuple containing the taxid (if found) and the taxonomic rank (e.g., 'species' or 'genus').
    """
    max_retries = 5  # Maximum number of retry attempts
    retry_delay = 2  # Initial delay between retries in seconds

    for attempt in range(1, max_retries + 1):
        try:
            search_term = f"{genus} {species}" if species else genus
            handle = Entrez.esearch(db="taxonomy", term=search_term, retmode="xml")
            record = Entrez.read(handle)
            handle.close()

            if record['IdList']:
                return record['IdList'][0], 'species' if species else 'genus'
            return None, 'unmatched'

        except Exception as e:
            if "HTTP 429" in str(e) or "HTTP 500" in str(e):
                if attempt < max_retries:
                    delay = retry_delay * (2 ** (attempt - 1)) + random.uniform(0, 1)
                    print(f"Error fetching taxid for {genus} {species} (Process ID: {process_id}) on attempt {attempt}. "
                          f"Retrying in {delay:.2f} seconds. Error: {e}")
                    time.sleep(delay)
                else:
                    print(f"Error fetching taxid for {genus} {species} (Process ID: {process_id}) after {max_retries} attempts. "
                          f"Giving up. Error: {e}")
                    return None, 'error'
            else:
                print(f"Error fetching taxid for {genus} {species} (Process ID: {process_id}): {e}")
                return None, 'error'

def resolve_taxid(order, family, genus, species, process_id):
    """
    Resolves the taxid based on the provided order, family, genus, and species.
    Starts with the highest level (order) and progresses down to species, replacing any
    higher-level matches with more specific ones.

    Args:
        order (str): The order of the organism.
        family (str): The family of the organism.
        genus (str): The genus of the organism.
        species (str): The species of the organism.
        process_id (str): The Process ID for the current sample (for debugging purposes).

    Returns:
        tuple: A tuple containing the most specific resolved taxid and the corresponding matched taxonomic rank.
    """
    taxid = None
    matched_rank = 'unmatched'

    # Ensure all inputs are strings
    order = str(order) if not pd.isna(order) else 'not collected'
    family = str(family) if not pd.isna(family) else 'not collected'
    genus = str(genus) if not pd.isna(genus) else 'not collected'
    species = str(species) if not pd.isna(species) else 'not collected'

    if order and order.lower() != 'not collected':
        taxid, matched_rank = fetch_taxid(order, None, process_id)
        if taxid:
            matched_rank = 'order'

    if family and family.lower() != 'not collected':
        family_taxid, family_rank = fetch_taxid(family, None, process_id)
        if family_taxid:
            taxid, matched_rank = family_taxid, 'family'

    if genus and genus.lower() != 'not collected':
        genus_taxid, genus_rank = fetch_taxid(genus, None, process_id)
        if genus_taxid:
            taxid, matched_rank = genus_taxid, 'genus'

    if genus and species and genus.lower() != 'not collected' and species.lower() != 'not collected':
        species_taxid, species_rank = fetch_taxid(genus, species, process_id)
        if species_taxid:
            taxid, matched_rank = species_taxid, 'species'

    return taxid, matched_rank

def standardise_date_format(date_str):
    """
    standardises date format to YYYY-MM-DD. If the date is missing or cannot be parsed, returns 'not collected'.

    Args:
        date_str (str): The date string to standardise.

    Returns:
        str: The standardised date string in YYYY-MM-DD format or 'not collected'.
    """
    if pd.isna(date_str) or date_str.strip() == '':
        return 'not collected'
    
    date_formats = [
        '%d-%b-%Y',   # e.g 16-Jun-1927
        '%d-%b-%y',   # e.g 16-Jun-27
        '%d-%m-%Y',   # e.g 16-06-1927
        '%d-%m-%y',   # e.g 16-06-27
        '%Y-%m-%d',   # e.g 1927-06-16
        '%m/%d/%Y',   # e.g 06/16/1927
        '%m/%d/%y',   # e.g 06/16/27
        '%d/%m/%Y',   # e.g 16/06/1927
        '%d/%m/%y'    # e.g 16/06/27
    ]
    
    date_str = date_str.strip()
    
    for fmt in date_formats:
        try:
            return datetime.strptime(date_str, fmt).strftime('%Y-%m-%d')
        except ValueError:
            continue
    
    return 'not collected'

def create_sample_metadata(input_dir):
    voucher_file, collection_file, specimen_file, taxonomy_file, lab_file = locate_files_in_directory(input_dir)
    print(f"Located files: voucher_file={voucher_file}, collection_file={collection_file}, specimen_file={specimen_file}, taxonomy_file={taxonomy_file}, lab_file={lab_file}")

    voucher_df, collection_df, specimen_df, taxonomy_df, lab_df = load_tsv_files(input_dir)

    print("Voucher Sample IDs:", voucher_df['Sample ID'].head().tolist())  
    print("Collection Sample IDs:", collection_df['Sample ID'].head().tolist())  
    print("Specimen Sample IDs:", specimen_df['Sample ID'].head().tolist())  
    print("Taxonomy Sample IDs:", taxonomy_df['Sample ID'].head().tolist())  
    print("Lab Sample IDs:", lab_df['Sample ID'].head().tolist())  

    sample_metadata_file = "sample_metadata.csv"
    sample_metadata_columns = [
        "Sample ID", "Process ID", "identified_by", "order", "family", "genus", "species", "taxid",
        "matched_rank", "collection date", "geographic_location", "geographic_location_locality", "latitude",
        "longitude", "collected_by", "habitat", "organism part", "sex", "lifestage", "specimen_voucher", 
        "collecting institution"
    ]

    def fill_missing(value):
        return 'not collected' if pd.isna(value) else value

    with open(sample_metadata_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=sample_metadata_columns)
        writer.writeheader()

        for sample_id in voucher_df['Sample ID']:
            order = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Order'].values
            family = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Family'].values
            genus = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Genus'].values
            species = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Species'].values
            process_id = lab_df.loc[lab_df['Sample ID'] == sample_id, 'Process ID'].values

            # Extract values or default to 'not collected'
            order = fill_missing(order[0] if len(order) > 0 else np.nan)
            family = fill_missing(family[0] if len(family) > 0 else np.nan)
            genus = fill_missing(genus[0] if len(genus) > 0 else np.nan)
            species = fill_missing(species[0] if len(species) > 0 else np.nan)
            process_id = fill_missing(process_id[0] if len(process_id) > 0 else np.nan)

            taxid, matched_rank = resolve_taxid(order, family, genus, species, process_id)

            row = {
                "Sample ID": sample_id,
                "Process ID": process_id,
                "identified_by": fill_missing(taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Identifier'].values[0] if 'Identifier' in taxonomy_df.columns else np.nan),
                "order": order,
                "family": family,
                "genus": genus,
                "species": species,
                "taxid": fill_missing(taxid if taxid else np.nan),
                "matched_rank": matched_rank,
                "collection date": standardise_date_format(
                    fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Collection Date'].values[0] if 'Collection Date' in collection_df.columns else np.nan)
                ),
                "geographic_location": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Country/Ocean'].values[0] if 'Country/Ocean' in collection_df.columns else np.nan),
                "geographic_location_locality": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Exact Site'].values[0] if 'Exact Site' in collection_df.columns else np.nan),
                "latitude": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Lat'].values[0] if 'Lat' in collection_df.columns else np.nan),
                "longitude": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Lon'].values[0] if 'Lon' in collection_df.columns else np.nan),
                "collected_by": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Collectors'].values[0] if 'Collectors' in collection_df.columns else np.nan),
                "habitat": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Habitat'].values[0] if 'Habitat' in collection_df.columns else np.nan),
                "organism part": fill_missing(specimen_df.loc[specimen_df['Sample ID'] == sample_id, 'Tissue Descriptor'].values[0] if 'Tissue Descriptor' in specimen_df.columns else np.nan),
                "sex": fill_missing(specimen_df.loc[specimen_df['Sample ID'] == sample_id, 'Sex'].values[0] if 'Sex' in specimen_df.columns else np.nan),
                "lifestage": fill_missing(specimen_df.loc[specimen_df['Sample ID'] == sample_id, 'Life Stage'].values[0] if 'Life Stage' in specimen_df.columns else np.nan),
                "specimen_voucher": fill_missing(voucher_df.loc[voucher_df['Sample ID'] == sample_id, 'Museum ID'].values[0] if 'Museum ID' in voucher_df.columns else np.nan),
                "collecting institution": fill_missing(voucher_df.loc[voucher_df['Sample ID'] == sample_id, 'Institution Storing'].values[0] if 'Institution Storing' in voucher_df.columns else np.nan),
            }
            writer.writerow(row)

    print(f"Sample metadata CSV file '{sample_metadata_file}' created successfully.")

def print_usage():
    """
    Prints usage information for the script.
    """
    usage_info = """
Usage: python metadata_script.py [input_directory]

[input_directory]: Directory containing the .tsv files required by the script.

The script expects the following TSV files in the specified directory:
- A file containing voucher data with 'Sample ID' as a key column.
- A file containing collection data with 'Sample ID' as a key column.
- A file containing specimen details with 'Sample ID' as a key column.
- A file containing taxonomy data with 'Sample ID' as a key column.
- A file containing lab data with 'Sample ID' as a key column.

Ensure all required TSV files are present and named appropriately.
"""
    print(usage_info)

if __name__ == "__main__":
    Entrez.email = "###&###.ac.uk"
    Entrez.api_key = "###"

    if len(sys.argv) == 2:
        input_dir = sys.argv[1]
        if os.path.isdir(input_dir):
            create_sample_metadata(input_dir)
        else:
            print(f"Error: '{input_dir}' is not a valid directory.")
            print_usage()
    else:
        print_usage()
