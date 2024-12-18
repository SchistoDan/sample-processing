"""
Sample Metadata Generator for Taxonomic Data

This script processes taxonomic and specimen data from multiple tsv files to generate
consolidated sample metadata. It handles taxonomic classification, geographic information,
collection details, and specimen characteristics.

Key Features:
- Processes multiple tsv input files (voucher, collection, specimen, taxonomy, lab, custom fields)
- Resolves taxonomic IDs using a ranked lineage reference file
- Handles missing data with standardised 'not collected' values
- Standardises date formats
- Supports parallel processing for improved performance
- Validates taxonomic lineages for consistency

Required Input Files:
- voucher.tsv: Contains voucher specimen information
- collection.tsv: Contains collection event data
- specimen.tsv: Contains specimen-specific information
- taxonomy.tsv: Contains taxonomic classification data
- lab.tsv: Contains laboratory processing information
- merged_custom_fields.tsv: Contains additional type status information
- rankedlineage.dmp: Reference file for taxonomic lineage resolution

Output:
- sample_metadata.csv: Consolidated metadata file containing all processed information

Usage:
    python 3_local_bold_processing.py <input_directory> <rankedlineage_path> <output_directory>

Arguments:
    input_directory: Directory containing all required tsv files
    rankedlineage_path: Path to the rankedlineage.dmp file
    output_directory: Directory where the output csv will be saved

Dependencies:
    - pandas
    - numpy
    - concurrent.futures (for parallel processing)
    - datetime
    - csv
    - os
    - sys
    - time

Note:
    All missing values in any field will be replaced with 'not collected' in the output.
    Dates are standardised to YYYY-MM-DD (BCDM) format.
"""

import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import time


def load_rankedlineage(rankedlineage_path):
    """
    Load taxonomic data with comprehensive indexing at all taxonomic ranks
    and improved species name handling
    
    Parameters:
    rankedlineage_path (str): Path to the rankedlineage.dmp file
    
    Returns:
    dict: Dictionary containing indexed taxonomic data at all ranks
    """
    print(f"Loading rankedlineage file from {rankedlineage_path}...")
    
    tax_data = {
        'by_scientific_name': {},
        'by_species': {},
        'by_genus': {},
        'by_family': {},
        'by_order': {},
        'by_class': {},
        'by_phylum': {}
    }
    
    line_count = 0
    species_indexed = 0
    
    with open(rankedlineage_path, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            parts = [p.strip() for p in line.strip().split('|')]
            if len(parts) < 10:
                continue
                
            taxid = parts[0]
            scientific_name = parts[1]
            species_field = parts[2]
            genus = parts[3]
            family = parts[4]
            order = parts[5]
            class_name = parts[6]
            phylum = parts[7]
            kingdom = parts[8]
            superkingdom = parts[9]
            
            # Extract species without duplicating genus
            species = None
            if species_field:
                # Remove genus if it's at the start of species field
                if genus and species_field.lower().startswith(genus.lower()):
                    species = species_field[len(genus):].strip()
                else:
                    species = species_field
            elif genus and scientific_name.startswith(genus):
                species = scientific_name[len(genus):].strip()
                
            # Validate species name
            if species and any(x in species.lower() for x in ['sp.', 'cf.', 'aff.', 'x ', 'subsp.', 'var.']):
                species = None
            
            lineage = {
                'taxid': taxid,
                'scientific_name': scientific_name,
                'species': species,
                'genus': genus,
                'family': family,
                'order': order,
                'class': class_name,
                'phylum': phylum,
                'kingdom': kingdom,
                'superkingdom': superkingdom
            }
            
            # Index by scientific name
            tax_data['by_scientific_name'][scientific_name.lower()] = lineage
            
            # Index by species if we have it - store without duplicating genus
            if species and genus:
                species_key = f"{genus.lower()} {species.lower()}"
                tax_data['by_species'][species_key] = lineage
                species_indexed += 1
            
            # Index by genus
            if genus:
                if genus.lower() not in tax_data['by_genus']:
                    tax_data['by_genus'][genus.lower()] = []
                tax_data['by_genus'][genus.lower()].append(lineage)
            
            # Index by family
            if family:
                if family.lower() not in tax_data['by_family']:
                    tax_data['by_family'][family.lower()] = []
                tax_data['by_family'][family.lower()].append(lineage)
            
            # Index by order
            if order:
                if order.lower() not in tax_data['by_order']:
                    tax_data['by_order'][order.lower()] = []
                tax_data['by_order'][order.lower()].append(lineage)
            
            # Index by class
            if class_name:
                if class_name.lower() not in tax_data['by_class']:
                    tax_data['by_class'][class_name.lower()] = []
                tax_data['by_class'][class_name.lower()].append(lineage)
            
            # Index by phylum
            if phylum:
                if phylum.lower() not in tax_data['by_phylum']:
                    tax_data['by_phylum'][phylum.lower()] = []
                tax_data['by_phylum'][phylum.lower()].append(lineage)
            
            line_count += 1
            if line_count % 500000 == 0:
                print(f"Processed {line_count} lines... ({species_indexed} species indexed)")
                print(f"Current index sizes:")
                print(f"  Scientific names: {len(tax_data['by_scientific_name'])}")
                print(f"  Species: {len(tax_data['by_species'])}")
                print(f"  Genera: {len(tax_data['by_genus'])}")
                print(f"  Families: {len(tax_data['by_family'])}")
                print(f"  Orders: {len(tax_data['by_order'])}")
                print(f"  Classes: {len(tax_data['by_class'])}")
                print(f"  Phyla: {len(tax_data['by_phylum'])}")
    
    print(f"\nFinished loading {line_count} taxonomic records.")
    print(f"Final index sizes:")
    print(f"  Scientific names: {len(tax_data['by_scientific_name'])}")
    print(f"  Species: {len(tax_data['by_species'])}")
    print(f"  Genera: {len(tax_data['by_genus'])}")
    print(f"  Families: {len(tax_data['by_family'])}")
    print(f"  Orders: {len(tax_data['by_order'])}")
    print(f"  Classes: {len(tax_data['by_class'])}")
    print(f"  Phyla: {len(tax_data['by_phylum'])}")
    
    return tax_data

def validate_against_higher_ranks(lineage, target_ranks, validation_level='family'):
    """
    Validate taxonomy with hierarchical fallback at different taxonomic levels.
    Returns true if a match is found at the specified level, without checking lower ranks.
    
    Parameters:
    lineage (dict): The lineage to validate
    target_ranks (dict): The target taxonomy to validate against
    validation_level (str): The taxonomic level to validate at ('family', 'order', 'class', or 'phylum')
    
    Returns:
    bool: True if validation passes at the specified level
    """
    print(f"\nValidating at {validation_level} level")
    
    # Define validation hierarchy (from lowest to highest rank)
    validation_ranks = {
        'family': ['family'],
        'order': ['order'],
        'class': ['class'],
        'phylum': ['phylum']
    }
    
    ranks_to_check = validation_ranks.get(validation_level, [])
    print(f"Checking ranks: {ranks_to_check}")
    
    # Only check the specified rank, not lower ranks
    for rank in ranks_to_check:
        if target_ranks.get(rank) and lineage.get(rank):
            target_value = target_ranks[rank].lower()
            lineage_value = lineage[rank].lower()
            print(f"Comparing {rank}: {lineage_value} vs {target_value}")
            
            if lineage_value == target_value:
                print(f"âœ“ Direct match at {rank}")
                return True
            
            print(f"âœ— Mismatch at {rank}")
            return False
            
    return False

def resolve_taxid(phylum, class_name, order, family, genus, species, tax_data):
    """
    Resolve taxonomic ID with hierarchical fallback and validation.
    Tries to match at each rank level and validates against higher ranks in order.
    
    Parameters:
    phylum (str): Phylum name
    class_name (str): Class name
    order (str): Order name
    family (str): Family name
    genus (str): Genus name
    species (str): Species name
    tax_data (dict): The taxonomic data dictionary from load_rankedlineage
    
    Returns:
    tuple: (taxid, matched_rank, lineage_string, is_mismatch)
    """
    target_ranks = {
        'phylum': phylum,
        'class': class_name,
        'order': order,
        'family': family,
        'genus': genus,
        'species': species
    }
    
    validation_levels = ['family', 'order', 'class', 'phylum']
    print(f"\n==========: Looking up: {genus} {species} :==========")
    print(f"Input taxonomy: {target_ranks}")
    
    # Try species-level match first
    if species and genus:
        # Remove genus from species if it's duplicated at the start
        species_name = species
        if species.lower().startswith(genus.lower()):
            species_name = species[len(genus):].strip()
            
        species_key = f"{genus.lower()} {species_name.lower()}"
        print(f"Trying species lookup with key: {species_key}")
        if species_key in tax_data['by_species']:
            print("Found species match")
            lineage = tax_data['by_species'][species_key]
            for level in validation_levels:
                if validate_against_higher_ranks(lineage, target_ranks, level):
                    return create_return_tuple(lineage, 'species', False)
    
    # Try genus-level match
    if genus:
        print(f"Trying genus-level match for: {genus}")
        if genus.lower() in tax_data['by_genus']:
            matches = tax_data['by_genus'][genus.lower()]
            print(f"Found {len(matches)} genus matches")
            for lineage in matches:
                # For genus match, try validating at each level in order
                for level in validation_levels:
                    if validate_against_higher_ranks(lineage, target_ranks, level):
                        return create_return_tuple(lineage, 'genus', False)
        else:
            print("No genus matches found")
    
    # Try matches at each higher rank
    for rank in validation_levels:
        rank_value = target_ranks.get(rank)
        if rank_value:
            print(f"Trying {rank}-level match for: {rank_value}")
            if rank_value.lower() in tax_data[f'by_{rank}']:
                matches = tax_data[f'by_{rank}'][rank_value.lower()]
                print(f"Found {len(matches)} {rank} matches")
                # For each rank, only validate at that specific rank
                for lineage in matches:
                    if validate_against_higher_ranks(lineage, target_ranks, rank):
                        return create_return_tuple(lineage, rank, False)
            else:
                print(f"No {rank} matches found")
    
    print("No valid matches found at any rank")
    return None, 'unmatched', None, False

def create_return_tuple(lineage, matched_rank, is_mismatch):
    """Helper function to create consistent return values"""
    lineage_string = ";".join([
        f"{rank}:{lineage.get(rank, '')}" 
        for rank in ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ])
    return lineage['taxid'], matched_rank, lineage_string, is_mismatch

def locate_files_in_directory(directory):
    """
    Function to locate the required TSV files in the specified directory
    """
    voucher_file = collection_file = specimen_file = taxonomy_file = lab_file = custom_fields_file = None
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
        elif 'merged_custom_fields' in filename.lower():
            custom_fields_file = os.path.join(directory, filename)
    
    if not all([voucher_file, collection_file, specimen_file, taxonomy_file, lab_file, custom_fields_file]):
        print("Error: Some required files are missing in the directory.")
    
    return voucher_file, collection_file, specimen_file, taxonomy_file, lab_file, custom_fields_file

def load_tsv_files(tsv_directory_path):
    """
    Load all required TSV files into pandas DataFrames
    """
    print(f"Loading TSV files from {tsv_directory_path}...")
    voucher_file, collection_file, specimen_file, taxonomy_file, lab_file, custom_fields_file = locate_files_in_directory(tsv_directory_path)
    
    voucher_df = pd.read_csv(voucher_file, sep='\t', usecols=['Sample ID', 'Museum ID', 'Institution Storing'])
    collection_df = pd.read_csv(collection_file, sep='\t', usecols=['Sample ID', 'Collection Date', 'Country/Ocean', 'Collectors', 'Habitat', 'Exact Site', 'Lat', 'Lon'])
    specimen_df = pd.read_csv(specimen_file, sep='\t', usecols=['Sample ID', 'Tissue Descriptor', 'Sex', 'Life Stage'])
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t', usecols=['Sample ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Identifier'])
    lab_df = pd.read_csv(lab_file, sep='\t', usecols=['Sample ID', 'Process ID'])
    
    # Modified to skip first row and use correct column name
    custom_fields_df = pd.read_csv(custom_fields_file, 
                                  sep='\t', 
                                  skiprows=1,  # Skip the first row
                                  usecols=['SampleID', 'Type Status'])
    
    # Rename the column to match the others
    custom_fields_df = custom_fields_df.rename(columns={'SampleID': 'Sample ID'})
    
    print("TSV files loaded successfully.")
    return voucher_df, collection_df, specimen_df, taxonomy_df, lab_df, custom_fields_df

def process_sample(sample_id, voucher_df, collection_df, specimen_df, taxonomy_df, lab_df, custom_fields_df, tax_tree, fill_missing, standardise_date_format):
    """
    Process a single sample and return its metadata
    """
    # Extract taxonomy data
    phylum = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Phylum'].values[0] if not taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Phylum'].empty else np.nan
    class_name = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Class'].values[0] if not taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Class'].empty else np.nan
    order = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Order'].values[0] if not taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Order'].empty else np.nan
    family = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Family'].values[0] if not taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Family'].empty else np.nan
    genus = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Genus'].values[0] if not taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Genus'].empty else np.nan
    species = taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Species'].values[0] if not taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Species'].empty else np.nan
    
    # Fill missing values
    phylum = fill_missing(phylum)
    class_name = fill_missing(class_name)
    order = fill_missing(order)
    family = fill_missing(family)
    genus = fill_missing(genus)
    species = fill_missing(species)
    
    # Resolve taxonomic ID
    taxid, matched_rank, lineage, is_mismatch = resolve_taxid(phylum, class_name, order, family, genus, species, tax_tree)
    
    # Return processed sample data
    return {
        "Sample ID": sample_id,
        "Process ID": fill_missing(lab_df.loc[lab_df['Sample ID'] == sample_id, 'Process ID'].values[0] if not lab_df.loc[lab_df['Sample ID'] == sample_id, 'Process ID'].empty else np.nan),
        "phylum": phylum,
        "class": class_name,
        "order": order,
        "family": family,
        "genus": genus,
        "species": species,
        "taxid": fill_missing(taxid if taxid else np.nan),
        "matched_rank": matched_rank,
        "lineage": fill_missing(lineage if lineage else np.nan),
        "lineage_mismatch": "Yes" if is_mismatch else "No",
        "identified_by": fill_missing(taxonomy_df.loc[taxonomy_df['Sample ID'] == sample_id, 'Identifier'].values[0] if 'Identifier' in taxonomy_df.columns else np.nan),
        "collection_date": standardise_date_format(
            fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Collection Date'].values[0] if 'Collection Date' in collection_df.columns else np.nan)
        ),
        "geographic_location": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Country/Ocean'].values[0] if 'Country/Ocean' in collection_df.columns else np.nan),
        "geographic_location_locality": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Exact Site'].values[0] if 'Exact Site' in collection_df.columns else np.nan),
        "latitude": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Lat'].values[0] if 'Lat' in collection_df.columns else np.nan),
		"longitude": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Lon'].values[0] if 'Lon' in collection_df.columns else np.nan),
        "collected_by": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Collectors'].values[0] if 'Collectors' in collection_df.columns else np.nan),
        "habitat": fill_missing(collection_df.loc[collection_df['Sample ID'] == sample_id, 'Habitat'].values[0] if 'Habitat' in collection_df.columns else np.nan),
        "organism_part": fill_missing(specimen_df.loc[specimen_df['Sample ID'] == sample_id, 'Tissue Descriptor'].values[0] if 'Tissue Descriptor' in specimen_df.columns else np.nan),
        "sex": fill_missing(specimen_df.loc[specimen_df['Sample ID'] == sample_id, 'Sex'].values[0] if 'Sex' in specimen_df.columns else np.nan),
        "lifestage": fill_missing(specimen_df.loc[specimen_df['Sample ID'] == sample_id, 'Life Stage'].values[0] if 'Life Stage' in specimen_df.columns else np.nan),
        "specimen_voucher": fill_missing(voucher_df.loc[voucher_df['Sample ID'] == sample_id, 'Museum ID'].values[0] if 'Museum ID' in voucher_df.columns else np.nan),
        "collecting_institution": fill_missing(voucher_df.loc[voucher_df['Sample ID'] == sample_id, 'Institution Storing'].values[0] if 'Institution Storing' in voucher_df.columns else np.nan),
        "type_status": fill_missing(custom_fields_df.loc[custom_fields_df['Sample ID'] == sample_id, 'Type Status'].values[0] if not custom_fields_df.loc[custom_fields_df['Sample ID'] == sample_id, 'Type Status'].empty or custom_fields_df.loc[custom_fields_df['Sample ID'] == sample_id, 'Type Status'].values[0] == '' else np.nan),
    }

def create_sample_metadata(input_dir, rankedlineage_path, output_dir):
    print(f"Starting metadata generation from input directory: {input_dir} and rankedlineage file: {rankedlineage_path}")
    
    # Load all required files
    voucher_df, collection_df, specimen_df, taxonomy_df, lab_df, custom_fields_df = load_tsv_files(input_dir)
    tax_tree = load_rankedlineage(rankedlineage_path)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    sample_metadata_file = os.path.join(output_dir, "sample_metadata.csv")
    sample_metadata_columns = [
        "Sample ID", "Process ID", "phylum", "class", "order", "family", "genus", 
        "species", "taxid", "matched_rank", "lineage", "lineage_mismatch", "collection_date", "geographic_location",
        "geographic_location_locality", "latitude", "longitude", "collected_by",
        "habitat", "organism_part", "sex", "lifestage", "specimen_voucher",
        "collecting_institution", "identified_by", "type_status"
    ]
    
    def fill_missing(value):
        return 'not collected' if pd.isna(value) else value

    def standardise_date_format(date_str):
        if pd.isna(date_str) or date_str.strip() == '':
            return 'not collected'
        date_formats = [
            '%d-%b-%Y', '%d-%b-%y', '%d-%m-%Y', '%d-%m-%y', '%Y-%m-%d',
            '%m/%d/%Y', '%m/%d/%y', '%d/%m/%Y', '%d/%m/%y'
        ]
        date_str = date_str.strip()
        for fmt in date_formats:
            try:
                return datetime.strptime(date_str, fmt).strftime('%Y-%m-%d')
            except ValueError:
                continue
        return 'not collected'

    print(f"Writing metadata to {sample_metadata_file}...")
    total_samples = len(voucher_df['Sample ID'])
    processed_samples = 0
    start_time = time.time()

    with open(sample_metadata_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=sample_metadata_columns)
        writer.writeheader()

        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=4) as executor:
            future_to_sample = {
                executor.submit(
                    process_sample, 
                    sample_id, 
                    voucher_df, 
                    collection_df, 
                    specimen_df, 
                    taxonomy_df, 
                    lab_df, 
                    custom_fields_df,
                    tax_tree, 
                    fill_missing, 
                    standardise_date_format
                ): sample_id 
                for sample_id in voucher_df['Sample ID']
            }
            
            for future in as_completed(future_to_sample):
                sample_id = future_to_sample[future]
                try:
                    data = future.result()
                    writer.writerow(data)
                    processed_samples += 1
                    
                    # Update progress
                    progress = (processed_samples / total_samples) * 100
                    elapsed_time = time.time() - start_time
                    estimated_total_time = elapsed_time / (processed_samples / total_samples)
                    remaining_time = estimated_total_time - elapsed_time
                    
                    print(f"\rProgress: {progress:.2f}% - Processed {processed_samples}/{total_samples} samples. "
                          f"Estimated time remaining: {remaining_time:.2f} seconds", end="")
                    
                    # Flush the output to ensure it's displayed immediately
                    sys.stdout.flush()
                    
                except Exception as exc:
                    print(f"\nError processing sample {sample_id}: {exc}")

    print(f"\nSample metadata CSV file '{sample_metadata_file}' created successfully.")
    print(f"Total processing time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        input_dir = sys.argv[1]
        rankedlineage_path = sys.argv[2]
        output_dir = sys.argv[3]
        if os.path.isdir(input_dir):
            create_sample_metadata(input_dir, rankedlineage_path, output_dir)
        else:
            print(f"Error: '{input_dir}' is not a valid directory.")
    else:
        print("Usage: python 3_local_bold_processing.py <input_directory> <rankedlineage_path> <output_directory>")
