"""
Sample Metadata Generator for Taxonomic Data

This script processes taxonomic and specimen data from multiple TSV files to generate
consolidated sample metadata. It handles taxonomic classification, geographic information,
collection details, and specimen characteristics.

Key Features:
- Processes multiple TSV input files (voucher, collection, specimen, taxonomy, lab, custom fields)
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
    input_directory: Directory containing all required TSV files
    rankedlineage_path: Path to the rankedlineage.dmp file
    output_directory: Directory where the output CSV will be saved

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


def determine_rank(parts):
    """
    Determine the rank of a tax_name based on which column first contains data
    parts[1] = tax_name
    parts[2] = species
    parts[3] = genus
    parts[4] = family
    parts[5] = order
    parts[6] = class
    parts[7] = phylum
    parts[8] = kingdom
    parts[9] = superkingdom
    """
    if parts[2].strip():  # If species column has data
        return 'subspecies'
    
    rank_map = {
        3: 'species',
        4: 'genus',
        5: 'family',
        6: 'order',
        7: 'class',
        8: 'phylum',
        9: 'kingdom',
        10: 'superkingdom'
    }
    
    for idx, value in enumerate(parts[3:], start=3):
        if value.strip():
            return rank_map.get(idx, 'unknown')
    
    return 'unknown'

def validate_higher_ranks(lineage, target_ranks):
    """
    Validate if lineage matches higher taxonomic ranks
    Returns True if all available ranks match
    """
    rank_order = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    
    for rank in rank_order:
        if lineage.get(rank) and target_ranks.get(rank):
            if lineage[rank].lower() != target_ranks[rank].lower():
                return False
    return True

def find_best_match(matches, target_ranks):
    """
    From multiple matches, find the one that best matches higher ranks
    """
    best_match = None
    best_match_score = -1
    
    rank_order = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    
    for match in matches:
        score = 0
        for rank in rank_order:
            if match.get(rank) and target_ranks.get(rank):
                if match[rank].lower() == target_ranks[rank].lower():
                    score += 1
                else:
                    # Penalise mismatches at lower ranks more heavily
                    score -= (len(rank_order) - rank_order.index(rank))
        
        if score > best_match_score:
            best_match_score = score
            best_match = match
    
    return best_match if best_match_score >= 0 else None

def search_by_rank(search_term, rank, tax_tree, target_ranks):
    """
    Search for a term at a specific rank, validating against higher ranks
    Returns (taxid, matched_rank, lineage) or (None, None, None)
    """
    matches = []
    
    # Search in the specific rank category if it exists
    if rank in tax_tree:
        for tax_name, lineage in tax_tree[rank].items():
            if tax_name.lower() == search_term.lower():
                matches.append(lineage)
    
    # If we found matches, validate against higher ranks
    if matches:
        best_match = find_best_match(matches, target_ranks)
        if best_match:
            return best_match['taxid'], rank, best_match
    
    return None, None, None

def fetch_taxid_from_local(genus, species, tax_tree, target_ranks):
    """
    Attempt to find taxid with hierarchical fallback and validation
    """
    # First try species-level match if species is provided
    if species:
        # Try exact species match under the specified genus
        if genus in tax_tree:
            genus_data = tax_tree[genus]
            species_matches = []
            
            # Check direct species entries
            for tax_name, lineage in genus_data.items():
                if isinstance(lineage, dict) and not 'subspecies' in lineage:
                    if tax_name.lower() == species.lower():
                        species_matches.append(lineage)
            
            # Check subspecies entries
            for species_entry in genus_data.values():
                if isinstance(species_entry, dict) and 'subspecies' in species_entry:
                    for subspecies_data in species_entry['subspecies'].values():
                        if subspecies_data['species'].lower() == species.lower():
                            species_matches.append(subspecies_data)
            
            if species_matches:
                best_match = find_best_match(species_matches, target_ranks)
                if best_match:
                    return best_match['taxid'], 'species', best_match
    
    # Try searching at each rank level with validation
    search_ranks = [
        ('genus', genus),
        ('family', target_ranks.get('family')),
        ('order', target_ranks.get('order')),
        ('class', target_ranks.get('class')),
        ('phylum', target_ranks.get('phylum'))
    ]
    
    for rank, search_term in search_ranks:
        if search_term:  # Only search if we have a term to search for
            taxid, matched_rank, lineage = search_by_rank(search_term, rank, tax_tree, target_ranks)
            if taxid:
                return taxid, matched_rank, lineage
    
    return None, 'unmatched', None

def resolve_taxid(phylum, class_name, order, family, genus, species, tax_tree):
    """
    Resolve taxonomic ID with hierarchical fallback and lineage retrieval
    """
    target_ranks = {
        'phylum': phylum,
        'class': class_name,
        'order': order,
        'family': family,
        'genus': genus,
        'species': species
    }
    
    taxid, matched_rank, lineage = fetch_taxid_from_local(genus, species, tax_tree, target_ranks)
    
    if not taxid:
        return None, 'unmatched', None, False
    
    # Create lineage string
    lineage_string = ";".join([f"{rank}:{lineage.get(rank, '')}" for rank in ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']])
    
    # Verify the match is consistent with provided higher taxonomy
    if lineage and not validate_higher_ranks(lineage, target_ranks):
        return taxid, matched_rank, lineage_string, True
    
    return taxid, matched_rank, lineage_string, False

def load_rankedlineage(rankedlineage_path):
    print(f"Loading rankedlineage file from {rankedlineage_path}...")
    tax_tree = {}
    line_count = 0
    
    with open(rankedlineage_path, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            parts = [p.strip() for p in line.strip().split('|')]
            taxid = parts[0]
            tax_name = parts[1]
            rank = determine_rank(parts)
            
            lineage = {
                'taxid': taxid,
                'tax_name': tax_name,
                'rank': rank,
                'species': parts[2] if parts[2].strip() else None,
                'genus': parts[3] if parts[3].strip() else None,
                'family': parts[4] if parts[4].strip() else None,
                'order': parts[5] if parts[5].strip() else None,
                'class': parts[6] if parts[6].strip() else None,
                'phylum': parts[7] if parts[7].strip() else None,
                'kingdom': parts[8] if parts[8].strip() else None,
                'superkingdom': parts[9] if parts[9].strip() else None
            }
            
            # Build tree structure based on available data
            if rank == 'subspecies':
                genus = parts[3]
                species_name = parts[2]
                tax_tree.setdefault(genus, {}).setdefault(species_name, {})['subspecies'] = {
                    tax_name: lineage
                }
            elif rank == 'species':
                genus = parts[3]
                tax_tree.setdefault(genus, {})[tax_name] = lineage
            else:
                tax_tree.setdefault(rank, {})[tax_name] = lineage
            
            line_count += 1
            if line_count % 500000 == 0:
                print(f"Processed {line_count} lines from rankedlineage.dmp...")
    
    print(f"Finished loading rankedlineage.dmp ({line_count} lines processed).")
    return tax_tree

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
