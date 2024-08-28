import sys
import os
import pandas as pd
import csv
from datetime import datetime
from Bio import Entrez
from pathlib import Path
import time




def usage():
    print("""
    Usage: python sample_processing.py path/to/BOLD/download/dir path/to/raw/read/dir path/to/output/dir/[output].csv
        path/to/BOLD/download/dir: Directory containing BOLD tsv files (voucher.tsv, collection_data.tsv, specimen_details.tsv, taxonomy.tsv, and lab.tsv).
        path/to/raw/read/dir: Parent directory with subdirs containing raw PE reads files.
        path/to/output/dir/[output].csv: Directory to create [output].csv, samples.csv and unique_taxids.txt within. Need to name [output].csv, whereas samples.csv and unique_taxids.txt are named as is.
    Example:
        python sample_processing.py /path/to/tsv_directory /path/to/raw_reads /path/to/output_file.csv
    """)




#Load selected fields from TSV files
def load_tsv_files(tsv_directory_path):
    voucher_df = pd.read_csv(os.path.join(tsv_directory_path, 'voucher.tsv'), sep='\t', usecols=['Sample ID', 'Museum ID', 'Institution Storing'])

    collection_df = pd.read_csv(os.path.join(tsv_directory_path, 'collection_data.tsv'), sep='\t', usecols=['Sample ID', 'Collection Date', 'Country/Ocean', 'Collectors', 'Habitat', 'Exact Site', 'Lat', 'Lon'])

    specimen_df = pd.read_csv(os.path.join(tsv_directory_path, 'specimen_details.tsv'), sep='\t', usecols=['Sample ID', 'Tissue Descriptor', 'Sex', 'Life Stage'])

    taxonomy_df = pd.read_csv(os.path.join(tsv_directory_path, 'taxonomy.tsv'), sep='\t', usecols=['Sample ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Identifier'])

    lab_df = pd.read_csv(os.path.join(tsv_directory_path, 'lab.tsv'), sep='\t', usecols=['Sample ID', 'Process ID'])
    
    return voucher_df, collection_df, specimen_df, taxonomy_df, lab_df




#Fetch taxonomic ids for each sample using taxonomic ranks
def fetch_taxid(taxon_name, taxonomic_rank, retries=4, delay=5):
    attempt = 0
    while attempt < retries:
        try:
            term = f"{taxon_name}[{taxonomic_rank}]"
            handle = Entrez.esearch(db="taxonomy", term=term)
            record = Entrez.read(handle)
            handle.close()
            if record['IdList']:
                return record['IdList'][0]
            else:
                return None
        except Exception as e:
            print(f"Error fetching TaxID for {taxon_name} ({taxonomic_rank}) on attempt {attempt + 1}: {e}")
            attempt += 1
            time.sleep(delay)
            if 'Search Backend failed' in str(e):
                break
    
    print(f"Failed to fetch TaxID for {taxon_name} ({taxonomic_rank}) after {retries} attempts.")
    return None





#convert collection dates to BCDM format
def standardize_date_format(date_str):
    if pd.isna(date_str) or date_str.strip() == '':
        return 'not collected'
    
    date_formats = [
        '%d-%b-%Y',   #e.g 16-Jun-1927
        '%d-%b-%y',   #e.g 16-Jun-27
        '%d-%m-%Y',   #e.g 16-06-1927
        '%d-%m-%y',   #e.g 16-06-27
        '%Y-%m-%d',   #e.g 1927-06-16
        '%m/%d/%Y',   #e.g 06/16/1927
        '%m/%d/%y',   #e.g 06/16/27
        '%d/%m/%Y',   #e.g 16/06/1927
        '%d/%m/%y'    #e.g 16/06/27
    ]
    
    date_str = date_str.strip()
    
    for fmt in date_formats:
        try:
            return datetime.strptime(date_str, fmt).strftime('%Y-%m-%d')
        except ValueError:
            continue
    
    return 'not collected'




#Format metadata as needed and successuively search from low to high taxonomic ranks for valid taxid
def format_metadata(row, unique_taxids):
    sample_id = row.get('Sample ID')
    order = row.get('Order')
    family = row.get('Family')
    genus = row.get('Genus')
    species = row.get('Species')

    taxid = None
    matched_rank = 'not collected'

    if species:
        species_taxid = fetch_taxid(species, 'species')
        if species_taxid:
            taxid = species_taxid
            matched_rank = "species"
            unique_taxids.add(species_taxid)

    if not taxid and genus:
        genus_taxid = fetch_taxid(genus, 'genus')
        if genus_taxid:
            taxid = genus_taxid
            matched_rank = "genus"
            unique_taxids.add(genus_taxid)

    if not taxid and family:
        family_taxid = fetch_taxid(family, 'family')
        if family_taxid:
            taxid = family_taxid
            matched_rank = "family"
            unique_taxids.add(family_taxid)

    if not taxid and order:
        order_taxid = fetch_taxid(order, 'order')
        if order_taxid:
            taxid = order_taxid
            matched_rank = "order"
            unique_taxids.add(order_taxid)

    collection_date = standardize_date_format(row.get('Collection Date'))

#Handle empty lat and lon fields
    def parse_float(value):
        try:
            return float(value)
        except (ValueError, TypeError):
            return None

    latitude = row.get('Lat')
    longitude = row.get('Lon')

    latitude = 'not collected' if pd.isna(latitude) or parse_float(latitude) is None else latitude
    longitude = 'not collected' if pd.isna(longitude) or parse_float(longitude) is None else longitude


#Handle missing or invalid values and populate with 'not colelcted'
    def handle_missing_or_empty(value):
        if pd.isna(value) or value.strip() == '':
            return 'not collected'
        value_str = str(value).strip().lower()
        if value_str in ['unrecoverable', 'not collected']:
            return 'not collected'
        return value

    geographic_location = handle_missing_or_empty(row.get('Country/Ocean'))

    geographic_location_locality = handle_missing_or_empty(row.get('Exact Site'))

    habitat = handle_missing_or_empty(row.get('Habitat'))

    organism_part = handle_missing_or_empty(row.get('Tissue Descriptor'))

    sex = handle_missing_or_empty(row.get('Sex'))

    lifestage = handle_missing_or_empty(row.get('Life Stage'))

    return {
        'Sample ID': sample_id,
        'Process ID': row.get('Process ID'),
        'identified_by': row.get('Identifier', 'not collected'),
        'order': order,
        'family': family,
        'genus': genus,
        'species': species,
        'taxid': taxid,
        'matched_rank': matched_rank,
        'collection date': collection_date,
        'geographic_location': geographic_location,
        'geographic_location_locality': geographic_location_locality,
        'latitude': latitude,
        'longitude': longitude,
        'collected_by': row.get('Collectors', 'not collected'),
        'habitat': habitat,
        'organism part': organism_part,
        'sex': sex,
        'lifestage': lifestage,
        'specimen_voucher': row.get('Museum ID', 'not collected'),
        'collecting institution': row.get('Institution Storing', 'not collected')
    }




#Merge data between input.tsv files based on Sample ID
def merge_data(voucher_df, collection_df, specimen_df, taxonomy_df, lab_df):
    unique_taxids = set()

    merged_df = pd.merge(lab_df[['Sample ID', 'Process ID']], taxonomy_df[['Sample ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Identifier']], on='Sample ID', how='left')
    merged_df = pd.merge(merged_df, collection_df, on='Sample ID', how='left')
    merged_df = pd.merge(merged_df, specimen_df, on='Sample ID', how='left')
    merged_df = pd.merge(merged_df, voucher_df, on='Sample ID', how='left')

    string_columns = ['Habitat', 'Sex', 'Life Stage', 'Tissue Descriptor']
    for col in string_columns:
        if col in merged_df.columns:
            merged_df[col] = merged_df[col].astype('object')

    merged_df = merged_df.apply(lambda row: pd.Series(format_metadata(row, unique_taxids)), axis=1)

    for col in merged_df.columns:
        if merged_df[col].dtype == 'object':
            merged_df[col].fillna('not collected', inplace=True)

    return merged_df, unique_taxids




#Create [output].csv
def write_to_csv(df, output_file):
    df.to_csv(output_file, index=False)
    print(f"Final output written to {output_file}")



#Locate R1 and R2 raw read files in subdirs of given parent dir
def find_files(folder_path):
    results = {}
    for root, dirs, files in os.walk(folder_path):
        for subfolder in dirs:
            r1_path, r2_path = None, None
            subfolder_path = os.path.join(root, subfolder)
            for filename in os.listdir(subfolder_path):
                if "_R1_" in filename:
                    r1_path = os.path.abspath(os.path.join(subfolder_path, filename))
                elif "_R2_" in filename:
                    r2_path = os.path.abspath(os.path.join(subfolder_path, filename))
            if r1_path and r2_path:
                subfolder_id = '-'.join(subfolder.split('-')[2:])
                results[subfolder_id] = (r1_path, r2_path)
    return results




#Create samples.csv with relevant fields
def write_samples_csv(df, read_paths, output_dir):
    samples_data = []
    for _, row in df.iterrows():
        process_id = row.get('Process ID')
        taxid = row.get('taxid')
        r1_path, r2_path = read_paths.get(process_id, (None, None))
        samples_data.append([process_id, r1_path, r2_path, taxid])
    
    samples_csv_path = os.path.join(output_dir, 'samples.csv')
    with open(samples_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Process ID', 'forward read', 'reverse read', 'taxid'])
        writer.writerows(samples_data)
    
    print(f"Samples CSV file written to {samples_csv_path}")




#Create unique_taxids.txt contianing dedupped taxids
def write_unique_taxids_txt(df, output_dir):
    unique_taxids = set(df['taxid'].dropna())
    unique_taxids_txt_path = os.path.join(output_dir, 'unique_taxids.txt')
    
    with open(unique_taxids_txt_path, 'w') as file:
        for taxid in sorted(unique_taxids):
            file.write(f"{taxid}\n")
    
    print(f"Unique taxids text file written to {unique_taxids_txt_path}")




#MAIN
def main(tsv_directory_path, raw_reads_directory_path, output_file):
    Entrez.email = "d.parsons@nhm.ac.uk"
    
#Load input.tsv files, merge, write to [output].csv, samples.csv, and unique_taxids.csv
    voucher_df, collection_df, specimen_df, taxonomy_df, lab_df = load_tsv_files(tsv_directory_path)

    merged_df, unique_taxids = merge_data(voucher_df, collection_df, specimen_df, taxonomy_df, lab_df)
    
    write_to_csv(merged_df, output_file)

    read_paths = find_files(raw_reads_directory_path)
    
    output_dir = os.path.dirname(output_file)

    write_samples_csv(merged_df, read_paths, output_dir)

    write_unique_taxids_txt(merged_df, output_dir)

    print("Processing complete.")




if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Error: Incorrect number of arguments.")
        usage()
        sys.exit(1)
    
    tsv_directory_path = sys.argv[1]
    raw_reads_directory_path = sys.argv[2]
    output_file = sys.argv[3]

    if not os.path.isdir(tsv_directory_path):
        print(f"Error: {tsv_directory_path} is not a valid directory.")
        sys.exit(1)
    
    if not os.path.isdir(raw_reads_directory_path):
        print(f"Error: {raw_reads_directory_path} is not a valid directory.")
        sys.exit(1)

    main(tsv_directory_path, raw_reads_directory_path, output_file)
