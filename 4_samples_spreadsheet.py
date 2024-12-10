"""
Sample Spreadsheet Generation Script
==================================

This script processes NGS raw read files and generates samples.csv, samples_nonproject.csv, 
and samples_types.csv files containing file paths and taxonomic IDs. It searches through a 
directory structure to find paired-end read files (R1 and R2) that match specific Process IDs, 
and associates them with corresponding taxonomic IDs and type status from a metadata file.

Key Features:
------------
- Recursively searches directories for paired FASTQ files (.fastq.gz)
- Matches Process IDs from sample metadata with raw read files
- Filters samples based on project codes and type status
- Excludes 'Undetermined' and 'NC' (Negative Control) samples
- Associates taxonomic IDs with each sample
- Generates three output files:
  * samples.csv: Contains samples with valid project codes
  * samples_nonproject.csv: Contains samples without project codes
  * samples_types.csv: Contains samples marked as types in metadata

Input Requirements:
-----------------
1. Directory Structure:
   - Parent directory containing subdirectories with raw read files
   - Raw read files should be paired (R1 and R2) and in .fastq.gz format
   - File names should contain "_R1_" and "_R2_" to identify pairs

2. Sample Metadata File:
   - CSV format
   - Must contain columns: 'Process ID', 'taxid', and 'type_status'
   - Process IDs should match those in the raw read directory names
"""

import os
import csv
import sys
import pandas as pd

# Define valid project codes
PROJECT_CODES = {
    'BSNHM', 'NHMCG', 'BGETR', 'BSUIO', 'BGLIB', 'BSNTN',
    'BGEGR', 'DTAUT', 'HMAUT', 'BGENL', 'DTKNU', 'BBIOP',
    'BHNHM', 'UNIFI', 'DTULO', 'MEAMP', 'MUSBA', 'BGSNL',
    'BGSNH', 'BGEPL', 'EBGEP', 'BSCRO', 'BIOSC', 'INVBG',
    'BCEMI', 'ILECA', 'ALPFU'
}

def get_directory_name(path):
    # Remove trailing slash if present
    path = path.rstrip(os.path.sep)
    # Get the last part of the path
    return os.path.basename(path)

def has_project_code(process_id):
    return any(code in process_id for code in PROJECT_CODES)

def extract_id_reads_taxid_type(parent_dir, metadata_df):
    """
    Extract file paths and metadata for each sample.
    
    Args:
        parent_dir (str): Parent directory containing sample subdirectories
        metadata_df (pd.DataFrame): DataFrame containing sample metadata
    
    Returns:
        dict: Dictionary with Process IDs as keys and tuples of 
             (r1_path, r2_path, taxid, type_status) as values
    """
    results = {}
    process_ids = set(metadata_df['Process ID'].astype(str))
    
    # Create dictionaries for quick lookup
    taxid_dict = dict(zip(metadata_df['Process ID'].astype(str), metadata_df['taxid']))
    type_status_dict = dict(zip(metadata_df['Process ID'].astype(str), metadata_df['type_status']))
    
    for root, dirs, files in os.walk(parent_dir):
        for subdir in dirs:
            if "Undetermined" in subdir or "NC" in subdir:
                continue
                
            subdir_path = os.path.join(root, subdir)
            process_id = None
            for pid in process_ids:
                if pid in subdir:
                    process_id = pid
                    break
            if not process_id:
                continue
                
            r1_path, r2_path = None, None
            for filename in os.listdir(subdir_path):
                if "_R1_" in filename and filename.endswith(".fastq.gz"):
                    r1_path = os.path.abspath(os.path.join(subdir_path, filename))
                elif "_R2_" in filename and filename.endswith(".fastq.gz"):
                    r2_path = os.path.abspath(os.path.join(subdir_path, filename))
                    
            taxid = taxid_dict.get(process_id, 'NA')
            type_status = type_status_dict.get(process_id, '')
            
            if r1_path and r2_path:
                results[process_id] = (r1_path, r2_path, taxid, type_status)
    
    return results

def write_filtered_results(results, parent_dir):
    """
    Writes raw read file paths and metadata to three separate CSV files.
    
    Args:
        results (dict): Dictionary with Process IDs as keys and tuples of 
                       (r1_path, r2_path, taxid, type_status) as values
        parent_dir (str): Parent directory path used to name output files
    
    Output Files:
        - samples_<parent_dir_name>.csv: Contains samples with valid project codes
        - samples_<parent_dir_name>_nonproject.csv: Contains samples without project codes
        - samples_<parent_dir_name>_types.csv: Contains samples marked as types
    """
    dir_name = get_directory_name(parent_dir)
    
    # Create output filenames
    project_filename = f'samples_{dir_name}.csv'
    nonproject_filename = f'samples_{dir_name}_nonproject.csv'
    types_filename = f'samples_{dir_name}_types.csv'
    
    # Separate results based on criteria
    project_samples = {}
    non_project_samples = {}
    type_samples = {}
    
    for pid, (r1_path, r2_path, taxid, type_status) in results.items():
        # Check if sample is a type specimen
        if isinstance(type_status, str) and 'type' in type_status.lower():
            type_samples[pid] = (r1_path, r2_path, taxid, type_status)
        
        # Separate by project code
        if has_project_code(pid):
            project_samples[pid] = (r1_path, r2_path, taxid, type_status)
        else:
            non_project_samples[pid] = (r1_path, r2_path, taxid, type_status)
    
    # Write project samples
    with open(project_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse', 'taxid', 'type_status'])
        for process_id, (r1_path, r2_path, taxid, type_status) in project_samples.items():
            writer.writerow([process_id, r1_path, r2_path, taxid, type_status])
    
    # Write non-project samples
    with open(nonproject_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse', 'taxid', 'type_status'])
        for process_id, (r1_path, r2_path, taxid, type_status) in non_project_samples.items():
            writer.writerow([process_id, r1_path, r2_path, taxid, type_status])
    
    # Write type samples
    with open(types_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse', 'taxid', 'type_status'])
        for process_id, (r1_path, r2_path, taxid, type_status) in type_samples.items():
            writer.writerow([process_id, r1_path, r2_path, taxid, type_status])
    
    return (len(project_samples), len(non_project_samples), len(type_samples),
            project_filename, nonproject_filename, types_filename)

def main(parent_dir, sample_metadata_file):
    # Read metadata file
    metadata_df = pd.read_csv(sample_metadata_file, low_memory=False)
    
    # Find files and extract metadata
    results = extract_id_reads_taxid_type(parent_dir, metadata_df)
    
    # Write filtered results to separate files
    (project_count, non_project_count, type_count,
     project_file, nonproject_file, types_file) = write_filtered_results(results, parent_dir)
    
    print(f"Processing complete:")
    print(f"- {project_count} samples written to '{project_file}'")
    print(f"- {non_project_count} samples written to '{nonproject_file}'")
    print(f"- {type_count} type specimens written to '{types_file}'")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        parent_dir = sys.argv[1]
        sample_metadata_file = sys.argv[2]
        main(parent_dir, sample_metadata_file)
    else:
        print(
        """
        Usage: python script.py [parent_directory] [sample_metadata_file]
        parent_directory: parent directory of subdirectories containing raw read pairs.
        sample_metadata_file: Path to the sample_metadata.csv file.
        Example: python script.py /path/to/raw_reads /path/to/sample_metadata.csv
        """
        )
