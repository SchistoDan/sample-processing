"""
Sample Spreadsheet Generation Script
==================================

This script processes NGS raw read files and generates samples.csv and samples_nonproject.csv files 
containing file paths and taxonomic IDs. It searches through a directory structure to find 
paired-end read files (R1 and R2) that match specific Process IDs, and associates them with 
corresponding taxonomic IDs from a metadata file. Process IDs are filtered based on project codes.

Key Features:
------------
- Recursively searches directories for paired FASTQ files (.fastq.gz)
- Matches Process IDs from sample metadata with raw read files
- Filters samples based on project codes
- Excludes 'Undetermined' and 'NC' (Negative Control) samples
- Associates taxonomic IDs with each sample
- Generates two output files:
  * samples.csv: Contains samples with valid project codes
  * samples_nonproject.csv: Contains samples without project codes

Input Requirements:
-----------------
1. Directory Structure:
   - Parent directory containing subdirectories with raw read files
   - Raw read files should be paired (R1 and R2) and in .fastq.gz format
   - File names should contain "_R1_" and "_R2_" to identify pairs

2. Sample Metadata File:
   - CSV format
   - Must contain columns: 'Process ID' and 'taxid'
   - Process IDs should match those in the raw read directory names

Output:
-------
Generates two CSV files:
1. samples.csv with columns:
   - ID: Process ID (containing a valid project code)
   - forward: Absolute path to R1 (forward) read file
   - reverse: Absolute path to R2 (reverse) read file
   - taxid: Taxonomic ID from metadata

2. samples_nonproject.csv with the same columns but containing samples
   without valid project codes

Usage:
------
python script.py [parent_directory] [sample_metadata_file]

Example:
python script.py /path/to/raw_reads /path/to/sample_metadata.csv

Dependencies:
------------
- pandas
- Python standard libraries (os, csv, sys)
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

def extract_id_reads_taxid(parent_dir, process_ids, taxid_dict):
    results = {}
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
            if r1_path and r2_path:
                results[process_id] = (r1_path, r2_path, taxid)
    return results

def write_filtered_results(results, parent_dir):
    """
    Writes raw read file paths and taxid to two separate CSV files based on project codes.
    
    Args:
        results (dict): A dictionary where keys are Process IDs and values are tuples containing:
            - Path to R1 (forward) read file
            - Path to R2 (reverse) read file
            - Taxonomic ID
        parent_dir (str): Parent directory path used to name output files
    
    Output Files:
        - samples_<parent_dir_name>.csv: Contains samples with valid project codes
        - samples_<parent_dir_name>_nonproject.csv: Contains samples without project codes
    """
    # Get directory name for file naming
    dir_name = get_directory_name(parent_dir)
    
    # Create output filenames
    project_filename = f'samples_{dir_name}.csv'
    nonproject_filename = f'samples_{dir_name}_nonproject.csv'
    
    # Separate results based on project codes
    project_samples = {pid: data for pid, data in results.items() if has_project_code(pid)}
    non_project_samples = {pid: data for pid, data in results.items() if not has_project_code(pid)}
    
    # Write project samples
    with open(project_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse', 'taxid'])
        for process_id, (r1_path, r2_path, taxid) in project_samples.items():
            writer.writerow([process_id, r1_path, r2_path, taxid])
    
    # Write non-project samples
    with open(nonproject_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse', 'taxid'])
        for process_id, (r1_path, r2_path, taxid) in non_project_samples.items():
            writer.writerow([process_id, r1_path, r2_path, taxid])
    
    return len(project_samples), len(non_project_samples), project_filename, nonproject_filename

def main(parent_dir, sample_metadata_file):
    # Read Process IDs and taxid from sample_metadata.csv
    metadata_df = pd.read_csv(sample_metadata_file, low_memory=False)  # Added low_memory=False
    process_ids = set(metadata_df['Process ID'].astype(str))
    taxid_dict = dict(zip(metadata_df['Process ID'].astype(str), metadata_df['taxid']))
    
    # Find files (forward/reverse reads)
    results = extract_id_reads_taxid(parent_dir, process_ids, taxid_dict)
    
    # Write filtered results to separate files
    project_count, non_project_count, project_file, nonproject_file = write_filtered_results(results, parent_dir)
    
    print(f"Processing complete:")
    print(f"- {project_count} samples written to '{project_file}'")
    print(f"- {non_project_count} samples written to '{nonproject_file}'")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        parent_dir = sys.argv[1]
        sample_metadata_file = sys.argv[2]
        main(parent_dir, sample_metadata_file)
    else:
        print(
        """
        Usage: python 2_sample_spreadsheet.py [parent_directory] [sample_metadata_file]
        parent_directory: parent directory of subdirectories containing raw read pairs.
        sample_metadata_file: Path to the sample_metadata.csv file.
        Example: python 2_samples_spreadsheet.py /path/to/raw_reads /path/to/sample_metadata.csv
        """
        )
