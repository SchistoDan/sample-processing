#!/usr/bin/env python3

"""
DNA Sequencing Sample Processing Script:

This script processes raw DNA sequencing read files, and matches them with sample metadata, 
supporting flexible directory structures and file naming patterns. 
It features progress tracking and detailed logging of the processing steps.

Key Features:
-------------
1. Flexible Directory Structure Support:
   - Can process files directly in the input directory
   - Can recursively search subdirectories
   - Automatically detects whether to use subdirectory searching based on input directory content

2. Process ID Detection:
   - Identifies process IDs in both filenames and directory paths
   - Position-agnostic matching (process ID can be anywhere in the string)
   - Supports various delimiter patterns (_, -, /)
   - Handles compound IDs (e.g., XE-4013-BGSNL100-23)

3. File Pattern Support:
   - Handles multiple fastq file extensions (.fastq, .fq, with/without .gz)
   - Case-insensitive extension matching
   - Flexible R1/R2 read pair identification
   - Automatically excludes 'Undetermined' and 'NC' (negative control) files

4. Progress Tracking and Logging:
   - Displays progress bars for file processing
   - Detailed logging of file detection and matching
   - Logs both to console and file for debugging

Output Files:
-------------
The script generates three categorised output files:
1. samples_<dirname>.csv: Samples with valid project codes
2. samples_<dirname>_nonproject.csv: Samples without project codes
3. samples_<dirname>_types.csv: Samples marked as types

Each output file contains:
- ID: Process ID of the sample
- forward: Absolute path to R1 (forward) read file
- reverse: Absolute path to R2 (reverse) read file
- taxid: Taxonomic ID from metadata
- type_status: Type specimen status from metadata

Usage:
------
python 4_sample_sheet.py [parent_directory] [sample_metadata_file]

Arguments:
    parent_directory: Directory containing FASTQ/FQ files (directly (flat structure) or in subdirectories (nested structure))
    sample_metadata_file: CSV file with sample metadata

Required Metadata CSV Columns:
    - Process ID: Unique identifier for each sample
    - taxid: Taxonomic identifier
    - type_status: Type specimen status

Example Directory Structures Supported:
1. Flat structure:
   /parent_dir/
   â”œâ”€â”€ BSNHM593-24_R1.fq
   â”œâ”€â”€ BSNHM593-24_R2.fq
   â””â”€â”€ ...

2. Nested structure:
   /parent_dir/XE-4013/
   â””â”€â”€ 20240906_LH00179_0123_A22CKGHLT4/
       â”œâ”€â”€ Sample_XE-4013-BGSNL096-23/
       â”‚   â”œâ”€â”€ BGSNL096-23_R1.fastq.gz
       â”‚   â””â”€â”€ BGSNL096-23_R2.fastq.gz
       â””â”€â”€ ...

Dependencies:
    - Python 3.x
    - pandas
    - tqdm (for progress bars)
    - Standard library modules: os, csv, sys, logging
"""



import os
import csv
import sys
import pandas as pd
import logging
from tqdm import tqdm


# Define BGE project codes
PROJECT_CODES = {
    'BSNHM', 'NHMCG', 'BGETR', 'BSUIO', 'BGLIB', 'BSNTN',
    'BGEGR', 'DTAUT', 'HMAUT', 'BGENL', 'DTKNU', 'BBIOP',
    'BHNHM', 'UNIFI', 'DTULO', 'MEAMP', 'MUSBA', 'BGSNL',
    'BGSNH', 'BGEPL', 'EBGEP', 'BSCRO', 'BIOSC', 'INVBG',
    'BCEMI', 'ILECA', 'ALPFU'
}


def setup_logging(log_filename):
    """Set up logging to both file and console"""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(logging.DEBUG)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


def get_directory_name(path):
    path = path.rstrip(os.path.sep)
    return os.path.basename(path)


def is_fastq_file(filename):
    """Check if a filename corresponds to a fastq file."""
    fastq_extensions = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
    filename_lower = filename.lower()
    return any(filename_lower.endswith(ext) for ext in fastq_extensions)


def find_process_id_in_string(string, process_ids, logger=None):
    """
    Find a process ID in a string.
    Handles IDs in filenames, directory names, and compound IDs.
    """
    if logger:
        logger.debug("Trying to find process ID in string: %s", string)
    
    # Clean up the string & try each process ID
    base_string = os.path.splitext(os.path.splitext(string)[0])[0]  
    
    for pid in process_ids:
        pid_str = str(pid)  # Ensure process ID is a string
        
        if not pid_str.strip():
            continue
            
        if pid_str in base_string:
            # Get the index where the process ID was found
            idx = base_string.index(pid_str)
            
            before = base_string[idx-1] if idx > 0 else None
            after = base_string[idx+len(pid_str)] if idx+len(pid_str) < len(base_string) else None
            
            # Process ID is valid if it's bounded by separators or string boundaries
            if (before is None or before in ['_', '-', '/', ' ']) and \
               (after is None or after in ['_', '-', '/', ' ']):
                if logger:
                    logger.debug(f"Found process ID: {pid_str}")
                    logger.debug(f"Position: {idx}")
                    logger.debug(f"Before: {before}, After: {after}")
                return pid_str
    
    if logger:
        logger.debug("No process ID found")
    return None


def extract_id_reads_taxid_type(parent_dir, metadata_df, logger):
    results = {}
    process_ids = set(metadata_df['Process ID'].astype(str))
    
    logger.debug("\nDEBUG: Number of Process IDs in metadata: %d", len(process_ids))
    logger.debug("DEBUG: First few Process IDs from metadata: %s", list(process_ids)[:5])
    
    taxid_dict = dict(zip(metadata_df['Process ID'].astype(str), metadata_df['taxid']))
    type_status_dict = dict(zip(metadata_df['Process ID'].astype(str), metadata_df['type_status']))
    
    # First check if there are fastq files in the input directory
    input_dir_files = os.listdir(parent_dir)
    fastq_files_in_input = [f for f in input_dir_files if is_fastq_file(f)]
    using_subdirs = len(fastq_files_in_input) == 0
    
    logger.debug("\nDEBUG: FASTQ files found in input directory: %d", len(fastq_files_in_input))
    logger.debug("DEBUG: Will search in subdirectories: %s", using_subdirs)
    
    if not using_subdirs:
        logger.debug("\nDEBUG: Processing files in input directory")
        r1_files = {}
        r2_files = {}
        
        for filename in tqdm(input_dir_files, desc="Processing input files"):
            if not is_fastq_file(filename) or "Undetermined" in filename or "NC" in filename:
                continue
            
            logger.debug("\nDEBUG: Processing file: %s", filename)
            process_id = find_process_id_in_string(filename, process_ids, logger)
            
            if not process_id:
                logger.debug("DEBUG: No process ID found in filename: %s", filename)
                continue
            
            logger.debug("DEBUG: Found process ID: %s", process_id)
                
            filepath = os.path.abspath(os.path.join(parent_dir, filename))
            logger.debug("File path: %s", filepath)
            if "R1" in filename:
                logger.debug("Found R1 file - adding to r1_files with process_id: %s", process_id)
                r1_files[process_id] = filepath
            elif "R2" in filename:
                logger.debug("Found R2 file - adding to r2_files with process_id: %s", process_id)
                r2_files[process_id] = filepath
            else:
                logger.debug("No R1/R2 pattern found in filename: %s", filename)
    else:
        logger.debug("\nDEBUG: Searching in subdirectories")
        r1_files = {}
        r2_files = {}
        
        # First, count total files for progress bar
        total_files = sum([len([f for f in files if is_fastq_file(f)])
                         for _, _, files in os.walk(parent_dir)])
        
        processed_files = 0
        pbar = tqdm(total=total_files, desc="Processing FASTQ files")
        
        for root, _, files in os.walk(parent_dir):
            fastq_files = [f for f in files if is_fastq_file(f)]
            if not fastq_files:
                continue
                
            logger.debug("\nDEBUG: Found %d FASTQ files in %s", len(fastq_files), root)
            
            for filename in fastq_files:
                if "Undetermined" in filename or "NC" in filename:
                    continue
                
                # Try to find process ID in both filename and directory path
                process_id = find_process_id_in_string(filename, process_ids)
                if not process_id:
                    dir_path = os.path.relpath(root, parent_dir)
                    process_id = find_process_id_in_string(dir_path, process_ids)
                
                if process_id:
                    logger.debug("DEBUG: Found process ID %s in %s", 
                               process_id, os.path.join(root, filename))
                    
                    filepath = os.path.abspath(os.path.join(root, filename))
                    if "_R1_" in filename:
                        r1_files[process_id] = filepath
                    elif "_R2_" in filename:
                        r2_files[process_id] = filepath
                
                processed_files += 1
                pbar.update(1)
        
        pbar.close()
    
    logger.debug("\nDEBUG: Found R1 files for Process IDs: %s", list(r1_files.keys()))
    logger.debug("DEBUG: Number of R1 files found: %d", len(r1_files))
    logger.debug("DEBUG: Found R2 files for Process IDs: %s", list(r2_files.keys()))
    logger.debug("DEBUG: Number of R2 files found: %d", len(r2_files))
    
    matched_pairs = set(r1_files.keys()) & set(r2_files.keys())
    logger.debug("\nDEBUG: Number of matched pairs found: %d", len(matched_pairs))
    if len(matched_pairs) > 0:
        logger.debug("DEBUG: First few matched pairs: %s", list(matched_pairs)[:5])
    
    for process_id in set(r1_files.keys()) & set(r2_files.keys()):
        results[process_id] = (r1_files[process_id], r2_files[process_id],
                             taxid_dict.get(str(process_id), 'NA'),
                             type_status_dict.get(str(process_id), ''))
    
    logger.debug("\nDEBUG: Final number of matched pairs: %d", len(results))
    return results


def write_filtered_results(results, parent_dir, logger):
    dir_name = get_directory_name(parent_dir)
    
    project_filename = f'samples_{dir_name}.csv'
    nonproject_filename = f'samples_{dir_name}_nonproject.csv'
    types_filename = f'samples_{dir_name}_types.csv'
    
    project_samples = {}
    non_project_samples = {}
    type_samples = {}
    
    for pid, (r1_path, r2_path, taxid, type_status) in results.items():
        if isinstance(type_status, str) and 'type' in type_status.lower():
            type_samples[pid] = (r1_path, r2_path, taxid, type_status)
        
        if any(code in pid for code in PROJECT_CODES):
            project_samples[pid] = (r1_path, r2_path, taxid, type_status)
        else:
            non_project_samples[pid] = (r1_path, r2_path, taxid, type_status)
    
    for filename, samples in [(project_filename, project_samples),
                            (nonproject_filename, non_project_samples),
                            (types_filename, type_samples)]:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['ID', 'forward', 'reverse', 'taxid', 'type_status'])
            for process_id, (r1_path, r2_path, taxid, type_status) in samples.items():
                writer.writerow([process_id, r1_path, r2_path, taxid, type_status])
    
    return (len(project_samples), len(non_project_samples), len(type_samples),
            project_filename, nonproject_filename, types_filename)


def main(parent_dir, sample_metadata_file):
    # Set up logging
    dir_name = get_directory_name(parent_dir)
    log_filename = f'samples_{dir_name}.log'
    logger = setup_logging(log_filename)
    
    logger.debug("DEBUG: Reading metadata file: %s", sample_metadata_file)
    metadata_df = pd.read_csv(sample_metadata_file, low_memory=False)
    logger.debug("DEBUG: Metadata shape: %s", metadata_df.shape)
    
    results = extract_id_reads_taxid_type(parent_dir, metadata_df, logger)
    
    (project_count, non_project_count, type_count,
     project_file, nonproject_file, types_file) = write_filtered_results(results, parent_dir, logger)
    
    logger.debug("\nProcessing complete:")
    logger.debug("- %d samples written to '%s'", project_count, project_file)
    logger.debug("- %d samples written to '%s'", non_project_count, nonproject_file)
    logger.debug("- %d type specimens written to '%s'", type_count, types_file)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        parent_dir = sys.argv[1]
        sample_metadata_file = sys.argv[2]
        main(parent_dir, sample_metadata_file)
    else:
        print(
        """
        Usage: python 4_sample_sheet.py [parent_directory] [sample_metadata_file]
        parent_directory: parent directory containing either:
          - subdirectories with raw read pairs, or
          - raw read pairs directly in the parent directory
        sample_metadata_file: Path to the sample_metadata.csv file.
        Example: python script.py /path/to/raw_reads /path/to/sample_metadata.csv
        """
        )
