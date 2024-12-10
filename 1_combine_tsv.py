"""Mixed Directory TSV Processor and Combiner

This script processes and combines TSV files from both:
- Unzipped TSV files in the input directory and its subdirectories
- TSV files contained within ZIP archives anywhere in the directory tree

The script maintains header consistency while combining files with matching names,
handling both zipped and unzipped sources uniformly.

Usage:
    python combine_tsv.py <output_dir> <input_dir>

Arguments:
    output_dir (str): Path where combined TSV files will be saved
    input_dir (str): Root directory containing both TSV and ZIP files to process
"""

import os
import sys
import zipfile
import tempfile
import shutil
import glob
import csv
import logging
import time
from typing import Set, Dict, List
from datetime import datetime

class ProgressTracker:
    def __init__(self):
        self.total_files = 0
        self.processed_files = 0
        self.processed_zips = 0
        self.error_count = 0
        self.file_counts: Dict[str, int] = {}  # Tracks rows per output file
        self.start_time = time.time()
        
    def increment_processed(self):
        self.processed_files += 1
        
    def increment_zip(self):
        self.processed_zips += 1
        
    def increment_error(self):
        self.error_count += 1
        
    def update_file_count(self, filename: str, rows: int):
        self.file_counts[filename] = self.file_counts.get(filename, 0) + rows
        
    def get_elapsed_time(self) -> str:
        elapsed = time.time() - self.start_time
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        return f"{minutes}m {seconds}s"
        
    def print_summary(self):
        logging.info("=" * 50)
        logging.info("Processing Summary:")
        logging.info("-" * 50)
        logging.info(f"Total files processed: {self.processed_files}")
        logging.info(f"ZIP archives processed: {self.processed_zips}")
        logging.info(f"Errors encountered: {self.error_count}")
        logging.info(f"Total processing time: {self.get_elapsed_time()}")
        logging.info("\nOutput files summary:")
        for filename, count in self.file_counts.items():
            logging.info(f"  {filename}: {count:,} rows")
        logging.info("=" * 50)

def setup_logging(output_dir: str):
    """Set up logging configuration."""
    log_file = os.path.join(output_dir, f"tsv_processor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def count_rows(reader) -> tuple[List, int]:
    """Count rows while preserving the data."""
    data = list(reader)
    return data, len(data)

def process_tsv_file(tsv_path: str, output_dir: str, processed_files: Set[str], progress: ProgressTracker) -> None:
    """
    Process a single TSV file, either creating a new output file or appending to existing one.
    
    Args:
        tsv_path: Path to the TSV file to process
        output_dir: Directory where processed files should be saved
        processed_files: Set of already processed file paths to avoid duplicates
        progress: ProgressTracker instance for monitoring progress
    """
    if tsv_path in processed_files:
        return
    
    filename = os.path.basename(tsv_path)
    output_file = os.path.join(output_dir, filename)
    
    try:
        with open(tsv_path, 'r', newline='') as infile:
            reader = csv.reader(infile, delimiter='\t')
            headers = next(reader)  # Read the header row
            
            # Count rows while preserving the data
            data, row_count = count_rows(reader)
            
            if not os.path.exists(output_file):
                logging.info(f"Creating new file: {filename}")
                with open(output_file, 'w', newline='') as outfile:
                    writer = csv.writer(outfile, delimiter='\t')
                    writer.writerow(headers)
                    writer.writerows(data)
            else:
                logging.info(f"Appending to existing file: {filename}")
                with open(output_file, 'a', newline='') as outfile:
                    writer = csv.writer(outfile, delimiter='\t')
                    writer.writerows(data)
            
            progress.update_file_count(filename, row_count)
            processed_files.add(tsv_path)
            progress.increment_processed()
            logging.info(f"Processed {tsv_path} ({row_count:,} rows)")
            
    except Exception as e:
        progress.increment_error()
        logging.error(f"Error processing {tsv_path}: {e}")

def process_zip(zip_path: str, output_dir: str, processed_files: Set[str], progress: ProgressTracker) -> None:
    """
    Process TSV files within a ZIP archive.
    
    Args:
        zip_path: Path to the ZIP file
        output_dir: Directory where processed files should be saved
        processed_files: Set of already processed file paths to avoid duplicates
        progress: ProgressTracker instance for monitoring progress
    """
    try:
        logging.info(f"Processing ZIP file: {zip_path}")
        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Process all extracted TSV files
            tsv_files = glob.glob(os.path.join(temp_dir, '**', '*.tsv'), recursive=True)
            for tsv_file in tsv_files:
                process_tsv_file(tsv_file, output_dir, processed_files, progress)
            
        progress.increment_zip()
        logging.info(f"Completed processing ZIP file: {zip_path}")
        
    except Exception as e:
        progress.increment_error()
        logging.error(f"Error processing ZIP file {zip_path}: {e}")

def process_directory(input_dir: str, output_dir: str) -> None:
    """
    Process all TSV and ZIP files in the input directory and its subdirectories.
    
    Args:
        input_dir: Root directory to process
        output_dir: Directory where processed files should be saved
    """
    progress = ProgressTracker()
    processed_files = set()
    
    # Count total files for progress tracking
    tsv_files = glob.glob(os.path.join(input_dir, '**', '*.tsv'), recursive=True)
    zip_files = glob.glob(os.path.join(input_dir, '**', '*.zip'), recursive=True)
    
    logging.info(f"Found {len(tsv_files)} TSV files and {len(zip_files)} ZIP files to process")
    
    # First, process all unzipped TSV files
    for tsv_file in tsv_files:
        process_tsv_file(tsv_file, output_dir, processed_files, progress)
    
    # Then process all ZIP files
    for zip_file in zip_files:
        process_zip(zip_file, output_dir, processed_files, progress)
    
    # Print summary
    progress.print_summary()

def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python combine_tsv.py <output_dir> <input_dir>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    input_dir = sys.argv[2]
    
    # Validate output directory
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error: Unable to create output directory {output_dir}. {e}")
            sys.exit(1)
    
    # Validate input directory
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)
    
    # Set up logging
    setup_logging(output_dir)
    
    # Log start of processing
    logging.info(f"Starting TSV processing")
    logging.info(f"Input directory: {input_dir}")
    logging.info(f"Output directory: {output_dir}")
    
    # Process the directory
    process_directory(input_dir, output_dir)

if __name__ == "__main__":
    main()
