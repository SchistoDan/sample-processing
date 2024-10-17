import os
import sys
import zipfile
import tempfile
import shutil
import glob
import csv

def process_zip(zip_file, output_dir):
    # Create temp dir for file extraction
    with tempfile.TemporaryDirectory() as temp_dir:
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
        
        # Process all extracted TSV files
        for tsv_file in glob.glob(os.path.join(temp_dir, '**', '*.tsv'), recursive=True):
            filename = os.path.basename(tsv_file)
            output_file = os.path.join(output_dir, filename)
            
            with open(tsv_file, 'r', newline='') as infile:
                reader = csv.reader(infile, delimiter='\t')
                headers = next(reader)  # Read the header row
                
                if not os.path.exists(output_file):
                    # If output file doesn't exist, write headers and all data
                    with open(output_file, 'w', newline='') as outfile:
                        writer = csv.writer(outfile, delimiter='\t')
                        writer.writerow(headers)
                        writer.writerows(reader)
                else:
                    # If output file exists, append data without headers
                    with open(output_file, 'a', newline='') as outfile:
                        writer = csv.writer(outfile, delimiter='\t')
                        writer.writerows(reader)



def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python combine_tsv.py <output_dir> <zip_files_dir>")
        sys.exit(1)

    output_dir = sys.argv[1]
    zip_dir = sys.argv[2]



    # Validate output directory
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error: Unable to create output directory {output_dir}. {e}")
            sys.exit(1)
    
    # Validate zip files directory
    if not os.path.isdir(zip_dir):
        print(f"Error: {zip_dir} is not a valid directory.")
        sys.exit(1)



    # Find all .zip files in the provided zip files directory
    zip_files = glob.glob(os.path.join(zip_dir, '*.zip'))
    if not zip_files:
        print(f"No zip files found in {zip_dir}.")
        sys.exit(1)

    # Process each zip file
    for zip_file in zip_files:
        print(f"Processing {zip_file}...")
        process_zip(zip_file, output_dir)

    print(f"TSV files have been combined in {output_dir}")



if __name__ == "__main__":
    main()
