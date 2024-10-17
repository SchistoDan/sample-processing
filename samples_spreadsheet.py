import os
import csv
import sys
import pandas as pd





def extract_id_reads_taxid(parent_dir, process_ids, taxid_dict):
    """
    Finds paired raw read files (R1 and R2) and corresponding Process IDs in the given directory and its subdirectories.
    Skips subdirs containing "Undetermined" or "NC".
    
    Args:
    parent_dir (str): Path to the parent dir containing subdirs with raw read files.
    process_ids (set): Set of valid Process IDs to match against.
    taxid_dict (dict): Dictionary of Process IDs to taxid.

    Returns:
    dict: A dictionary where keys are Process IDs and values are tuples of paths to R1, R2 files, and taxid.
    """
    results = {}
    for root, dirs, files in os.walk(parent_dir):
        for subdir in dirs:
            if "Undetermined" in subdir or "NC" in subdir:
                continue  # Skip subdirs containing "Undetermined" or "NC"
                
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

            taxid = taxid_dict.get(process_id, 'NA')  # Get taxid from dictionary or 'NA' if not found

            if r1_path and r2_path:
                results[process_id] = (r1_path, r2_path, taxid)

    return results





def write_to_csv(results, output_filename):
    """
    Writes raw read file paths and taxid to .csv file.

    Args:
    results (dict): A dictionary where keys are Process IDs and values are tuples of paths to R1, R2 files, and taxid.
    output_filename (str): Path to the output .csv file.
    """
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse', 'taxid'])  
        for process_id, (r1_path, r2_path, taxid) in results.items():
            writer.writerow([process_id, r1_path, r2_path, taxid]) 







def main(parent_dir, sample_metadata_file):
    """
    Main function processing raw read files and generates samples.csv with taxid.

    Args:
    parent_dir (str): Directory containing raw read files.
    sample_metadata_file (str): Path to the sample metadata CSV file.
    """

    #Read Process IDs and taxid from sample_metadata.csv
    metadata_df = pd.read_csv(sample_metadata_file)
    process_ids = set(metadata_df['Process ID'].astype(str))
    taxid_dict = dict(zip(metadata_df['Process ID'].astype(str), metadata_df['taxid']))


    #Find files (forward/reverse reads)
    results = extract_id_reads_taxid(parent_dir, process_ids, taxid_dict)
    

    #Write to samples.csv
    samples_file = "samples.csv"
    write_to_csv(results, samples_file)
    

    print(f"Samples CSV file '{samples_file}' created successfully.")





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

