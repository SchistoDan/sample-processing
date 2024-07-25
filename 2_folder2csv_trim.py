import os
import csv
import sys

def find_files(folder_path):
#initialize dir to store results
    results = {}

#walk through subdirs looking for _R1_ and _R2_ files
    for root, _, files in os.walk(folder_path):
        r1_path, r2_path = None, None
        for filename in files:
            if "_R1_" in filename:
                r1_path = os.path.abspath(os.path.join(root, filename))
            elif "_R2_" in filename:
                r2_path = os.path.abspath(os.path.join(root, filename))

        if r1_path and r2_path:
            subfolder_name = os.path.basename(root)

#clean subdir name by removing first 15 characters
            subfolder_name = subfolder_name[15:]
            results[subfolder_name] = (r1_path, r2_path)

    return results



def write_to_csv(results, output_filename):
#write R1 nd R2 file paths to csv
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'forward', 'reverse'])
        for subfolder, (r1_path, r2_path) in results.items():
            writer.writerow([subfolder, r1_path, r2_path])



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 1_folder2csv_trim.py /path/to/dir/raw_seqs.fq.gz")
        sys.exit(1)

    main_folder_path = sys.argv[1]

#get base dir name of main_folder_path
    base_dir = os.path.basename(os.path.normpath(main_folder_path))
    output_filename = f"{base_dir}_read_paths.csv"

    files_info = find_files(main_folder_path)
    write_to_csv(files_info, output_filename)
    print(f"CSV file '{output_filename}' created successfully.")
