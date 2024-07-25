import sys
import pandas as pd

def combine_tsv_files(lab_file, taxonomy_file, output_file, unmatched_file):
    #Read the data from the two TSV files
    lab_df = pd.read_csv(lab_file, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')

    #Merge the dataframes on the 'Sample ID' column using an inner join & only take Sample ID and Process ID from lab.tsv and taxonomy info from taxonomy.tsv
    merged_df = pd.merge(lab_df[['Sample ID', 'Process ID']], taxonomy_df[['Sample ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']], on='Sample ID', how='inner')

    #Find unmatched Sample IDs
    unmatched_ids_lab = lab_df[~lab_df['Sample ID'].isin(taxonomy_df['Sample ID'])]
    unmatched_ids_taxonomy = taxonomy_df[~taxonomy_df['Sample ID'].isin(lab_df['Sample ID'])]

    #Combine unmatched Sample IDs into a single dataframe
    unmatched_ids = pd.concat([unmatched_ids_lab['Sample ID'], unmatched_ids_taxonomy['Sample ID']]).unique()

    #Write the merged dataframe to a new CSV file
    merged_df.to_csv(output_file, index=False)
    print(f"Combined file written to {output_file}")

    #Write unmatched Sample IDs to a text file
    with open(unmatched_file, 'w') as txt_file:
        for sample_id in unmatched_ids:
            txt_file.write(f"{sample_id}\n")
    print(f"Unmatched Sample IDs written to {unmatched_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python combine_tsv.py /path/to/lab.tsv /path/to/taxonomy.tsv BOLD_output.csv unmatched_sampleids.txt")
    else:
        lab_file, taxonomy_file, output_file, unmatched_file = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
        combine_tsv_files(lab_file, taxonomy_file, output_file, unmatched_file)
