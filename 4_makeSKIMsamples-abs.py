import sys
import os
import pandas as pd




def merge_for_skim2mt_in(input_file1, input_file2):
#Read data from csv files
    try:
        df1 = pd.read_csv(input_file1)  #File 1 = [run]_read_paths.csv
        df2 = pd.read_csv(input_file2)  #File 2 = sample2taxid_out.csv
    except Exception as e:
        print(f"Error reading input files: {e}")
        sys.exit(1)



#Debug - print first few rows of input dfs
    print("DataFrame 1 (df1) head:\n", df1.head())
    print("DataFrame 2 (df2) head:\n", df2.head())

#Check columns available in each df
    print("Columns in df1:", df1.columns)
    print("Columns in df2:", df2.columns)



#Make sure expected columns are present
    if 'ID' not in df1.columns:
        print("Error: 'ID' column not found in DataFrame 1")
        sys.exit(1)
    if 'Process ID' not in df2.columns:
        print("Error: 'Process ID' column not found in DataFrame 2")
        sys.exit(1)



#Merge dfs on matching values in "ID" column of df1 and the "Process ID" column of df2 using inner join (meaning only matching values will be included in the merged_df)
    merged_df = pd.merge(df1, df2, left_on="ID", right_on="Sample ID", how="inner")



#Debug - print merged df
    print("Merged DataFrame head:\n", merged_df.head())



#Select only desired columns from 2 csv files
    selected_columns = ["Process ID", "forward", "reverse", "taxid"]
    merged_df = merged_df[selected_columns]



#Debug - print selected columns of merged df
    print("Selected Columns DataFrame head:\n", merged_df.head())



#Rename 'Process ID' column to 'ID' (nneded for skim2mito input)
    merged_df = merged_df.rename(columns={"Process ID": "ID"})



#Debug - print df after converting paths
    print("DataFrame after path conversion head:\n", merged_df.head())



#get dir and filename from input_file1, create output filename based on the first part of input_file1 before '_' and add _samples.csv
    output_dir = os.path.dirname(input_file1)
    base_filename = os.path.basename(input_file1).split('_')[0]
    output_filename = f"{base_filename}_samples.csv"
    output_path = os.path.join(output_dir, output_filename)



#write merged df to output csv
    merged_df.to_csv(output_path, index=False)
    print(f"Output written to {output_path}")



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 3_makeSKIMsamples.py input_file1.csv input_file2.csv")
        print("input_file1 = [run_dir]_read_paths.csv input_file2 = sample2taxid_out.csv")
    else:
        input_file1, input_file2 = sys.argv[1], sys.argv[2]
        merge_for_skim2mt_in(input_file1, input_file2)

    print(f"skim2mito input csv created successfully.")

