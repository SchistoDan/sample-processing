Collection of scripts for processing BOLD-downloaded sample metadata, collects other useful sample-related information, and generates input files for downstream processes ([uploading trimmed reads to ENA](https://github.com/bge-barcoding/ena-read-upload), input into [skim2mito](https://github.com/o-william-white/skim2mito) and [MGE](https://github.com/bge-barcoding/MitoGeneExtractor-BGE), and [requesting taxid creation from ENA](https://github.com/bge-barcoding/ena-taxid-creation)). Written by Ben Price and Dan Parsons @ NHMUK.
- For help, see usage information and docstrings within each script.

## 1_combine_tsv.py ##
- Merges TSV files from multiple zipped folders downloaded from BOLD.
  - **`usage: python combine_tsv.py <output_dir> <zip_files_dir>`**
  - `output_dir: Directory where combined .tsv files will be stored.`
  - `zip_files_dir: Directory containing the zipped files to be merged.
    - Recommended directory structure:parent_directory/
```├── parent_file1.txt
├── collection_data.tsv
├── lab.tsv
├── merged_custom_fields.tsv
├── specimen_details.tsv
├── tags.tsv
├── taxonomy.tsv
├── voucher.tsv
│
├── subdirectory.zip
│   ├── collection_data.tsv
│   ├── lab.tsv
│   ├── merged_custom_fields.tsv
│   ├── specimen_details.tsv
│   ├── tags.tsv
│   ├── taxonomy.tsv
│   └── voucher.tsv
│
└── subdirectory2.zip
    ├──  collection_data.tsv
    ├── lab.tsv
    ├── merged_custom_fields.tsv
    ├── specimen_details.tsv
    ├── tags.tsv
    ├── taxonomy.tsv
    └── voucher.tsv````
 
## 2_download_taxonomy.sh ##
- Downloads the [newst NCBI taxonomy (i.e. new_taxdump)](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/).
- If the taxdump was previously downloaded to the same directory, the script will replace that taxdump with the newst version. The tar.gz download will also be removed post-extraction.
  - **`usage: bash/sbatch/srun download_taxonomy.sh <output_file> <extract_dir>`**
  - `output_file = Name of NCBI taxdump (e.g. new_taxdump.tar.gz)`
  - `extract_dir = Directory to output taxdump files to (e.g. new_taxdump). If given a relative or absolute path, new directories will be created.`
 


## 3_local_bold_processing.py ##
- Merges relevant sample metadata from BOLD .tsv files, and resolves the taxonomic ID using the hierarchical structure and checks for mismatches between BOLD and NCBI taxonomy. Outputs sample_metadata.csv containing fields below:
  - Sample ID
  - Process ID
  - BOLD taxonomic ranks (phylum->species)
  - taxid
  - matched_rank (taxonimic rank the taxid corresponds to)
  - lineage (full NCBI lineage for taxid)
  - lineage_mismatch (did the major BOLD taxonomic ranks match the fetched NCBI lineage)
  - BOLD sample metadata (Identifier (identified_by), Collection Date (collection_date), Geographic Location (Country/Ocean) (geographic_location), Exact Site (geographic_location_locality), Latitude (latitude), Longitude (longitude), Collected By (collected_by), Habitat (habitat), Tissue Descriptor (organism_part), Sex (sex), Life stage (lifestage), Museum ID (specimen_voucher), Institution Storing (collecting_institution).
- **`usage: python local_bold_processing.py <input_dir> <rankedlineage_path> <output_dir>`**
  - `input_dir = Directory containing BOLD-downloaded sample metadata (.tsv files).`
  - `rankedlineage_path = Path to NCBI taxonomic hierarchy/lineage (new_taxdump file).`
  - `output_dir = Directory to output sample_metadata.tsv to.`
 
If metadata was not collected for a particular sample, 'not collected' output to field as required by [ToL ENA sample registration checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000053).

**Example samples_metadata.csv**
| Sample ID | Process ID  | Phylum | Class | Order | Family | Subfamily | Genus | Species  | taxid | matched_rank | specimen_voucher | lifestage | collection_date | geographic_location | geographic_location_locality | latitude | longitude | collected_by | habitat | identified_by | collecting_institution | organism_part | sex |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |--- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |
| BGE_0001_A01  | BSNHM001-24 | Arthropoda | Insecta | Trichoptera | Apataniidae | Apataniinae | Apatania | Apatania stylata | 177658 | genus | 'Museum ID' | adult | YYYY-MM-DD | France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Whole | M |
| BGE_0001_A02 | BSNHM002-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | Agapetus | Agapetus iridipennis | 177627 | genus | 'Museum ID' | adult | YYYY-MM-DD | Switzerland | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | Naturalis | not collected | F |
| BGE_0001_A03 | BSNHM003-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Hydropsychidae | Diplectrona | Diplectrona meridionalis | 177860 | genus | 'Museum ID' | adult | YYYY-MM-DD |  France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Leg | M |
 
## 4_samples_spreadsheet.py ##
Script to generate samples.csv requried to run MGE and skim2mito piplines. Script also screen BGE project samples using the 27 BGE project codes.
- **`usage: python 2_samples_spreadsheet.py [path/to/raw/read/dir] [path/to/samples_metadata.csv]`**
**- `path/to/raw/read/dir: Path to parent directory of subdirectories containing raw PE read files.`
- `path/to/output/dir/samples_metadata.csv`: Directory containing sample_metadta.csv file.`
- outputs two CSV files containing ID (Process ID), forward (absolute path to R1 read (fastq.gz), reverse (absolute path to R2 read (fastq.gz) and taxid to current directory.
  - samples_[parent_dir_name].csv = contains all samples containing BGE project codes in their Process ID's.
  - samples_[parent_dir_name]_nonproject.csv = contains all samples not containing BGE project codes.

## Miro workflow ##
![image](![image](https://github.com/user-attachments/assets/dd960f04-2957-4ef1-bd58-d5457ecba777)
