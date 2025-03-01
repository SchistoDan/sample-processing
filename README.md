Collection of scripts for processing and parsing necessary BOLD-downloaded sample metadata, and generating input files for downstream processes (e.g. ([uploading trimmed reads to ENA](https://github.com/bge-barcoding/ena-read-upload), input into [skim2mito](https://github.com/o-william-white/skim2mito) and [MGE](https://github.com/bge-barcoding/MitoGeneExtractor-BGE), and [requesting taxid creation from ENA](https://github.com/bge-barcoding/ena-taxid-creation))). 
- Written by Ben Price and Dan Parsons @ NHMUK.
- For help, see usage information and docstrings within each script.
- Requires sample_processing conda environment to be activated. Environment with all necessary dependencies can be created from sample_processing.yaml file in this repo.

## 1_combine_tsv.py ##
- Merges TSV files from multiple zipped folders downloaded from BOLD.
  - **`usage: python combine_tsv.py <output_dir> <zip_files_dir>`**
  - `output_dir: Directory where combined .tsv files will be stored.`
  - `zip_files_dir: Directory containing the zipped files to be merged.
    - Recommended directory structure:
```
parent_directory/
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
    └── voucher.tsv
```
 
## 2_download_taxonomy.sh ##
- Downloads the [newst NCBI taxonomy (i.e. new_taxdump)](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/).
- If the taxdump was previously downloaded to the same directory, the script will replace that taxdump with the newst version. The tar.gz download will also be removed post-extraction.
  - **`usage: srun download_taxonomy.sh <output_file> <extract_dir>`**
  - `output_file = Name of NCBI taxdump (e.g. <date>_taxdump.tar.gz)`
  - `extract_dir = Directory to output taxdump files to (e.g. <date>_taxdump). If given a relative or absolute path, new directories will be created.`
 


## 3_local_bold_processing.py ##
- Merges relevant sample metadata from BOLD .tsv files, and resolves the taxonomic ID using the hierarchical structure and checks for mismatches between BOLD and NCBI taxonomy. Outputs sample_metadata.csv containing fields below:
  - Sample ID
  - Process ID
  - BOLD taxonomic ranks (phylum->species)
  - taxid
  - matched_rank (taxonimic rank the taxid corresponds to)
  - lineage (full NCBI lineage for taxid)
  - lineage_mismatch (did the major BOLD taxonomic ranks match the fetched NCBI lineage)
  - BOLD sample metadata (Identifier (identified_by), Collection Date (collection_date), Geographic Location (Country/Ocean) (geographic_location), Exact Site (geographic_location_locality), Latitude (latitude), Longitude (longitude), Collected By (collected_by), Habitat (habitat), Tissue Descriptor (organism_part), Sex (sex), Life stage (lifestage), Museum ID (specimen_voucher), Institution Storing (collecting_institution), Type Status (type_status).
- **`usage: python local_bold_processing.py <input_dir> <rankedlineage_path> <output_dir>`**
  - `input_dir = Directory containing BOLD-downloaded sample metadata (.tsv files).`
  - `rankedlineage_path = Path to NCBI taxonomic hierarchy/lineage (<date>_taxdump/rankedlineage.dmp file).`
  - `output = Directory to output sample_metadata.csv to. Provide name of output .csv file`
 
If metadata was not collected for a particular sample, 'not collected' output to field as required by [ToL ENA sample registration checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000053).

**Example samples_metadata.csv**
| Sample ID | Process ID  | Phylum | Class | Order | Family | Subfamily | Genus | Species  | taxid | matched_rank | specimen_voucher | lifestage | collection_date | geographic_location | geographic_location_locality | latitude | longitude | collected_by | habitat | identified_by | collecting_institution | organism_part | sex | type_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |--- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- | --- |
| BGE_0001_A01  | BSNHM001-24 | Arthropoda | Insecta | Trichoptera | Apataniidae | Apataniinae | Apatania | Apatania stylata | 177658 | genus | 'Museum ID' | adult | YYYY-MM-DD | France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Whole | M | type |
| BGE_0001_A02 | BSNHM002-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | Agapetus | Agapetus iridipennis | 177627 | genus | 'Museum ID' | adult | YYYY-MM-DD | Switzerland | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | Naturalis | not collected | F | Paratype |
| BGE_0001_A03 | BSNHM003-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Hydropsychidae | Diplectrona | Diplectrona meridionalis | 177860 | genus | 'Museum ID' | adult | YYYY-MM-DD |  France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Leg | M | no |

- This can also be run using the 3_local_bold_proccessing.sh script in this repo to run it on a slurm cluster. This will greatly speed up creation of the sample_metadata.csv output file.

## 4_samples_spreadsheet.py ##
Script to generate samples.csv requried to run MGE and skim2mito piplines. Script also screen BGE project samples using the 27 BGE project codes.
- **`usage: python 2_samples_spreadsheet.py [path/to/raw/read/dir] [path/to/samples_metadata.csv]`**
- `path/to/raw/read/dir: Path to parent directory with 'flat' or 'nested' structure (see below/docstring for more example) containing raw PE read files.`
- `path/to/output/dir/samples_metadata.csv`: Directory containing sample_metadata.csv file.`
- outputs three CSV files containing ID (Process ID), forward (absolute path to R1 read (fastq.gz), reverse (absolute path to R2 read (fastq.gz) and taxid to current directory. samples_[parent_dir_name]_types.csv also contains type_status field, for reference.
  - samples_[parent_dir_name].csv = contains all samples containing BGE project codes in their Process ID's.
  - samples_[parent_dir_name]_nonproject.csv = contains all samples not containing BGE project codes.
  - samples_[parent_dir_name]_types.csv = contains all samples (project or non-project) that have 'type' in Type Status field.
```
Example Directory Structures Supported (where BSNHM593-24 is the process ID):
1. Flat structure:
   /parent_dir/
   ├── BSNHM593-24_R1.fq
   ├── BSNHM593-24_R2.fq
   └── ...

2. Nested structure:
   /parent_dir/XE-4013/
   └── 20240906_LH00179_0123_A22CKGHLT4/
       ├── Sample_XE-4013-BSNHM593-24/
       │   ├── BSNHM593-24_R1.fastq.gz
       │   └── BSNHM593-24_R2.fastq.gz
       └── ...
```

## Miro workflow ##
![image](![image](https://github.com/user-attachments/assets/dd960f04-2957-4ef1-bd58-d5457ecba777)
