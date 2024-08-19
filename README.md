Series of python scripts to take BOLD download and paths to raw reads and create sample submission CSV for input into [skim2mito](https://github.com/o-william-white/skim2mito)



## 1_combine_tsv.py
Takes BOLD download (lab.tsv and taxonomy.tsv), merges two .tsv files based on Sample ID, and parses Sample ID, Process ID, and taxonomic rankings (Phylum, Class, order, family, Genus, and Species) to .csv file.

**usage: python 1_combine_tsv.py [/path/to/lab.tsv] [/path/to/taxonomy.tsv] [/path/to/BOLD_output.csv] [/path/to/unmatched_sampleids.txt]**
- lab.tsv = contains Process ID and Sample ID.
- taxonomy.tsv = contains Sample ID and taxonomic rankings.
- [BOLD_output].csv = can be named anything.
- [unmatched_sampleids].txt = Can be named anything. Should always be empty.




## 2_folder2csv_trim.py
Takes path to directory containing raw reads (.fq.gz) and populates .csv with sample names (ID), R1 read path (forward) and R2 read path (reverse).

**usage: python 2_folder2csv_trim.py [/path/to/dir/raw_seqs.fq.gz]**
- /path/to/dir/raw_seqs.fq.gz = path to parent directory containing raw reads in .fq.gz format.
  
**output**
- [trimmed_parent_dir_name]_read_paths.csv

| ID  | forward | reverse |
| --------- | --------- | --------- |
| BGE_0001_A01  | path/to/R1.fq.gz  | path/to/R2.fq.gz |
| BGE_0001_A02 | path/to/R1.fq.gz  | path/to/R2.fq.gz |
| BGE_0001_A03 | path/to/R1.fq.gz | path/to/R2.fq.gz |





## 3_sample2taxid.py
Takes BOLD_output.csv containing sample ID, Process ID, and taxonomic ranks based on BOLD morphological identification (Phylum, class, Order, Family, Subfamily, Tribe, Genus, Species, Subspecies), grabs taxonomic ID (taxid) using NCBI Entrez API, determines taxonomic rank of taxid, parses taxid and matched_rank to .csv, and grabs other fields from BOLD downlaod that are required for downstream sample registration (generation of sample accession numbers) using ENA Webin.

**usage: python 3_sample2taxid.py [path/to/BOLD_output.csv] [path/to/directory_containing_tsv_files] [sample2taxid_out.csv]**
- path/to/[BOLD_output].csv = path to directory containing [BOLD_output].csv output by '1_combine_tsv.py'.
- path/to/directory_containing_tsv_files = path to directory containing BOLD container dowlonad (.tsv files - taxonomy, lab, collection_data, voucher, specimen_details).
- sample2taxid_out.csv = can be named anuthing. Contains parsed fields below.
- [BOLD_output]_unique_taxids.txt = list of deuplicated taxids.

| Sample ID | Process ID  | Phylum | Class | Order | Family | Subfamily | Tribe | Genus | Species | Subspecies | taxid | matched_rank | specimen_voucher | lifestage | collection date | geographic_location (country and/or sea) | geographic_location (region and locality) | latitude | longitude | collected_by | habitat | identified_by | collection institution | organism part | sex |
| --- | --- |--- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |--- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |
| BGE_0001_A01  | BSNHM001-24 | Arthropoda | Insecta | Trichoptera | Apataniidae | Apataniinae | | Apatania | Apatania stylata | | 177658 | genus | 'Museum ID' | adult | YYYY-MM-DD | France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Whole | M |
| BGE_0001_A02 | BSNHM002-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | | Agapetus | Agapetus iridipennis | | 177627 | genus | 'Museum ID' | adult | YYYY-MM-DD | Switzerland | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | Naturalis | not collected | F |
| BGE_0001_A03 | BSNHM003-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Hydropsychidae | | Diplectrona | Diplectrona meridionalis | | 177860 | genus | 'Museum ID' | adult | YYYY-MM-DD |  France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Leg | M |

- If data not collected for sample, 'not collected' output to field.




## 4_makeSKIMsample.py
Takes Process ID, forward and reverse read paths, and taxid, and parses them to .csv ready for input into skim2mito.
**4_makeSKIMsample-rel.py** = outputs relative paths to reads (currently needed by skim2mito - needs to be changed to abs).
**4_makeSKIMsample-abs.py** = outputs absolute paths to reads (needed by MitoGeneExtractor).


**usage: python 4_makeSKIMsamples.py [input_file1.csv] [input_file2.csv]**
- input_file1.csv = [trimmed_parent_dir_name]_read_paths.csv
- input_file1.csv = [sample2taxid].csv from '2_sample2taxid.py'.
- [trimmed_parent_dir_name]_samples.csv = 'samples.csv' for input into skim2mito and MGE (via snakemake).

| ID | forward | reverse | taxid |
| --------- | --------- |--------- | --------- |
| BSNHM001-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |


## TO DO ##
- Stitch these scripts together so that three input files can be given:
  - lab.tsv
  - taxonomy.tsv
  - /path/to/dir/raw_seqs.fq.gz
- It would still be necessary for all 'intermediate' files to be generated, however, as some are used for downstream processes (e.g. sample2taxid output.csv).
