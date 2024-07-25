Series of python scripts to take BOLD-downloaded taxonomy and paths to raw reads and create sample submission CSV for input into [skim2mito](https://github.com/o-william-white/skim2mito)



## 1_combine_tsv.py
Takes BOLD DB download, merges two .tsv files based on Sample ID, and parses Sample ID, Process ID, and taxonomic rankings (Phylum, Class, order, family, Genus, and Species) to .csv file

**usage: python 1_combine_tsv.py [/path/to/lab.tsv] [/path/to/taxonomy.tsv] [/path/to/BOLD_output.csv] [/path/to/unmatched_sampleids.txt]**
- lab.tsv = contains Process ID and Sample ID.
- taxonomy.tsv = contains Sample ID and taxonomic rankings.
- BOLD_output.csv = can be named anything.
- unmatched_sampleids.txt = Can be named anything. Should always be empty.



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
Takes BOLD taxonomy.tsv output containing sample ID and taxonomic ranks based on morphological identification (Phylum, class, Order, Family, Subfamily, Tribe, Genus, Species, Subspecies), grabs taxonomic ID (taxid) using NCBI Entrez API, determines taxonomic rank of taxid (matched_rank), and parses taxid and matched_rank to .csv.

**usage: python 3_sample2taxid.py [path/to/BOLD_output.csv] [output.csv]**
- path/to/BOLD_output.csv = path to directory containing sample taxonomy information downloaded from BOLD (must specify filename).
- output.csv = name of .csv file containing parsed fields (see below).

| Process ID  | Phylum | Class | Order | Family | Subfamily | Tribe | Genus | Species | Subspecies | taxid | matched_rank |
| --------- | --------- |--------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- |
| BGE_0001_A01  | Arthropoda | Insecta | Trichoptera | Apataniidae | Apataniinae | | Apatania | Apatania stylata | | 177658 | genus |
| BGE_0001_A02 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | | Agapetus | Agapetus iridipennis | | 177627 | genus |
| BGE_0001_A03 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | | Agapetus | Agapetus incertulus | | 3084599 | species |

Also outputs a .txt file containing a list of deuplicated taxids.

## 4_makeSKImsample.py
Takes Process ID, forward and reverse read paths, and taxid and parses them to .csv ready for input into skim2mito

**usage: python 4_makeSKIMsamples.py [input_file1.csv] [input_file2.csv]**
- input_file1.csv = [trimmed_parent_dir_name]_read_paths.csv
- input_file1.csv = output.csv from 2_sample2taxid.py

**output**
- [trimmed_parent_dir_name]_samples.csv

| ID | forward | reverse | taxid |
| --------- | --------- |--------- | --------- |
| BGE_0001_A01  | relative/path/to/R1.fq.gz | relative/path/to/R2.fq.gz | 177658 |
| BGE_0001_A02 | relative/path/to/R1.fq.gz | relative/path/to/R2.fq.gz | 177627 |
| BGE_0001_A03 | relative/path/to/R1.fq.gz | relative/path/to/R2.fq.gz | 3084599 |
