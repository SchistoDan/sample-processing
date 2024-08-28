Python script that takes BOLD container download (.tsv files) and paths to raw reads, fetches taxonomic ID (taxid) using taxonomic ranks for each sample, merges relevant fields into output ([output].csv) from BOLD-downlaoded files needed for downstream [upload of reads to ENA](https://github.com/bge-barcoding/ena-read-upload), creates analysis pipeline sample submission CSV (samples.csv) for input into [skim2mito](https://github.com/o-william-white/skim2mito) and [MGE](https://github.com/bge-barcoding/MitoGeneExtractor-BGE), and compiles list of deduplicated taxid's (unique_taxids.txt)



## sample_processing.py ##
**usage: python sample_processing.py [path/to/BOLD/download/dir] [path/to/raw/read/dir] [path/to/output/dir/[output].csv]**
- path/to/BOLD/download/dir: path to directory containing .tsv files (voucher.tsv, collection_data.tsv, specimen_details.tsv, taxonomy.tsv, and lab.tsv) downloaded from BOLD DB.
- path/to/raw/read/dir: Path to directory containing subdirectories with raw PE reads files. Give parent directory of subdirectories as input here.
- path/to/output/dir/[output].csv: Directory to create [output].csv, samples.csv and unique_taxids.txt within. Need to name [output].csv, whereas samples.csv and unique_taxids.txt are named as is.



**Example samples.csv**
| ID  | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM001-24  | abs/path/to/R1.fq.gz  | abs/path/to/R2.fq.gz | 1234 |
| BSNHM001-24 | abs/path/to/R1.fq.gz  | abs/path/to/R2.fq.gz | 5678 |
| BSNHM001-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz |  91011 |



**Example [output].csv**
| Sample ID | Process ID  | Phylum | Class | Order | Family | Subfamily | Tribe | Genus | Species | Subspecies | taxid | matched_rank | specimen_voucher | lifestage | collection date | geographic_location (country and/or sea) | geographic_location (region and locality) | latitude | longitude | collected_by | habitat | identified_by | collection institution | organism part | sex |
| --- | --- |--- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |--- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |
| BGE_0001_A01  | BSNHM001-24 | Arthropoda | Insecta | Trichoptera | Apataniidae | Apataniinae | | Apatania | Apatania stylata | | 177658 | genus | 'Museum ID' | adult | YYYY-MM-DD | France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Whole | M |
| BGE_0001_A02 | BSNHM002-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | | Agapetus | Agapetus iridipennis | | 177627 | genus | 'Museum ID' | adult | YYYY-MM-DD | Switzerland | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | Naturalis | not collected | F |
| BGE_0001_A03 | BSNHM003-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Hydropsychidae | | Diplectrona | Diplectrona meridionalis | | 177860 | genus | 'Museum ID' | adult | YYYY-MM-DD |  France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Leg | M |

- If data not collected for sample, 'not collected' output to field required in [ToL ENA sample registration checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000053).
