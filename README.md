Python scripts that take BOLD container download (.tsv files) and paths to raw reads, fetches taxonomic ID (taxid) from local NCBI tax dump using taxonomic ranks for each sample, merges relevant fields into output ([output].csv) from BOLD-downlaoded files needed for downstream [upload of reads to ENA](https://github.com/bge-barcoding/ena-read-upload), creates analysis pipeline sample submission CSV (samples.csv) for input into [skim2mito](https://github.com/o-william-white/skim2mito) and [MGE](https://github.com/bge-barcoding/MitoGeneExtractor-BGE), and requesting taxid creation from ENA.



## 1_sample_processing.py ##
Script compiles fields requried for downstream uses (e.g. submission of reads to ENA) from .tsv files downloaded from BOLD. 
**usage: python sample_processing.py [path/to/BOLD/download/dir]**
- path/to/BOLD/download/dir: Path to directory containing .tsv files (voucher.tsv, collection_data.tsv, specimen_details.tsv, taxonomy.tsv, and lab.tsv)
- outputs 'sample_metadata.csv' (see below for example).


## 2_samples_spreadsheet.py ##
Script to generate samples.csv requried to run MGE and skim2mito piplines.
**usage: python 2_samples_spreadsheet.py [path/to/raw/read/dir] [path/to/samples_metadata.csv]**
- path/to/raw/read/dir: Path to parent directory of subdirectories containing raw PE read files.
- path/to/output/dir/samples_metadata.csv: Directory containing sample_metadta.csv file.
- outputs 'samples.csv' (see below for example).



**Example samples.csv**
| ID  | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM001-24  | abs/path/to/R1.fq.gz  | abs/path/to/R2.fq.gz | 1234 |
| BSNHM001-24 | abs/path/to/R1.fq.gz  | abs/path/to/R2.fq.gz | 5678 |
| BSNHM001-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz |  91011 |



**Example samples_metadata.csv**
| Sample ID | Process ID  | Phylum | Class | Order | Family | Subfamily | Tribe | Genus | Species | Subspecies | taxid | matched_rank | specimen_voucher | lifestage | collection date | geographic_location (country and/or sea) | geographic_location (region and locality) | latitude | longitude | collected_by | habitat | identified_by | collection institution | organism part | sex |
| --- | --- |--- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |--- | --- | --- | --- | --- | --- | --- | --- | ---| --- | --- | --- |
| BGE_0001_A01  | BSNHM001-24 | Arthropoda | Insecta | Trichoptera | Apataniidae | Apataniinae | | Apatania | Apatania stylata | | 177658 | genus | 'Museum ID' | adult | YYYY-MM-DD | France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Whole | M |
| BGE_0001_A02 | BSNHM002-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Agapetinae | | Agapetus | Agapetus iridipennis | | 177627 | genus | 'Museum ID' | adult | YYYY-MM-DD | Switzerland | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | Naturalis | not collected | F |
| BGE_0001_A03 | BSNHM003-24 | Arthropoda | Insecta | Trichoptera | Glossosomatidae | Hydropsychidae | | Diplectrona | Diplectrona meridionalis | | 177860 | genus | 'Museum ID' | adult | YYYY-MM-DD |  France | not collected | lat (DD) | lon (DD) | 'Collectors' | not collected | not collected | NHMMUK | Leg | M |

- If data not collected for sample, 'not collected' output to field required in [ToL ENA sample registration checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000053).
