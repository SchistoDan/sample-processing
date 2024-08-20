import csv
import os
from datetime import datetime
from Bio import Entrez
from pathlib import Path
import sys
import time




def load_voucher_data(voucher_tsv_file):
    voucher_dict = {}
    with open(voucher_tsv_file, mode='r') as voucherfile:
        reader = csv.DictReader(voucherfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            museum_id = row.get('Museum ID', 'not collected')
            coll_inst = row.get('Institution Storing', 'not collected')

            voucher_dict[sample_id] = {
                'specimen_voucher': museum_id,
                'collecting institution': coll_inst
            }
    return voucher_dict




def load_collection_data(collection_tsv_file):
    collection_dict = {}
    with open(collection_tsv_file, mode='r') as collectionfile:
        reader = csv.DictReader(collectionfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            collection_date = row.get('Collection Date', 'not collected')
            country_ocean = row.get('Country/Ocean', 'not collected')
            collected_by = row.get('Collectors', 'not collected')
            habitat = row.get('Habitat', 'not collected')
            exact_site = row.get('Exact Site', 'not collected')
            latitude = row.get('Lat', 'not collected')
            longitude = row.get('Lon', 'not collected')

#Convert dates to BCDM format (YYYY-MM-DD)
            if collection_date != 'not collected':
                try:
                    collection_date = datetime.strptime(collection_date, "%d-%b-%Y").strftime("%Y-%m-%d")
                except ValueError:
                    try:
                        collection_date = datetime.strptime(collection_date, "%d-%b-%y").strftime("%Y-%m-%d")
                    except ValueError:
                        collection_date = "not collected"
            else:
                collection_date = "not collected"


#Replace 'Unrecoverable', empty, or missing values in fields with 'not collected'
            if not country_ocean or country_ocean == 'Unrecoverable':
                country_ocean = 'not collected'

            if not exact_site:
                exact_site = 'not collected'

            if not latitude:
                latitude = 'not collected'

            if not longitude:
                longitude = 'not collected'

            if not collected_by:
                collected_by = 'not collected'

            if not habitat:
                habitat = 'not collected'


            collection_dict[sample_id] = {
                'collection date': collection_date,
                'geographic_location': country_ocean,
                'geographic_location_locality': exact_site,
                'latitude': latitude,
                'longitude': longitude,
                'collected_by': collected_by,
                'habitat': habitat,
            }
    return collection_dict




def load_specimen_details(specimen_tsv_file):
    specimen_dict = {}
    with open(specimen_tsv_file, mode='r') as specimenfile:
        reader = csv.DictReader(specimenfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            organism_part = row.get('Tissue Descriptor', 'not collected')
            sex = row.get('Sex', 'not collected')
            life_stage = row.get('Life Stage', 'not collected')


#Replace empty or missing values in fields with 'not collected'
            if not organism_part:
                 organism_part = 'not collected'

            if not life_stage:
                 life_stage = 'not collected'

            if not sex:
                 sex = 'not collected'


            specimen_dict[sample_id] = {
                'organism part': organism_part,
                'sex': sex,
                'lifestage': life_stage
            }
    return specimen_dict




def load_taxonomy_data(taxonomy_tsv_file):
    taxonomy_dict = {}
    with open(taxonomy_tsv_file, mode='r') as taxonomyfile:
        reader = csv.DictReader(taxonomyfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            identifier = row.get('Identifier', 'not collected')
            taxonomy_dict[sample_id] = {
                'identified_by': identifier
            }
    return taxonomy_dict





def locate_files_in_directory(directory_path):
    files = os.listdir(directory_path)
    voucher_file = collection_file = specimen_file = taxonomy_file = None


    for file in files:
        if file.endswith("voucher.tsv"):
            voucher_file = os.path.join(directory_path, file)
        elif file.endswith("collection_data.tsv"):
            collection_file = os.path.join(directory_path, file)
        elif file.endswith("specimen_details.tsv"):
            specimen_file = os.path.join(directory_path, file)
        elif file.endswith("taxonomy.tsv"):
            taxonomy_file = os.path.join(directory_path, file)


    if not voucher_file or not collection_file or not specimen_file or not taxonomy_file:
        missing_files = []
        if not voucher_file:
            missing_files.append("voucher.tsv")
        if not collection_file:
            missing_files.append("collection_data.tsv")
        if not specimen_file:
            missing_files.append("specimen_details.tsv")
        if not taxonomy_file:
            missing_files.append("taxonomy.tsv")
        
        raise FileNotFoundError(f"Missing required files: {', '.join(missing_files)}")

    return voucher_file, collection_file, specimen_file, taxonomy_file





def fetch_taxid(taxon_name, taxonomic_rank, retries=5, delay=3):
    attempt = 0
    while attempt < retries:
        try:
            term = f"{taxon_name}[{taxonomic_rank}]"
            handle = Entrez.esearch(db="taxonomy", term=term)
            record = Entrez.read(handle)
            handle.close()
            if record['IdList']:
                return record['IdList'][0]
            else:
                return None
        except Exception as e:
            print(f"Error fetching TaxID for {taxon_name} ({taxonomic_rank}) on attempt {attempt + 1}: {e}")
            attempt += 1
            time.sleep(delay)  

#If all retries fail, return None
    print(f"Failed to fetch TaxID for {taxon_name} ({taxonomic_rank}) after {retries} attempts.")
    return None






def main(input_file, directory_path, output_file):
    Entrez.email = "d.parsons@nhm.ac.uk"  # Add your email for NCBI API usage

#Locate BOLD files in provided directory & load
    voucher_tsv_file, collection_tsv_file, specimen_tsv_file, taxonomy_tsv_file = locate_files_in_directory(directory_path)

    print(f"Located files: voucher_tsv_file={voucher_tsv_file}, collection_tsv_file={collection_tsv_file}, specimen_tsv_file={specimen_tsv_file}, taxonomy_tsv_file={taxonomy_tsv_file}")


    voucher_dict = load_voucher_data(voucher_tsv_file)
    collection_dict = load_collection_data(collection_tsv_file)
    specimen_dict = load_specimen_details(specimen_tsv_file)
    taxonomy_dict = load_taxonomy_data(taxonomy_tsv_file)
    

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + [
            "taxid", "matched_rank", "specimen_voucher", "lifestage", "collection date",
            "geographic location (country and/or sea)", "geographic location (region and locality)",
            "geographic location (latitude)", "geographic location (longitude)",
            "collected_by", "habitat", "identified_by", "collecting institution",
            "organism part", "sex"
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        unique_taxids = set()  #Store unique TaxIDs. May be useful downstream


        for row in reader:
            sample_id = row.get("Sample ID")
            species_name = row.get("Species")
            genus_name = row.get("Genus")
            family_name = row.get("Family")
            order_name = row.get("Order")

            print(f"Processing sample_id={sample_id}, species_name={species_name}, genus_name={genus_name}, family_name={family_name}, order_name={order_name}")

            matched_rank = None


#Successively go through taxonomic ranks and grab taxid
            if species_name:
                species_taxid = fetch_taxid(species_name, 'species')
                if species_taxid:
                    row["taxid"] = species_taxid
                    matched_rank = "species"
                    unique_taxids.add(species_taxid)

            if not matched_rank and genus_name:
                genus_taxid = fetch_taxid(genus_name, 'genus')
                if genus_taxid:
                    row["taxid"] = genus_taxid
                    matched_rank = "genus"
                    unique_taxids.add(genus_taxid)

            if not matched_rank and family_name:
                family_taxid = fetch_taxid(family_name, 'family')
                if family_taxid:
                    row["taxid"] = family_taxid
                    matched_rank = "family"
                    unique_taxids.add(family_taxid)

            if not matched_rank and order_name:
                order_taxid = fetch_taxid(order_name, 'order')
                if order_taxid:
                    row["taxid"] = order_taxid
                    matched_rank = "order"
                    unique_taxids.add(order_taxid)

            row["matched_rank"] = matched_rank


#Retrieve specimen_voucher, lifestage, and collection institution from voucher_dict
            voucher_info = voucher_dict.get(sample_id, {})
            row["specimen_voucher"] = voucher_info.get('specimen_voucher', "not collected")
            row["collecting institution"] = voucher_info.get('collecting institution', "not collected")


#Retrieve and format collection data from collection_dict
            collection_info = collection_dict.get(sample_id, {})
            row["collection date"] = collection_info.get('collection date', "not collected")
            row["geographic location (country and/or sea)"] = collection_info.get('geographic_location', "not collected")
            row["geographic location (region and locality)"] = collection_info.get('geographic_location_locality', "not collected")
            row["geographic location (latitude)"] = collection_info.get('latitude', "not collected")
            row["geographic location (longitude)"] = collection_info.get('longitude', "not collected")
            row["collected_by"] = collection_info.get('collected_by', "not collected")
            row["habitat"] = collection_info.get('habitat', "not collected")
            row["identified_by"] = collection_info.get('identified_by', "not collected")


#Retrieve organism part and sex from specimen_dict
            specimen_info = specimen_dict.get(sample_id, {})
            row["organism part"] = specimen_info.get('organism part', "not collected")
            row["sex"] = specimen_info.get('sex', "not collected")
            row["lifestage"] = specimen_info.get('lifestage', "not collected")


#Retrieve identified_by from taxonomy_dict
            taxonomy_info = taxonomy_dict.get(sample_id, {})
            row["identified_by"] = taxonomy_info.get('identified_by', "not collected")


            writer.writerow(row)
            print(f"Wrote row for sample_id={sample_id}")


#Write unique TaxIDs to a separate file
    input_filename = Path(input_file).stem  # Get the filename without extension
    unique_taxids_file = f"{input_filename}_unique_taxids.txt"
    with open(unique_taxids_file, 'w') as unique_file:
        for taxid in unique_taxids:
            unique_file.write(f"{taxid}\n")

    print(f"CSV file created successfully with data from voucher, collection, specimen, and taxonomy TSV files.")




if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""
        Usage: python your_script.py path/to/BOLD_output.csv path/to/directory_containing_BOLd_tsv_files [output].csv
        """)
        sys.exit(1)

    input_file = sys.argv[1]
    directory_path = sys.argv[2]
    output_file = sys.argv[3]

    main(input_file, directory_path, output_file)

