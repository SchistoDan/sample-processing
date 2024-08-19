import csv
import os
from datetime import datetime
from Bio import Entrez
from pathlib import Path
import sys




def load_voucher_data(voucher_tsv_file):
    voucher_dict = {}
    with open(voucher_tsv_file, mode='r') as voucherfile:
        reader = csv.DictReader(voucherfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            museum_id = row.get('Museum ID', 'not collected')
            life_stage = row.get('Life Stage', 'not collected')
            coll_inst = row.get('Institution Storing', 'not collected')

            voucher_dict[sample_id] = {
                'specimen_voucher': museum_id,
                'lifestage': life_stage,
                'collection_ institution': coll_inst
            }
    return voucher_dict




def load_collection_data(collection_tsv_file):
    collection_dict = {}
    with open(collection_tsv_file, mode='r') as collectionfile:
        reader = csv.DictReader(collectionfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            collection_date = row.get('Collection date', 'not collected')
            country_ocean = row.get('Country/Ocean', 'not collected')
            collected_by = row.get('Collectors', 'not collected')
            habitat = row.get('Habitat', 'not collected')
            identified_by = row.get('Identifier', 'not collected')
            exact_site = row.get('region and locality', 'not collected')
            latitude = row.get('Lat', 'not collected')
            longitude = row.get('Lon', 'not collected')

#Convert date to BCDM format (YYYY-MM-DD)
            try:
                collection_date = datetime.strptime(collection_date, "%d-%m-%Y").strftime("%Y-%m-%d")
            except ValueError:
#If dates are malformed then they are not currently handled gracefully
                collection_date = "not collected"  

            collection_dict[sample_id] = {
                'collection date': collection_date,
                'geographic_location': country_ocean,
                'geographic_location_locality': exact_site,
                'latitude': latitude,
                'longitude': longitude,
                'collected_by': collected_by,
                'habitat': habitat,
                'identified_by': identified_by
            }
    return collection_dict





def load_specimen_details(specimen_tsv_file):
    specimen_dict = {}
    with open(specimen_tsv_file, mode='r') as specimenfile:
        reader = csv.DictReader(specimenfile, delimiter='\t')
        for row in reader:
            sample_id = row['Sample ID']
            notes = row.get('Notes', 'not collected')
            sex = row.get('Sex', 'not collected')

            specimen_dict[sample_id] = {
                'organism part': notes,
                'sex': sex
            }
    return specimen_dict





def locate_files_in_directory(directory_path):
    files = os.listdir(directory_path)
    voucher_file = collection_file = specimen_file = None

    for file in files:
        if file.endswith("voucher.tsv"):
            voucher_file = os.path.join(directory_path, file)
        elif file.endswith("collection_data.tsv"):
            collection_file = os.path.join(directory_path, file)
        elif file.endswith("specimen_details.tsv"):
            specimen_file = os.path.join(directory_path, file)

    if not voucher_file or not collection_file or not specimen_file:
        raise FileNotFoundError("One or more required files (voucher.tsv, collection_data.tsv, specimen_details.tsv) are missing in the specified directory.")

    return voucher_file, collection_file, specimen_file





def fetch_taxid(taxon_name, taxonomic_rank):
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
        print(f"Error fetching TaxID for {taxon_name} ({taxonomic_rank}): {e}")
        return None




def main(input_file, directory_path, output_file):
    Entrez.email = "user@email.ac.uk"  #Add your email for NCBI API usage

#Locate BOLD files in provided directory & load
    voucher_tsv_file, collection_tsv_file, specimen_tsv_file = locate_files_in_directory(directory_path)

    voucher_dict = load_voucher_data(voucher_tsv_file)
    collection_dict = load_collection_data(collection_tsv_file)
    specimen_dict = load_specimen_details(specimen_tsv_file)

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + [
            "taxid", "matched_rank", "specimen_voucher", "lifestage", "collection date",
            "geographic_location (country and/or sea)", "geographic_location (region and locality)",
            "latitude", "longitude",
            "collected_by", "habitat", "identified_by", "collection institution",
            "organism part", "sex"
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        unique_taxids = set()  # Store unique TaxIDs

        for row in reader:
            sample_id = row.get("Sample ID")
            species_name = row.get("Species")
            genus_name = row.get("Genus")
            family_name = row.get("Family")
            order_name = row.get("Order")

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
            row["lifestage"] = voucher_info.get('lifestage', "not collected")
            row["collection institution"] = voucher_info.get('collection institution', "not collected")


#Retrieve and format collection data from collection_dict
            collection_info = collection_dict.get(sample_id, {})
            row["collection date"] = collection_info.get('collection date', "not collected")
            row["geographic_location (country and/or sea)"] = collection_info.get('geographic_location', "not collected")
            row["geographic_location (region and locality)"] = collection_info.get('geographic_location_locality', "not collected")
            row["latitude"] = collection_info.get('latitude', "not collected")
            row["longitude"] = collection_info.get('longitude', "not collected")
            row["collected_by"] = collection_info.get('collected_by', "not collected")
            row["habitat"] = collection_info.get('habitat', "not collected")
            row["identified_by"] = collection_info.get('identified_by', "not collected")


#Retrieve organism part and sex from specimen_dict
            specimen_info = specimen_dict.get(sample_id, {})
            row["organism part"] = specimen_info.get('organism part', "not collected")
            row["sex"] = specimen_info.get('sex', "not collected")

            writer.writerow(row)


#Write unique TaxIDs to a separate file
    input_filename = Path(input_file).stem  # Get the filename without extension
    unique_taxids_file = f"{input_filename}_unique_taxids.txt"
    with open(unique_taxids_file, 'w') as unique_file:
        for taxid in unique_taxids:
            unique_file.write(f"{taxid}\n")

    print(f"CSV file created successfully with data from voucher, collection, and specimen tsv files.")




if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""
        Usage: python 2_sample2taxid.py path/to/BOLD_output.csv path/to/directory_containing_tsv_files [output].csv
        """)
        sys.exit(1)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
