import sys
import csv
from Bio import Entrez
from pathlib import Path



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



def main(input_file, output_file):
    Entrez.email = "D.parsons@nhm.ac.uk"  #add your email for NCBI API usage

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)  #read input as a dictionary
        fieldnames = reader.fieldnames + ["taxid", "matched_rank"]  #add taxid and matched_rank columns
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        unique_taxids = set()  #store unique TaxIDs

        for row in reader:
            species_name = row.get("Species")
            genus_name = row.get("Genus")
            family_name = row.get("Family")
            order_name = row.get("Order")

            matched_rank = None  #initialize matched rank


#successively go htrough tax ranks and grab taxid
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
            writer.writerow(row)

#write unique TaxIDs to a separate file
    input_filename = Path(input_file).stem  # Get the filename without extension
    unique_taxids_file = f"{input_filename}_unique_taxids.txt"
    with open(unique_taxids_file, 'w') as unique_file:
        for taxid in unique_taxids:
            unique_file.write(f"{taxid}\n")



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 2_sample2taxid.py path/to/BOLD_output.csv sample2taxid.csv")
    else:
        main(sys.argv[1], sys.argv[2])

    print(f"CSV file created successfully.")

