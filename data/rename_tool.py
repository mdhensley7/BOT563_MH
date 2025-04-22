#This tool will change the label of a sequence from the accession number to the species name
#Be sure to check to make sure the names of the .csv file are correct and labelled accordingly
#Place this code in the directory where the .csv file is located
import csv
from Bio import SeqIO

# Paths to your files
csv_file = "/Users/hensley/Desktop/BOT563/BOT563_MH/data/3loci_24samp_working/am3_loci.csv"
alignment_file = "/Users/hensley/Desktop/BOT563/BOT563_MH/data/3loci_24samp_working/fasta/tef1a_clustal_trim_9.fasta"  # Replace with your alignment file path
output_file = "/Users/hensley/Desktop/BOT563/BOT563_MH/data/3loci_24samp_working/alignments/tef1a/tef1a_clustal_spname_9.fasta"  # Output file with renamed sequences

# Step 1: Create a mapping of accession numbers to species names
accession_to_species = {}

with open(csv_file, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for column in reader:
        species_name = column["species"].strip()
        tef1a_accession = column["tef1a"].strip()
        accession_to_species[tef1a_accession] = species_name

# Step 2: Rename sequences in the alignment file
with open(alignment_file, "r") as infile, open(output_file, "w") as outfile:
    records = SeqIO.parse(infile, "fasta")
    renamed_records = []
    for record in records:
        if record.id in accession_to_species:
            # Rename the sequence header to the species name
            record.id = accession_to_species[record.id]
            record.description = ""  # Clear the description
        renamed_records.append(record)
    # Write the renamed sequences to the output file
    SeqIO.write(renamed_records, outfile, "fasta")

print(f"Renamed sequences saved to {output_file}")