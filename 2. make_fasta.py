import json
import ijson

# Load selected UniProt IDs into a set
with open('selected_uniprot_ids.txt', 'r') as id_file:
    selected_ids = set(line.strip() for line in id_file)

data_path = "/mnt/gemini/data/ramith/CMU-project/data/raw"

# Open the output FASTA file in write mode
with open(f"output_sequences.fasta", 'w') as fasta_out:
    # First, stream through the main uniprot2seq.json file (13 GB)

    c1 = 0
    with open(f"{data_path}/data/uniprot2seq.json", 'r') as json_file:
        # Assumes the JSON structure is a dictionary mapping IDs to sequences
        # Adjust the prefix in kvitems() if your structure is different.
        for uniprot_id, sequence in ijson.kvitems(json_file, ''):
            if uniprot_id in selected_ids:
                fasta_out.write(f">{uniprot_id}\n{sequence}\n")
                selected_ids.remove(uniprot_id)  # Remove to avoid duplicate processing
                c1+=1

    print(f"{c1} written from uniprot2seq.json")

    c2 = 0
    # Then process the uniref50.jsonl file (19 GB) for any IDs not found above
    with open(f"{data_path}/data/uniref50.jsonl", 'r') as jsonl_file:
        for line in jsonl_file:
            record = json.loads(line)
            # Adjust these keys according to your JSONL file's schema
            uniprot_id, sequence = next(iter(record.items()))
            
            if uniprot_id in selected_ids:
                fasta_out.write(f">{uniprot_id}\n{sequence}\n")
                selected_ids.remove(uniprot_id)
                c2+=1

    print(f"{c2} written from uniref50.jsonl")

print(len(selected_ids))