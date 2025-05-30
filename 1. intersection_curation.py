import json
import ijson

def stream_keys(filepath):
    """Generator that yields keys from a top-level JSON object."""
    with open(filepath, 'rb') as f:
        # This iterates over key-value pairs at the root of the JSON object.
        for key, _ in ijson.kvitems(f, ''):
            yield key

def stream_jsonl_keys(filepath):
    """Generator that yields keys from each JSON object in a JSONL file."""
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():  # Skip any empty lines
                try:
                    obj = json.loads(line)
                    if isinstance(obj, dict):
                        for key in obj.keys():
                            yield key
                    else:
                        # If the JSON object is not a dict, skip or handle accordingly.
                        continue
                except json.JSONDecodeError as e:
                    print(f"Error decoding JSON: {e}")
                    continue
                
data_path = "/mnt/gemini/data/ramith/CMU-project/data/raw"

protein2ligand_id  = set(stream_keys(f"{data_path}/data/protein2ligand_id.json"))
uniprot2seq  = set(stream_keys(f"{data_path}/data/uniprot2seq.json"))
uniprot2text = set(stream_keys(f"{data_path}/data/uniprot2text.json"))

print("Protein to Ligand IDs:", len(protein2ligand_id))
print("UniProt to Sequence:", len(uniprot2seq))
print("UniProt to Text:", len(uniprot2text))

# get uniprot ids that have 1) text, 2) associated ligand, and 3) sequence
ligand_seq_text =  protein2ligand_id.intersection(uniprot2seq,uniprot2text)
print("Ligand, Sequence, and Text intersection:", len(ligand_seq_text))

# get intersection of sequence and text
seq_text = uniprot2seq.intersection(uniprot2text)
print("Sequence and Text intersection:", len(seq_text))

# get intersection of ligand and text
ligand_text = protein2ligand_id.intersection(uniprot2text)
print("Ligand and Text intersection:", len(ligand_text))

# Check if all keys in uniprot2seq are also in uniprot2text
diff_seq_text = uniprot2seq.difference(uniprot2text)
if diff_seq_text:
    print(f"Elements in seq not in text: {len(diff_seq_text)}")
else:
    print("All elements of seq appear in text")
    
# Check if all keys in uniprot2seq are also in uniprot2text
diff_ligand_text = protein2ligand_id.difference(uniprot2text)
if diff_ligand_text:
    print(f"Elements in seq not in text: {len(diff_ligand_text)}")
else:
    print("All elements of seq appear in text")


# Count the number of elements in ligand and text but not in seq
result = protein2ligand_id.intersection(uniprot2text).difference(uniprot2seq)
if result:
    print(f"Elements in ligand and text but not in seq: {len(result)}")
else:
    print("All elements of ligand appear in text")
    
# Bring in more sequences from uniref50
uniref50 = set(stream_jsonl_keys(f"{data_path}/data/uniref50.jsonl"))
print("Uniref50 keys:", len(uniref50))

intersection2 = protein2ligand_id.intersection(uniref50, uniprot2text)
print("Biggest intersection after adding Uniref50:", len(intersection2))

print("Regular Intersection:", len(ligand_seq_text))
print("Intersection with Uniref50:", len(intersection2))

## let's calculate the net value of these uniprot ids
total_uniprot_ids = intersection2.union(ligand_seq_text)

print("Total uniprot ids after adding Uniref50:", len(total_uniprot_ids))

# Specify the file path
file_path = 'selected_uniprot_ids.txt'

# Open the file in write mode
with open(file_path, 'w') as file:
    # Write each ID to the file, each on a new line
    for id in total_uniprot_ids:
        file.write(f"{id}\n")

print(f"IDs have been saved to {file_path}")