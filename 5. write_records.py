import ijson
import json
import csv
from collections import defaultdict
from typing import Set, Dict, Optional, Generator, Tuple
from transformers import AutoTokenizer

data_path = "/mnt/gemini/data/ramith/CMU-project/data/raw"
output_path = "/mnt/gemini/data/ramith/CMU-project/data/raw/data/output"

# --- Tokenizer Initialization (Outside any function) ---
tokenizer = AutoTokenizer.from_pretrained("microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext")
special_tokens_dict = {
    'additional_special_tokens': ["<FUNCTION>", "</FUNCTION>"]
}
tokenizer.add_special_tokens(special_tokens_dict)
# --- (Rest of your import statements) ---

def load_json_file(filepath: str) -> dict:
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")

def load_cluster_ids(filepath: str) -> Set[str]:
    try:
        with open(filepath, 'r') as f:
            return set(line.strip() for line in f)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")

def load_ligand_ids(filepath: str) -> Set[str]:
    try:
        with open(filepath, 'r') as f:
            return set(line.strip() for line in f)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")

def get_uniprot_ids_from_clusters(cluster_ids: Set[str], cluster_mapping_file: str) -> Set[str]:
    uniprot_ids = set()
    try:
        with open(cluster_mapping_file, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:
                if row and row[0] in cluster_ids:
                    uniprot_ids.add(row[1])
    except FileNotFoundError:
        raise FileNotFoundError(f"Cluster mapping file not found: {cluster_mapping_file}")
    return uniprot_ids

def partial_load(json_path: str, needed_uniprot_ids: Set[str]) -> Dict[str, any]:
    result = {}
    try:
        with open(json_path, 'rb') as f:
            for key, value in ijson.kvitems(f, ''):
                if key in needed_uniprot_ids:
                    result[key] = value
    except FileNotFoundError:
        raise FileNotFoundError(f"JSON file not found: {json_path}")
    except Exception as e:  # Catch other potential parsing errors
        print(f"Error during partial loading of {json_path}: {e}")
    return result

def strip_sdf(ligand_id: str) -> str:
    if ligand_id.endswith('.sdf'):
        return ligand_id[:-4]
    return ligand_id

def tokenize_text(text: str, tokenizer: AutoTokenizer) -> Dict[str, any]: #tokenizer added
    """
    Tokenizes the input text using the provided tokenizer.

    Args:
        text: The text to tokenize.
        tokenizer: The tokenizer to be used

    Returns:
        A dictionary containing the original text, tokenized text (as strings),
        and token IDs (as numbers).
    """
    prefixed_text = f"The function of the target protein is <FUNCTION> {text} </FUNCTION>" #prefix added
    tokens = tokenizer.tokenize(prefixed_text)
    token_ids = tokenizer.encode(prefixed_text, add_special_tokens=True)
    return {
        "original_text": text,
        "tokenized_text": tokens,
        "token_ids": token_ids,
    }

def generate_records(protein2ligand_dict: Dict[str, list], ligand2smiles_dict: Dict[str, str],
                     uniprot2seq_dict: Dict[str, str], uniprot2text_dict: Dict[str, str],
                     allowed_uniprot_ids: Set[str], allowed_ligand_ids: Optional[Set[str]] = None) -> Generator[dict, None, None]:
    for uniprot_id, ligand_list in protein2ligand_dict.items():
        if uniprot_id not in allowed_uniprot_ids:
            continue

        seq = uniprot2seq_dict.get(uniprot_id, None)
        text = uniprot2text_dict.get(uniprot_id, None)
        if not seq or not text:
            continue

        # --- Tokenization ---
        tokenized_data = tokenize_text(text, tokenizer) #tokenizer added

        for lig_id in ligand_list:
            clean_lig_id = strip_sdf(lig_id)
            if allowed_ligand_ids is not None and clean_lig_id not in allowed_ligand_ids:
                continue

            smiles = ligand2smiles_dict.get(clean_lig_id, None)
            if smiles is None:
                continue

            yield {
                "input": {
                    "function_original": tokenized_data["original_text"], #original text
                    "function_tokens":   tokenized_data["tokenized_text"], #token strings
                    "function_token_ids":tokenized_data["token_ids"], #token ids
                    "ligand": smiles,
                    "ligand_id": clean_lig_id
                },
                "output": seq,
                "uniprot_id": uniprot_id
            }

def write_dataset_jsonl(record_generator: Generator[dict, None, None], out_path: str):
    try:
        with open(out_path, 'w') as f:
            for rec in record_generator:
                f.write(json.dumps(rec))
                f.write('\n')
    except Exception as e:
        print(f"Error writing to {out_path}: {e}")

# --- Main Script Execution --- 
# 1. Load Data and Cluster/Ligand IDs
cluster_mapping_file = f"{data_path}/data/clusterRes_cluster.tsv"  # Path to your cluster mapping file

train_cluster_ids = load_cluster_ids("train_clusters.txt")
val_cluster_ids = load_cluster_ids("val_clusters.txt")
test_cluster_ids = load_cluster_ids("test_clusters.txt")

val_ligand_ids = load_ligand_ids("ligands_val.txt")
test_ligand_ids = load_ligand_ids("ligands_test.txt")

# Get UniProt IDs from the cluster IDs
train_uniprot_ids = get_uniprot_ids_from_clusters(train_cluster_ids, cluster_mapping_file)
val_uniprot_ids = get_uniprot_ids_from_clusters(val_cluster_ids, cluster_mapping_file)
test_uniprot_ids = get_uniprot_ids_from_clusters(test_cluster_ids, cluster_mapping_file)

protein2ligand = load_json_file(f"{data_path}/data/protein2ligand_id.json")
ligand2smiles = load_json_file(f"{data_path}/data/ligand2smiles.json")
# print(f"Loaded ligand2smiles of {len(ligand2smiles):,} ligands.") No need to print this every time.

# Get *all* needed uniprot IDs, for loading sequence/text data
all_needed_uniprot_ids = train_uniprot_ids | val_uniprot_ids | test_uniprot_ids

uniprot2seq  = partial_load(f"{data_path}/data/uniprot2seq.json", all_needed_uniprot_ids)  # Efficient loading!
uniprot2text = partial_load(f"{data_path}/data/uniprot2text.json", all_needed_uniprot_ids)


# ### Save test and validation splits that have seen ligands during training

train_allowed_ligands = set(ligand2smiles.keys()) - val_ligand_ids - test_ligand_ids  #Correct way to get allowed training ligands

# Restrict protein2ligand to just the validation proteins
val_protein2ligand = {k: v for k, v in protein2ligand.items() if k in val_uniprot_ids}

print(f"total unique protein_ids {len(val_uniprot_ids)}")

val_record_gen = generate_records(
    val_protein2ligand,
    ligand2smiles,
    uniprot2seq,
    uniprot2text,
    val_uniprot_ids,
    train_allowed_ligands   # here we want the ligands to be stuff we saw during training
)
write_dataset_jsonl(val_record_gen, f"{output_path}/val_dataset_seen_ligands_tokenized.jsonl") #file name changed
print("Done writing validation dataset [seen ligands].")


# Restrict protein2ligand
test_protein2ligand = {k: v for k, v in protein2ligand.items() if k in test_uniprot_ids}

print(f"total unique protein_ids {len(test_uniprot_ids)}")

test_record_gen = generate_records(
    test_protein2ligand,
    ligand2smiles,
    uniprot2seq,
    uniprot2text,
    test_uniprot_ids,
    train_allowed_ligands # here we want the ligands to be stuff we saw during training
)
write_dataset_jsonl(test_record_gen, f"{output_path}/test.dataset_seen_ligands_tokenized.jsonl") #file name changed
print("Done writing test dataset [seen ligands].")


# ### Save unseen ligands in train for val and test set 

# --- Validation Data ---
# Restrict protein2ligand to just the validation proteins
val_protein2ligand = {k: v for k, v in protein2ligand.items() if k in val_uniprot_ids}

val_record_gen = generate_records(
    val_protein2ligand,
    ligand2smiles,
    uniprot2seq,
    uniprot2text,
    val_uniprot_ids,
    val_ligand_ids   # Only allow validation ligands
)
write_dataset_jsonl(val_record_gen, f"{output_path}/val_dataset_tokenized.jsonl") #file name changed
print("Done writing validation dataset.")


# --- Test Data ---
# Restrict protein2ligand
test_protein2ligand = {k: v for k, v in protein2ligand.items() if k in test_uniprot_ids}

test_record_gen = generate_records(
    test_protein2ligand,
    ligand2smiles,
    uniprot2seq,
    uniprot2text,
    test_uniprot_ids,
    test_ligand_ids  # Only allow test ligands
)
write_dataset_jsonl(test_record_gen, f"{output_path}/test.dataset_both_unseen_tokenized.jsonl") #file name changed
print("Done writing test dataset.")

# 2. Generate and Write Datasets

# --- Training Data ---
#  * Exclude val and test ligands
train_allowed_ligands = set(ligand2smiles.keys()) - val_ligand_ids - test_ligand_ids  #Correct way to get allowed training ligands
# Restrict protein2ligand to just the training proteins we care about
train_protein2ligand = {k: v for k, v in protein2ligand.items() if k in train_uniprot_ids}

train_record_gen = generate_records(
    train_protein2ligand,
    ligand2smiles,
    uniprot2seq,
    uniprot2text,
    train_uniprot_ids,
    train_allowed_ligands
)
write_dataset_jsonl(train_record_gen, f"{output_path}/train_dataset_tokenized.jsonl") #file name changed
print("Done writing training dataset.")


len(train_allowed_ligands)