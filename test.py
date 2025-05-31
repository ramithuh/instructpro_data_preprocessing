def load_uniprot_ids(path):
    """
    Read a text file of UniProt IDs, one per line,
    strip whitespace, and return as a set.
    """
    with open(path, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def compare_uniprot_files(file1, file2):
    ids1 = load_uniprot_ids(file1)
    ids2 = load_uniprot_ids(file2)

    print(f"Loaded {len(ids1)} IDs from {file1}")
    print(f"Loaded {len(ids2)} IDs from {file2}")
    print("Comparing UniProt IDs...")
    
    if ids1 == ids2:
        print("Both files contain exactly the same UniProt IDs.")
    else:
        only_in_file1 = ids1 - ids2
        only_in_file2 = ids2 - ids1

        print("Differences found:")
        if only_in_file1:
            print(f"  IDs only in {file1}: {sorted(only_in_file1)}")
        if only_in_file2:
            print(f"  IDs only in {file2}: {sorted(only_in_file2)}")

if __name__ == "__main__":
    
    file1 = "selected_uniprot_ids_taurus.txt"
    file2 = "/mnt/gemini/data/ramith/CMU-project/data/raw/data/selected_uniprot_ids.txt"
    
    var = "ligands_test.txt"
    
    file1 = f"/home/ramith/instructpro_data_preprocessing/{var}"
    file2 = f"/mnt/gemini/data/ramith/CMU-project/data/raw/data/processed/{var}"
    
    var1 = "test_dataset_tokenized.jsonl"
    var2 = "test.dataset_both_unseen_tokenized.jsonl"
    
    file1 = f"/mnt/gemini/data/ramith/CMU-project/data/raw/data/output/{var1}"
    file2 = f"/mnt/gemini/data/ramith/CMU-project/data/splits_with_instruction/has_ligand/{var2}"
    compare_uniprot_files(file1, file2)
