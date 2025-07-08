import csv
import json
from typing import List, Tuple, Set

def load_allow_list(filepath: str) -> Set[Tuple[str, str]]:
    """
    Loads the (protein_id, ligand_id) tuples from a TSV file into a set
    for very fast lookups.

    Args:
        filepath: The path to the allow list file (e.g., 'allow_list.tsv').

    Returns:
        A set of tuples, e.g., {('P12345', 'CHEBI_123'), ...}
    """
    print(f"Loading allow list from {filepath}...")
    allowed_set = set()
    try:
        with open(filepath, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) # Skip the header row
            for row in reader:
                if len(row) == 2:
                    allowed_set.add(tuple(row))
    except FileNotFoundError:
        print(f"Error: Allow list file not found at {filepath}")
    except Exception as e:
        print(f"Error reading allow list file: {e}")
    
    print(f"Loaded {len(allowed_set)} items into the allow list set.")
    return allowed_set

def filter_jsonl_by_allow_list(input_jsonl: str, output_jsonl: str, allow_list_path: str):
    """
    Reads a large JSONL file and writes a new, downsampled version based on an allow list.

    Args:
        input_jsonl: Path to the original large .jsonl file.
        output_jsonl: Path where the new downsampled .jsonl file will be saved.
        allow_list_path: Path to the .tsv file containing the items to keep.
    """
    # Load the IDs into a set for O(1) average time complexity lookups
    allowed_set = load_allow_list(allow_list_path)
    
    if not allowed_set:
        print("Allow list is empty. No output file will be generated.")
        return

    print(f"Starting to filter {input_jsonl}...")
    lines_kept = 0
    total_lines = 0

    try:
        with open(input_jsonl, 'r') as infile, open(output_jsonl, 'w') as outfile:
            for line in infile:
                total_lines += 1
                try:
                    data = json.loads(line)
                    protein_id = data.get('uniprot_id')
                    ligand_id = data.get('input', {}).get('ligand_id')

                    # Check if the (protein, ligand) pair is in our allowed set
                    if (protein_id, ligand_id) in allowed_set:
                        # If it is, write the original line to the new file
                        outfile.write(line)
                        lines_kept += 1
                except (json.JSONDecodeError, KeyError):
                    # Skip malformed lines
                    continue
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_jsonl}")
        return
        
    print("\n--- Filtering Complete ---")
    print(f"Total lines read: {total_lines:,}")
    print(f"Lines kept: {lines_kept:,}")
    print(f"New downsampled file saved to: {output_jsonl}")


if __name__ == '__main__':
    # --- This is the main execution block for your NEW filtering script ---
    
    # 1. Define your file paths
    # The original, large file you want to filter
    original_file = '/mnt/gemini/data/ramith/CMU-project/data/splits_with_instruction/has_ligand/test.dataset_both_unseen_tokenized.jsonl' 
    
    # The allow list you generated with your other script
    allow_list_file = 'allow_list.tsv' 
    
    # The name for your final, small, downsampled output file
    output_file = 'test_unseen_downsampled.jsonl'

    # 2. Run the filtering process
    filter_jsonl_by_allow_list(
        input_jsonl=original_file,
        output_jsonl=output_file,
        allow_list_path=allow_list_file
    )