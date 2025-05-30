# instructpro_data_preprocessing

### Step 1: Downlaod the raw data

Download these files from the following links:
```bash
wget https://static.ramith.io/scientificLLM/data/ligand2smiles.json
wget https://static.ramith.io/scientificLLM/data/protein2ligand_id.json
wget https://static.ramith.io/scientificLLM/data/uniprot2seq.json
wget https://static.ramith.io/scientificLLM/data/uniprot2text.json
wget https://static.ramith.io/scientificLLM/data/uniref50.jsonl
```

### Step 2: Run the Analysis script to see the number of uniprot IDs that have all three: sequence, text, and associated ligands

Running this code will also save the uniprot IDs that have all three: sequence, text, and associated ligands to a file called `selected_uniprot_ids.txt`. This file will be later used for mmseqs2 to generate the clusters.

```bash
(py311) ramith@taurus:~/instructpro_data_preprocessing$ python 1.\ intersection_curation.py 
Protein to Ligand IDs: 12235657
UniProt to Sequence: 33617397
UniProt to Text: 35584366
Ligand, Sequence, and Text intersection: 7680607
Sequence and Text intersection: 33617397
Ligand and Text intersection: 7997705
All elements of seq appear in text
Elements in ligand but not in text: 4237952
Elements in ligand and text but not in seq: 317098
Uniref50 keys: 64363428
Biggest intersection after adding Uniref50: 273568
Regular Intersection: 7680607
Intersection with Uniref50: 273568
Total uniprot ids after adding Uniref50: 7718246
IDs have been saved to selected_uniprot_ids.txt
```