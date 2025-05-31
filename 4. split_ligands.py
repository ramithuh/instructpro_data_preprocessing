import random
random.seed(42)
import json
import numpy as np
from helpers import * 
import plotly.io as pio
pio.renderers.default = "colab" 

data_path = "/mnt/gemini/data/ramith/CMU-project/data/raw"

def load_json_file(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

protein2ligand = load_json_file(f"{data_path}/data/protein2ligand_id.json")

# Load selected UniProt IDs into a set
with open(f"{data_path}/data/selected_uniprot_ids.txt", 'r') as id_file:
    selected_ids = set(line.strip() for line in id_file)


# # Curated 7.7 Million protein analysis
# 

# ## Counting number of Protein IDs

# ### 💡 Number of unique protein ids 

# number of unique uniprot_ids 
needed_uniprot_ids = selected_ids.intersection(set(protein2ligand.keys()))
print(f"We selected {len(needed_uniprot_ids):,} proteins that have both text and seqs")

filtered_protein2ligand = {}

for uniprot_id in protein2ligand.keys():
    if(uniprot_id in selected_ids):
        filtered_protein2ligand[uniprot_id] = protein2ligand[uniprot_id]

del protein2ligand

print(list(filtered_protein2ligand.keys())[0])

print(filtered_protein2ligand['O33839'], filtered_protein2ligand['P32234'])


# ## Counting number of ligands

all_ligands = []
total_pairings = 0

protein_counts_for_ligand = {}

for idx, current_protein in enumerate(filtered_protein2ligand.keys()):
    # extend list
    ligands_for_curr_protein = filtered_protein2ligand[current_protein]

    for ligand in ligands_for_curr_protein: #loop through each ligand
        if ligand not in protein_counts_for_ligand:
            protein_counts_for_ligand[ligand] = 1
        else:
            protein_counts_for_ligand[ligand] += 1
    
    all_ligands.extend(ligands_for_curr_protein)
    total_pairings += len(ligands_for_curr_protein)


num_unique_ligands = len(set(all_ligands))

assert total_pairings == len(all_ligands)


# ### 💡 Total uniprot_id <---> ligand pairings

len(all_ligands) #total pairings


# ### 💡 number of unique ligands

num_unique_ligands


# Sort by count in descending order
sorted_ligands = sorted(protein_counts_for_ligand.items(), key=lambda x: x[1], reverse=True)


for idx, (ligand, count) in enumerate(sorted_ligands):

    if(idx % 50 == 0):
        print(ligand, count)


SAMPLING_RANGES = []

increment = 10000
for i in range(1,50000,increment):
    SAMPLING_RANGES.append(((i,i+increment),4))


print(SAMPLING_RANGES)

print(sorted_ligands[0])

val_selected_ligands  = []
test_selected_ligands = []

for (low, high), desired_count in SAMPLING_RANGES:
    # Filter clusters that fall in the [low, high) range
    matching = [(l_name, cl) for (l_name, cl) in sorted_ligands if cl>=low and cl<high]

    # override : take all
    desired_count = len(matching)
    
    if len(matching) < desired_count:
        print(f"Warning: asked for {desired_count} clusters in range [{low}, {high}), "
              f"but only {len(matching)} available. therefore..")

        if(len(matching) < desired_count - 2):
            continue
        else:
            print(f" ==> but we can sample {desired_count - 2}")
            desired_count = desired_count - 2

    
    chosen = random.sample(matching, desired_count)
    
    val_selected_ligands.extend(chosen)#[:desired_count//2])
    test_selected_ligands.extend(chosen)#[desired_count//2:])



len(test_selected_ligands), len(val_selected_ligands), 


cluster_mapping_file = f"{data_path}/data/clusterRes_cluster.tsv"

train_cluster = load_cluster_split("train_clusters.txt")
val_cluster   = load_cluster_split("val_clusters.txt")
test_cluster  = load_cluster_split("test_clusters.txt")


# train_cluster_details , train_cluster_counts  = get_uniprot_ids_of_cluster(train_cluster , cluster_mapping_file)
val_cluster_details , val_cluster_counts  = get_uniprot_ids_of_cluster(val_cluster , cluster_mapping_file)
test_cluster_details, test_cluster_counts = get_uniprot_ids_of_cluster(test_cluster, cluster_mapping_file)


np.array([counts for _, counts in val_cluster_counts]).sum()
np.array([counts for _, counts in test_cluster_counts]).sum()


def load_json_file(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

protein2ligand = load_json_file(f"{data_path}/data/protein2ligand_id.json")

def get_ligands_in_cluster(cluster_details):

    all_ligands = []
    total_pairings = 0
    
    protein_counts_for_ligand = {}

    ligands_for_protein = {}
    
    for idx, cluster_info in enumerate(cluster_details.items()):
        # extend list
        cluster_proteins = cluster_info[1]

        for current_protein in cluster_proteins:
            
            ligands_for_curr_protein = protein2ligand[current_protein]
            ligands_for_protein[current_protein] = ligands_for_curr_protein
            
            for ligand in ligands_for_curr_protein: #loop through each ligand
                if ligand not in protein_counts_for_ligand:
                    protein_counts_for_ligand[ligand] = 1
                else:
                    protein_counts_for_ligand[ligand] += 1
            
            all_ligands.extend(ligands_for_curr_protein)
            total_pairings += len(ligands_for_curr_protein)
    
    
    num_unique_ligands = len(set(all_ligands))

    return set(all_ligands), protein_counts_for_ligand, num_unique_ligands, ligands_for_protein


def count_overlap_(mmseqs_cluster_random_split_ligands, random_ligand_split):

    random_ligand_split = set([ln for ln, c in random_ligand_split])
    #print(random_ligand_split)
    
    intersection = random_ligand_split.intersection(mmseqs_cluster_random_split_ligands)

    return intersection


val_all_ligands, val_protein_counts_for_ligand, _, val_ligands_for_protein  = get_ligands_in_cluster(val_cluster_details)
test_all_ligands, test_protein_counts_for_ligand, _, test_ligands_for_protein = get_ligands_in_cluster(test_cluster_details)


test_ligands_for_protein['A0A830GVA3']


len(test_all_ligands)


val_overlap = count_overlap_(val_all_ligands, val_selected_ligands)
val_overlap


test_overlap = count_overlap_(test_all_ligands, test_selected_ligands)
test_overlap


incorrect_overlap = val_overlap.intersection(test_overlap)


final_val_ligands = val_overlap - incorrect_overlap
final_test_ligands = test_overlap - incorrect_overlap


def pretty_filter(ligands_for_protein, ref_ligands):

    new_ligands_for_protein = {}
    for prot_name, ligands in ligands_for_protein.items():

        if(len(set(ligands).intersection(ref_ligands))>0):
            new_ligands_for_protein[prot_name]  = ligands
    return new_ligands_for_protein



print("Final Validation Ligands:", final_val_ligands)
print("Final Test Ligands:", final_test_ligands)


# Save the final ligand sets to files
with open("ligands_val.txt", "w") as f_val:
    for ligand_id in sorted(list(final_val_ligands)): # Sort for consistent order
        f_val.write(ligand_id[:-4] + "\n")
print(f"Saved {len(final_val_ligands)} ligands to ligands_val.txt")

with open("ligands_test.txt", "w") as f_test:
    for ligand_id in sorted(list(final_test_ligands)): # Sort for consistent order
        f_test.write(ligand_id[:-4] + "\n")
        
print(f"Saved {len(final_test_ligands)} ligands to ligands_test.txt")

# pretty_filter(test_ligands_for_protein, test_overlap - incorrect_overlap)


# import pandas as pd
# import plotly.express as px

# # Create a DataFrame from sorted_ligands
# df = pd.DataFrame(sorted_ligands, columns=['ligand', 'count'])

# # Create a bar plot (see counts per ligand)
# fig_bar = px.bar(df, x='ligand', y='count', log_y=True,
#                  title=f"Associated Protein Counts of the {num_unique_ligands} ligands",
#                  labels={'ligand': 'Ligand', 'count': 'Protein Count'},
#                  color_discrete_sequence=["blue"]
#                 )
# fig_bar.update_layout(xaxis_tickangle=-45)
# fig_bar.show()
# fig_bar.write_html("/home/ruh/www/static/scientificLLM/filtered_ligand_histogram.html")