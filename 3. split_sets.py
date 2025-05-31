import csv
import random
import plotly.io as pio
import plotly.express as px
pio.renderers.default = "colab"
from collections import defaultdict

# Set the random seed
random.seed(42)

# Create a dictionary to hold sets of sequence IDs for each cluster
clusters = defaultdict(set)

data_path = "/mnt/gemini/data/ramith/CMU-project/data/raw"

with open(f"{data_path}/data/clusterRes_cluster.tsv", "r") as tsvfile:
    reader = csv.reader(tsvfile, delimiter="\t")
    for row in reader:
        if row:
            cluster_id = row[0]
            sequence_id = row[1]
            clusters[cluster_id].add(sequence_id)

# Create a list of (cluster_id, count) tuples
cluster_counts = [(cluster_id, len(seq_set)) for cluster_id, seq_set in clusters.items()]

# Sort clusters by count descending (top clusters)
top_clusters = sorted(cluster_counts, key=lambda x: x[1], reverse=True)[:10]

# Sort clusters by count ascending (least clusters)
least_clusters = sorted(cluster_counts, key=lambda x: x[1])[:10]

len(cluster_counts)

print("Top 10 clusters:")
for cluster_id, count in top_clusters:
    print(f"Cluster {cluster_id} has {count} sequences.")

print("\nLeast 10 clusters:")
for cluster_id, count in least_clusters:
    print(f"Cluster {cluster_id} has {count} sequences.")


singleton_count = sum(1 for _, count in cluster_counts if count == 1)
print("Number of clusters with one sequence:", singleton_count)



# Extract just the counts from the (cluster_id, count) tuples
counts = [count for _, count in cluster_counts]

# Create a histogram with, say, 50 bins (adjust nbins as needed)
fig = px.histogram(
    x=counts, 
    nbins=50000,
    log_y=True,
    labels={'x': 'Number of Sequences per Cluster', 'y': 'Count of Clusters'},
    title='Histogram of Sequences per Cluster'
)

# fig.show()
fig.write_html("cluster_histogram.html")

SAMPLING_RANGES = [
    ((1, 501), 40),          # from clusters of size [1..100), pick 15
    ((501, 1001), 20),        # from clusters of size [101..200), pick 10
    ((1001, 2501), 10),         # from clusters of size [201..300), pick 5
    ((2501, float("inf")), 2),         # from clusters of size [201..300), pick 5
]

def in_range(size, low, high):
    """
    Returns True if 'size' is in [low, high), or [low, ∞) if high == float('inf').
    """
    if high == float('inf'):
        return size >= low
    else:
        return (size >= low) and (size < high)


# sample clusters for each range
selected_clusters = []  # will hold (cluster_id, cluster_size) for all chosen

for (low, high), desired_count in SAMPLING_RANGES:
    # Filter clusters that fall in the [low, high) range
    matching = [(cid, sz) for (cid, sz) in cluster_counts if in_range(sz, low, high)]
    
    if len(matching) < desired_count:
        print(f"Warning: asked for {desired_count} clusters in range [{low}, {high}), "
              f"but only {len(matching)} available. Sampling all of them.")
        desired_count = len(matching)
    
    chosen = random.sample(matching, desired_count)
    selected_clusters.extend(chosen)

val_clusters = []
test_clusters = []

start_idx = 0
for (low, high), desired_count in SAMPLING_RANGES:
    
    end_idx = start_idx + desired_count
    bin_slice = selected_clusters[start_idx:end_idx]
    start_idx = end_idx
    
    half = desired_count // 2
    val_clusters.extend(bin_slice[:half])
    test_clusters.extend(bin_slice[half:])

print(f"val_clusters has {len(val_clusters)} clusters, test_clusters has {len(test_clusters)} clusters.")

# ## Save the clusters

with open("val_clusters.txt", "w") as f_val:
    for (cid, sz) in val_clusters:
        f_val.write(cid + "\n")

with open("test_clusters.txt", "w") as f_test:
    for (cid, sz) in test_clusters:
        f_test.write(cid + "\n")


selected_cids = {cid for cid, _ in val_clusters + test_clusters}
all_cids = set(clusters.keys())
train_cids = all_cids - selected_cids

with open("train_clusters.txt", "w") as f_train:
    for cid in train_cids:
        f_train.write(cid + "\n")

print(f"train_clusters has {len(train_cids)} clusters.")