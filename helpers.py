import csv
import pandas as pd
import plotly.express as px
from collections import defaultdict

def load_cluster_split(file_name):
    cluster = set()
    with open(file_name, "r") as f:
        for line in f:
            cid = line.strip()
            if cid:
                cluster.add(cid)
    return cluster

def get_uniprot_ids_of_cluster(ref_clusters: set, cluster_mapping_file: str, delimiter: str = "\t"):
    """
    Extracts and returns UniProt IDs for the given set of cluster IDs from
    a tabular file (e.g. clusterRes_cluster.tsv).

    Assumes:
    - The file has at least two columns: [cluster_id, uniprot_id, ...]
    - No header row.
    - Each row corresponds to (cluster_id, uniprot_id).

    :param clusters: Set of cluster IDs you want to filter on.
    :param cluster_mapping_file: Path to the TSV/CSV file containing the cluster-to-UniProt mapping.
    :param delimiter: Delimiter used in the cluster mapping file. Defaults to tab.
    :return: List of UniProt IDs (duplicates may appear if multiple rows match).
    """
    cluster_details = defaultdict(set)
    
    with open(cluster_mapping_file, "r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        for row in reader:
            if row:
                cluster_id = row[0]
                sequence_id = row[1]

                if(cluster_id in ref_clusters):
                    cluster_details[cluster_id].add(sequence_id)

    # Create a list of (cluster_id, count) tuples
    cluster_counts = [(cluster_id, len(seq_set)) for cluster_id, seq_set in cluster_details.items()]

    return cluster_details, cluster_counts

def plot_cluster_distribution(cluster_counts, text = ""):
    import plotly.express as px

    # Extract just the counts from the (cluster_id, count) tuples
    counts = [count for _, count in cluster_counts]
    
    # Create a histogram with, say, 50 bins (adjust nbins as needed)
    fig = px.histogram(
        x=counts, 
        nbins=100000,
        log_y=True,
        labels={'x': f"Number of Sequences per {text} cluster", 'y': 'Count of Clusters'},
        title=f"Histogram of Sequences per {text} cluster"
    )
    
    return fig

def plot_cluster_distributions(cluster_counts_list, labels_list):
    """
    Create an overlaid histogram comparing the distribution of cluster sizes
    from multiple cluster_counts sources.
    """

    # Concatenate data from each dataset into a single DataFrame
    df = pd.DataFrame(columns=["cluster_id", "count", "source"])
    for i, cluster_counts in enumerate(cluster_counts_list):
        temp_df = pd.DataFrame(cluster_counts, columns=["cluster_id", "count"])
        temp_df["source"] = labels_list[i]
        df = pd.concat([df, temp_df], ignore_index=True)

    # Plot the combined data as an overlaid histogram
    fig = px.histogram(
        df,
        x="count",
        color="source",
        barmode="overlay",  # overlay all histograms
        nbins=50,           # adjust as needed
        log_y=True,         # log scale on y-axis
        labels={
            "count": "Number of Sequences per Cluster",
            "source": "Cluster Source"
        },
        title="Histogram of Sequences per Cluster (Comparison)"
    )
    fig.update_layout(legend_title_text="Dataset")

    return fig
