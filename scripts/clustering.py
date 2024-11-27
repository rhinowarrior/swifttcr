"""Script to cluster based on Cluspro

Greedy clustering
Author: Yannick Aarts
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5540229/pdf/nihms883869.pdf
Cluspro paper:
We calculate IRMSD values for each pair 
among the 1000 structures, and find the structure that has the highest number of neighbors
within 9Å IRMSD radius. The selected structure will be defined as the center of the first
cluster, and the structures within the 9Å IRMSD neighborhood of the center will constitute
the first cluster. The members of this cluster are then removed, and we select the structure
with the highest number of neighbors within the 9Å IRMSD radius among the remaining
structures as the next cluster center. These neighbors will form the next cluster. Up to 30
clusters are generated in this manner

Example usage: python3 python_programs/clustering.py pipeline/ 7
Iterates over directories with _60_ms or KB.
"""
import csv
from typing import List, Tuple
from pathlib import Path

def read_irmsd_values(file_path: str) -> List[Tuple[str, str, float]]:
    """Read the irmsd values from file and returns a list of tuples.

    Args:
        file_path: path to the pairwise irmsd values csv file.
        
    Returns:
        irmsd_values: list of tuples with (m1, m2, irmsd).
    """
    irmsd_values = []
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            irmsd_values.append((row[0], row[1], float(row[2])))
    return irmsd_values


def create_dict(irmsd_values, threshold=9):
    """Create dictionary with key: model_name, and value [(model2, irmsd)].
    Threshold sets the maximum irmsd values included in the dictionary. (Default 9A as done in Cluspro)
    
    Args:
        irmsd_values: List of tuples (m1, m2, irmsd).
        threshold: Maximum irmsd value to include in the dictionary. (Default 9)

    Returns:
        irmsd_dict: Dictionary with key: model_name, and value [(model2, irmsd)].
    """
    irmsd_dict = {}
    for (m1, m2, value) in irmsd_values:
        if value <= threshold:
            if m1 in irmsd_dict.keys():
                irmsd_dict[m1].append((m2, value))
            else:
                irmsd_dict[m1] = [(m2, value)]
            if m2 in irmsd_dict.keys():
                irmsd_dict[m2].append((m1, value))
            else:
                irmsd_dict[m2] = [(m1, value)]
    return irmsd_dict


def cluster(irmsd_dict, nr_of_clusters=100):
    """Greedy clustering of largest clusters based on irmsd.
    Keeps track of visited models with set(). Iterates over all keys and selects the largest cluster. 
    Then removes (hides with set()) those models from the dataset.
    
    Args:
        irmsd_dict: Dictionary with key: model_name, and value [(model2, irmsd)].
        nr_of_clusters: Number of clusters to output. (Default 100)
        
    Returns:
        clusters: List of tuples with (model_name, nr_of_neighbors, [members])
    """
    
    visited = set()
    clusters = []
    while len(clusters) < nr_of_clusters:
        count_dict = {}
        for key1 in irmsd_dict.keys():
            neighbors = 0
            if key1 not in visited:
                for (model, _) in irmsd_dict[key1]:
                    if model not in visited:
                        neighbors += 1
                count_dict[key1] = neighbors
        # Get the largest cluster.
        if count_dict:
            (model, count) = max(count_dict.items(), key=lambda x: x[1])
        else:
            return clusters
        
        visited.add(model)  # Hide center model from dataset.
        members = []
        for (member, _) in irmsd_dict[model]:  # Hide all members from dataset.
            visited.add(member)
            members.append(member.split("_")[1].strip(".pdb"))
        clusters.append((model, count, members))
    return clusters


def clustering_main(input_file, treshold=9, output_file=None,nr_of_clusters=100):
    """Main function for clustering based on irmsd values.
    
    Args:
        input_file: Path to the pairwise irmsd values csv file.
        treshold: Maximum irmsd value to include in the dictionary. (Default 9)
        nr_of_clusters: Number of clusters to output. (Default 100)
        output_file: Optional path to save the clustering output. If None, output is printed only.
    
    Returns:
        clusters: List of tuples with (model_name, nr_of_neighbors, [members])
    """
    irmsd_values = read_irmsd_values(input_file)
    irmsd_dict = create_dict(irmsd_values, treshold)
    clusters = cluster(irmsd_dict, nr_of_clusters)

    output_lines = []
    for model, neighbors, _ in clusters:
        line = f"Cluster center: {model} with {neighbors} neighbors."
        output_lines.append(line)
    output_lines.append(f"Number of clusters found: {len(clusters)}")
    
    # Handle output writing
    if output_file:
        with open(output_file, 'w') as f:
            f.write('\n'.join(output_lines))
        print(f"Clustering output written to {output_file}")
    else:
        print('\n'.join(output_lines))
    
    return clusters