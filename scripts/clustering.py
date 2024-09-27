"""Script to cluster based on Cluspro

Greedy clustering
Author: Jan Aarts
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
from sys import argv
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


def create_dict(irmsd_values, threshold = 9):
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


def cluster(irmsd_dict, nr_of_clusters = 100):
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
    while len(clusters)< nr_of_clusters:#nr of clusters
        count_dict = {}
        for key1 in irmsd_dict.keys():
            neighbors = 0
            if not key1 in visited:
                for (model, _) in irmsd_dict[key1]:
                    if not model in visited:
                        neighbors += 1
                count_dict[key1] = neighbors
        #get the largest cluster.
        if len(count_dict.items())>0:
            (model, count) = sorted(count_dict.items(),key=lambda x:x[1], reverse=True)[0]
        else:
            return clusters
        
        visited.add(model)#hide center model from dataset.
        members = []
        for (member, _) in irmsd_dict[model]:#hide all members from dataset.
            visited.add(member)
            members.append(member.split("_")[1].strip(".pdb"))
        clusters.append((model, count, members))
    return clusters


def clustering_main(input_file, treshold = 9,directory = None):
    """Main function for clustering based on irmsd values.
    
    Args:
        input_file: Path to the pairwise irmsd values csv file.
        directory: Path to the directory with the irmsd values. (Default None)
    
    Returns:
        clusters: List of tuples with (model_name, nr_of_neighbors, [members])
    """
    if directory == None:
        irmsd_values = read_irmsd_values(input_file)
    # If directory is given, iterate over all files in the directory. but is never given maybe remove this.
    elif directory != None:
        pipeline_dir = argv[1]
        p = Path(pipeline_dir)
        for f in p.iterdir():
            if "_60" in f.stem or "KB" in f.stem or "_120" in f.stem  or "_TCR_" in f.stem or "6mpp_1"and f.is_dir():
                try:
                    fname = "clustering_{}.txt".format(int(argv[2]))
                    fh = open(str(Path(f, fname)), 'w')
                    memberlistname = "member_list_{}.txt".format(int(argv[2]))
                    fm = open(str(Path(f, memberlistname)), 'w')
                    irmsd_values = read_irmsd_values(str(Path(f, 'irmsd.csv')))
                    irmsd_dict = create_dict(irmsd_values, int(argv[2]))
                    if irmsd_dict.keys():
                        clusters = cluster(irmsd_dict)
                        count_dict ={}
                        for key in irmsd_dict.keys():
                            neighbors = 0
                        for (model, _) in irmsd_dict[key]:
                                neighbors += 1
                        count_dict[key] = neighbors
                        for model, neighbors, members in clusters:
                            outline = "Cluster center: {} with {} neighbors.".format(model, neighbors)
                            model_i = model.strip(".pdb").split("_")[1]
                            member_final = "{} {}".format(model_i, members)
                            print(outline)
                            fh.write(outline + "\n")
                            fm.write(member_final + "\n")
                            
                    else:
                        clusters = []
                    finalline = "Number of clusters found: {}".format( len(clusters))
                    fh.write(finalline)
                    fh.close()
                    fm.close()
                    print(finalline)
                except FileNotFoundError:
                    print("File not found")
    else:
        exit(0)
    # If no directory is given, return the clusters.
    irmsd_dict = create_dict(irmsd_values, treshold)
    clusters = cluster(irmsd_dict)
    count_dict ={}
    for key in irmsd_dict.keys():
            neighbors = 0
            for (model, _) in irmsd_dict[key]:
                    neighbors += 1
            count_dict[key] = neighbors
    for model, neighbors, members in clusters:
        print("Cluster center: {} with {} neighbors.".format(model, neighbors))
    print("Number of clusters found: ", len(clusters))