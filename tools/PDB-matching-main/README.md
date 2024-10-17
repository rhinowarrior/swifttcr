# PDB-matching

## Overview
**PDB-matching** is a tool developed to map a structure against a reference structure. After this the RMSD can be calculated over the mapped structures.

## Features
- Make a mapping of the reference an target structure that are put into a new directory.
- Calculating RMSD over the mapped files

## Getting Started

To get started with PDB-matching, follow these steps

1. **Clone or Download** this repository
2. **Download** The neccesary packages

## Running PDB-matching With Clustering File
Use the following command to map the decoy structures with the reference structures using the clustering files.

```bash
python3 calc_with_cluster_files/map_pdbs.py input_directory reference_structures_directory output_directory name_clustering_file number_of_cores 
```

After this command has run you can use the following command to calculate the rmsd on the output directory of the mapped code.

```bash
python3 calc_with_cluster_files/calculate_rmsd.py output_directory_from_map_pdbs input_directory name_clustering_file output_file_name number_of_cores
```

Using the bash script `calc_rmsd_cluster_files.sh` it is possible to run both these scripts in succesion. To do that this format should be used.

```bash
bash calc_rmsd_cluster_files.sh input_directory reference_structures_directory output_directory_map_pdbs name_clustering_file number_of_cores output_file_name"
```


## Running PDB-matching With Clustering File
Use the following command to map the decoy strucutres with the reference structures for al the files in the input directory.

```bash
python3 calc_all_models/map_all_pdbs.py input_directory reference_structure_directory output_directory number_of_cores
```

After this command has run you can use the following command to calculate the rmsd of all the files in the outputdir.

```bash
python3 calc_all_models/calculate_rmsd_all.py output_directory_from_map_pdbs output_file_name number_of_cores
```

Using the bash script `calc_rmsd_all.sh` it is possible to run both these scripts in succesion. To do that this format should be used.

```bash
bash calc_rmsd_all.sh input_directory reference_structures_directory output_directory_map_pdbs number_of_cores output_file_name
```

## Dependecies
* Python 3.9.12
* [pdb2sql: 0.5.3](https://github.com/DeepRank/pdb2sql)
* [Biopython: 1.79](https://biopython.org/)
