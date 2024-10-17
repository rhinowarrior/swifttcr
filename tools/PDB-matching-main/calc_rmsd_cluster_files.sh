#This script runs the two Python scripts map_pdbs.py and calculate_rmsd.py in #sequence. to calculate the RMSD between the structures in the input directory #and the reference structures in the reference_structures_directory. The output #is saved in the output_file_name.

#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_directory> <reference_structures_directory> <output_directory_map_pdbs> <name_clustering_file> <number_of_cores> <output_file_name>"
    exit 1
fi

# Assign arguments to variables
input_directory="$1"
reference_structures_directory="$2"
output_directory_map_pdbs="$3"
name_clustering_file="$4"
number_of_cores="$5"
output_file_name="$6"

# Run the first Python command (map_pdbs.py)
python3 calc_with_cluster_files/map_pdbs.py "$input_directory" "$reference_structures_directory" "$output_directory_map_pdbs" "$name_clustering_file" "$number_of_cores"

# Check if the first command succeeded
if [ $? -eq 0 ]; then
    echo "map_pdbs.py completed successfully."

    # Run the second Python command (calculate_rmsd.py) without output_directory_rmsd
    python3 calc_with_cluster_files/calculate_rmsd.py "$output_directory_map_pdbs" "$input_directory" "$name_clustering_file" "$output_file_name" "$number_of_cores"

    # Check if the second command succeeded
    if [ $? -eq 0 ]; then
        echo "calculate_rmsd.py completed successfully."
    else
        echo "Error running calculate_rmsd.py"
    fi
else
    echo "Error running map_pdbs.py. Aborting."
fi
