# Calculates the RMSD between all models in a directory and a set of reference structures.

#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_directory> <reference_structures_directory> <output_directory_map_pdbs> <name_clustering_file> <number_of_cores> <output_file_name>"
    exit 1
fi

# Assign arguments to variables
input_directory="$1"
reference_structures_directory="$2"
output_directory_map_pdbs="$3"
number_of_cores="$4"
output_file_name="$5"

# Run the first Python command (map_pdbs.py)
python3 calc_all_models/map_all_pdbs.py "$input_directory" "$reference_structures_directory" "$output_directory_map_pdbs" "$number_of_cores"

# Check if the first command succeeded
if [ $? -eq 0 ]; then
    echo "map_pdbs.py completed successfully."

    # Run the second Python command (calculate_rmsd.py) without output_directory_rmsd
    python3 calc_all_models/calculate_rmsd_all.py "$output_directory_map_pdbs" "$output_file_name" "$number_of_cores"

    # Check if the second command succeeded
    if [ $? -eq 0 ]; then
        echo "calculate_rmsd.py completed successfully."
    else
        echo "Error running calculate_rmsd.py"
    fi
else
    echo "Error running map_pdbs.py. Aborting."
fi
