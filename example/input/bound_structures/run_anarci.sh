#!/bin/bash

# Directory containing PDB files
pdb_dir="/home/nils/swifttcr/example/input/bound_structures"  # Change this to the directory you want to search

# Python script that runs ANARCI (modify this if your ANARCI script is located elsewhere)
anarci_script="/home/nils/swifttcr/tools/ANARCI_master/Example_scripts_and_sequences/ImmunoPDB.py"  # Path to the script where immunopdb_main is defined

# Loop over all files matching *_b_DE.pdb
for pdb_file in "$pdb_dir"/*_b_DE.pdb; do
    # Get the base name of the file (without the extension)
    base_name=$(basename "$pdb_file" .pdb)

    # Define the output file name
    output_file="$pdb_dir/${base_name}_ANARCI.pdb"

    # Run the ANARCI script with the input PDB file and output file
    python3 "$anarci_script" -i "$pdb_file" -o "$output_file" -r tr
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $pdb_file"
    else
        echo "Failed to process $pdb_file"
    fi
done
