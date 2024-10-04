#!/bin/bash

for file_A in *_A.pdb; do
    base_name="${file_A%_A.pdb}"  # Extract the base name by removing the _A.pdb suffix
    file_D="${base_name}_D.pdb"   # Assume the corresponding B file exists
    
    if [[ -f "$file_D" ]]; then   # Check if file_B exists
        output="${base_name}_merged.pdb"  # Set the output file name
        pdb_merge "$file_A" "$file_D" > "$output"  # Run the merge command
	output_n="${base_name}_merged_n.pdb"
        grep 'ATOM ' "$output" > "$output_n"
        echo "Merged $file_A and $file_D into $output"
    else
        echo "File $file_D not found, skipping $file_A"
    fi
done
