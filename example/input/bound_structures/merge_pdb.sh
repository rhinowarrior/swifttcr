#!/bin/bash

# Define paths
INPUT_DIR="/home/nils/swifttcr/example/input/bound_structures"
OUTPUT_DIR="/home/nils/swifttcr/example/input/bound_structures"


# Iterate over each PDB file in the input directory
for FILE in "$INPUT_DIR"/*_ANARCI.pdb; do
    if [ -f "$FILE" ]; then
        # Extract the base name without the _ANARCI.pdb suffix
        BASENAME=$(basename "$FILE" _ANARCI.pdb)
        
        # Define the output file names
        OUTPUT_E="$OUTPUT_DIR/${BASENAME}_Eshifted.pdb"  # For chain E (shifted)
        OUTPUT_D="$OUTPUT_DIR/${BASENAME}_D.pdb"          # For chain D (normal)
        MERGED_NAME="$OUTPUT_DIR/${BASENAME}_merged.pdb"   # Final merged output

        # Process chain E (shifted to chain D)
        echo "Processing chain E from $FILE"
        command_shift="pdb_tidy \"$FILE\" | pdb_selchain -E | pdb_shiftres -2000 | pdb_chain -D > \"$OUTPUT_E\""
        echo "Running command for chain E: $command_shift"
        eval $command_shift
        
        # Check if OUTPUT_E is empty
        if [ ! -s "$OUTPUT_E" ]; then
            echo "Warning: OUTPUT_E ($OUTPUT_E) is empty!"
        fi

        # Process chain D (normal residue numbering)
        echo "Processing chain D from $FILE"
        command_lig="pdb_tidy \"$FILE\" | pdb_selchain -D > \"$OUTPUT_D\""
        echo "Running command for chain D: $command_lig"
        eval $command_lig
        
        # Check if OUTPUT_D is empty
        if [ ! -s "$OUTPUT_D" ]; then
            echo "Warning: OUTPUT_D ($OUTPUT_D) is empty!"
        fi

        # Merge chain D and chain E
        echo "Merging chain D ($OUTPUT_D) and chain E ($OUTPUT_E) into $MERGED_NAME"
        merge_command="pdb_merge \"$OUTPUT_D\" \"$OUTPUT_E\" | pdb_sort > \"$MERGED_NAME\""
        eval $merge_command

        # Check if merged file is created successfully
        if [ -f "$MERGED_NAME" ]; then
            echo "Merged files into $MERGED_NAME"
        else
            echo "Error: Merging failed for $OUTPUT_D and $OUTPUT_E"
        fi

        # Clean up temporary files
        rm "$OUTPUT_D" "$OUTPUT_E"
        echo "Removed temporary files: $OUTPUT_D and $OUTPUT_E"
    else
        echo "No files matching *_ANARCI.pdb found in $INPUT_DIR"
    fi
done