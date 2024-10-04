import os
import subprocess

# Define the directory containing your PDB files
pdb_directory = '/home/nils/swifttcr/example/input/bound_structures'  # Change this to your directory
pdb_files = os.listdir(pdb_directory)

# Create a dictionary to hold the pairs
file_pairs = {}

# Find files and group them by the prefix before the first underscore
for file in pdb_files:
    if file.endswith('.pdb'):
        prefix = file.split('_')[0]  # Get the prefix before the first underscore
        if prefix not in file_pairs:
            file_pairs[prefix] = []
        file_pairs[prefix].append(file)

print(prefix)

# Merge files based on pairs found
for prefix, files in file_pairs.items():
    if prefix + '_b_A.pdb' in files:
        file1 = os.path.join(pdb_directory, prefix + '_b_A.pdb')
    if prefix + '_b_DE_merged.pdb' in files:
        file2 = os.path.join(pdb_directory, prefix + '_b_DE_merged.pdb')
    output_file = os.path.join(pdb_directory, f'{prefix}_merged.pdb')

        # Call pdb_merge
    print(f'Merging {file1} and {file2} into {output_file}')
    command = f'pdb_merge {file1} {file2} > {output_file}'
    
    os.system(command)
