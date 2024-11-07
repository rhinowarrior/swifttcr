"""
Name: postfilter.py
Function: This script is used to filter piper results based on distance between atoms.
Date: 08-06-2024
Author: Yannick Aarts

Inputs: Ft file with piper results, prm file with rotations, distance restraint file.
"""

import os
import time
import math
import json
import numpy as np
from scipy.spatial.transform import Rotation as R
from multiprocessing import Pool
from pdb2sql import pdb2sql

# Cache for storing PDB coordinates at the model level
LOCAL_COORDINATE_DICTIONARY = {}

def post_filter_main(output_dir, ft_file, rot_file, res_file, receptor, ligand, outfilename, num_cores):
    """Main function to filter piper results based on distance restraints."""
    start = time.time_ns()
    
    # Build full paths for input and output files
    ft_file = os.path.join(output_dir, ft_file)
    rot_file = os.path.join(output_dir, rot_file)
    res_file = os.path.join(output_dir, res_file)
    outfilename = os.path.join(output_dir, outfilename)

    # Parse the ft, rot, and res files
    ft_dict = parse_ft_file(ft_file)
    rot_dict = parse_rot_file(rot_file)
    restraints = parse_res_file(res_file)

    # Filter the ft file based on the restraints using multiprocessing
    indices_to_keep = post_filter(ft_dict, rot_dict, restraints, receptor, ligand, num_cores)

    # Write the filtered ft file to the output file
    filter_file_by_indices(ft_file, outfilename, indices_to_keep)
    
    stop = time.time_ns()
    print("Time to run: ", stop - start)


def post_filter(ft_dict, rot_dict, restraints, receptor, ligand, num_cores):
    """Apply restraints to receptor and rotated + translated ligand using multiprocessing."""
    # Combine tasks into batches of reasonable size
    batch_size = max(10, len(ft_dict) // (num_cores * 2))  # Increase batch size for larger tasks
    tasks = create_batch_tasks(ft_dict, rot_dict, restraints, receptor, ligand, batch_size)

    # Dynamically set the number of processes based on available cores
    num_processes = num_cores  # Use user-defined number of cores directly
    
    # Calculate chunksize to avoid large memory loads
    chunksize = max(1, len(tasks) // num_processes)

    with Pool(processes=num_processes) as pool:
        # Here, we use `chunksize` to control memory usage across processes
        results = pool.map(check_batch_restraints, tasks, chunksize=chunksize)
    
    # Flatten results and keep indices that satisfy restraints
    indices_to_keep = [i for batch_result in results for i in batch_result]
    return indices_to_keep


def create_batch_tasks(ft_dict, rot_dict, restraints, receptor, ligand, batch_size):
    """Create batch tasks for multiprocessing to reduce overhead."""
    # Group tasks into batches
    tasks = []
    current_batch = []
    
    for idx, (i, translation) in ft_dict.items():
        current_batch.append((idx, translation, rot_dict[idx], restraints, receptor, ligand))
        
        if len(current_batch) >= batch_size:
            tasks.append(current_batch)
            current_batch = []
    
    if current_batch:
        tasks.append(current_batch)  # Add the last batch if it has remaining tasks
    
    return tasks


def check_batch_restraints(task_batch):
    """Process a batch of restraint checks."""
    indices_to_keep = []
    for key, translation, rot_mat, restraints, receptor, ligand in task_batch:
        if check_restraints(restraints, rot_mat, translation, receptor, ligand):
            indices_to_keep.append(key)
    return indices_to_keep


def check_restraints(restraints, rot_mat, translation, receptor, ligand):
    """Check if restraints are satisfied for a single transformation."""
    groups = restraints['groups']
    total_required = restraints['required']
    total = 0
    for group in groups:
        required_group = group['required']
        valid = 0
        group_res = group['restraints']
        for group_r in group_res:
            rec_chain = group_r['rec_chain']
            rec_resid = group_r['rec_resid']
            lig_chain = group_r['lig_chain']
            lig_resid = group_r['lig_resid']
            dmax = group_r['dmax']
            dmin = group_r['dmin']
            coords_receptor = get_cached_coords(receptor, [rec_chain], [rec_resid])
            coords_ligand = get_cached_coords(ligand, [lig_chain], [lig_resid])
            if validate(coords_receptor, coords_ligand, rot_mat, translation, dmin, dmax):
                valid += 1
        if valid >= required_group:
            total += 1
    return total >= total_required


def get_cached_coords(model, chains, residues):
    """Get coordinates of CA atoms in pdb file using cache."""
    cache_key = (model, tuple(chains), tuple(residues))
    
    if cache_key not in LOCAL_COORDINATE_DICTIONARY:
        coords = get_pdb_coords(model, chains, residues)
        LOCAL_COORDINATE_DICTIONARY[cache_key] = coords
    
    return LOCAL_COORDINATE_DICTIONARY[cache_key]


def get_pdb_coords(model, chains, residues):
    """Get coordinates of CA atoms in pdb file using pdb2sql."""
    if model not in LOCAL_COORDINATE_DICTIONARY:
        pdb = pdb2sql(model)
        LOCAL_COORDINATE_DICTIONARY[model] = pdb
    pdb = LOCAL_COORDINATE_DICTIONARY[model]
    
    # Query for CA atoms in the specific chains and residues
    xyz = pdb.get('x,y,z', chainID=chains, resSeq=residues, name=['CA'])
    return xyz


def validate(xyz1, xyz2, r1, translation, dmin, dmax):
    """Validate if restraints are satisfied."""
    if not xyz1 or not xyz2:
        return True
    mean_xyz1 = np.mean(np.array(xyz1), axis=0)
    mean_xyz2 = np.mean(np.array(xyz2), axis=0)
    rotation_matrix = np.array(r1).reshape(3, 3)
    r = R.from_matrix(rotation_matrix)
    transformed_coords = r.apply(mean_xyz2) + np.array(translation)
    distance = coord_distance(mean_xyz1, transformed_coords)
    return dmin <= distance <= dmax


def coord_distance(coord1, coord2):
    """Calculate distance between two coordinates."""
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)))


def parse_ft_file(ft_file):
    """Parse ft file to get rotation index and translation (x y z)"""
    ft_dict = {}
    with open(ft_file, 'r') as f:
        for i, line in enumerate(f):
            values = line.strip().split('\t')
            if len(values) >= 4:
                index = int(values[0])
                x_translation = float(values[1])
                y_translation = float(values[2])
                z_translation = float(values[3])
                ft_dict[index] = (i, (x_translation, y_translation, z_translation))
    return ft_dict


def parse_rot_file(rot_file):
    """Parse rot file to get rotation matrix."""
    rot_dict = {}
    with open(rot_file, 'r') as f:
        for line in f:
            values = line.strip().split()
            if len(values) >= 10:
                index = int(values[0])
                matrix_values = list(map(float, values[1:]))
                rot_dict[index] = matrix_values
    return rot_dict


def parse_res_file(res_file):
    """Parse restraints file."""
    with open(res_file, 'r') as f:
        json_line = f.read()
    restraints = json.loads(json_line)
    return restraints


def filter_file_by_indices(input_file, output_file, indices):
    """Filter lines in input file based on indices and write to output file."""
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line_num, line in enumerate(infile, start=1):
                if line_num in indices:
                    outfile.write(line)
    except FileNotFoundError:
        print(f"File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
