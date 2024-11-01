"""
Name: postfilter.py
Function: This script is used to filter piper results based on distance between atoms.
Date: 08-06-2024
Author: Yannick Aarts

Inputs: Ft file with piper results, prm file with rotations, distance restraint file.
"""

import time
import math
import json
from pdb2sql import pdb2sql
import numpy as np
from scipy.spatial.transform import Rotation as R

LOCAL_COORDINATE_DICTIONARY = {}#local dictionary to avoid duplicate lookups in pdb files.

def post_filter_main(output_dir, ft_file, rot_file, res_file, receptor, ligand, outfilename):
    """Main function to filter piper results based on distance restraints.

    Args:
        output_dir (str): The path to the output directory
        ft_file (str): The path to the ft files from piper
        rot_file (str): The path to the prm file with rotations
        res_file (str): The path to the distance restraint file
        receptor (str): The path to the receptor pdb file
        ligand (str): The path to the ligand pdb file
        outfilename (str): The name of the output file
    """
    start = time.time_ns()
    ft_file = output_dir + ft_file
    rot_file = rot_file
    res_file = res_file
    receptor = receptor
    ligand = ligand
    outfilename = output_dir + outfilename

    # Parse the ft, rot and res files
    ft_dict = parse_ft_file(ft_file)
    rot_dict = parse_rot_file(rot_file)
    restraints = parse_res_file(res_file)
    
    # Filter the ft file based on the restraints
    index_to_keep = post_filter(ft_dict, rot_dict, restraints, receptor, ligand)

    # Write the filtered ft file to the output file
    filter_file_by_indices(ft_file, outfilename, index_to_keep)
    stop = time.time_ns()
    print("Time to run: ", stop-start)


def filter_file_by_indices(input_file, output_file, indices):
    """Filter lines in input file based on indices and write to output file.

    Args:
        input_file (str): path to input file
        output_file (str): path to output file
        indices (list): list of indices to keep
        
    Raises:
        FileNotFoundError: If the input file is not found
        Exception: If an error occurs while reading or writing
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Enumerate the lines in the input file and write the lines with the indices to the output file
            for line_num, line in enumerate(infile, start=1):
                if line_num in indices:
                    outfile.write(line)
    except FileNotFoundError:
        print(f"File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


def parse_ft_file(ft_file):
    """ Parse ft file to get rotation index and translation (x y z)

    Args:
        ft_file (str): path to ft file
        
    Returns:
        dict: with key: index, value: translation
    """
    ft_dict = {}
    with open(ft_file, 'r') as f:
        for i, line in enumerate(f):
            values = line.strip().split('\t')
            # Check if the line contains at least 4 values
            if len(values) >= 4:
                index = int(values[0])
                x_translation = float(values[1])
                y_translation = float(values[2])
                z_translation = float(values[3])
                
                # Store the translations in the dictionary with the index as the key
                ft_dict[index] = (i, (x_translation, y_translation, z_translation))

    return ft_dict


def parse_rot_file(rot_file):
    """Parse rot file to get rotation matrix.
    
    Args:
        rot_file (str): path to rot file

    Returns: 
        dict: with key: index, value: rot_matrix
    """
    rot_dict = {}
    with open(rot_file, 'r') as f:
        for line in f:
            # Split each line into values using space as the delimiter
            values = line.strip().split()
            
            # Check if the line contains at least 10 values (1 index + 9 matrix values)
            if len(values) >= 10:
                index = int(values[0])
                matrix_values = list(map(float, values[1:]))
            
                # Store the rotation matrix in the dictionary with the index as the key
                rot_dict[index] = matrix_values
    return rot_dict


def parse_res_file(res_file):
    """Parse restraints file 

    Args:
        res_file (str): path to restraints file

    Returns:
        json object: restraints
    """
    f = open(res_file, 'r')
    json_line = f.read()
    f.close()
    restraints = json.loads(json_line)

    return restraints


def get_pdb_coords(model, chains, residues):
    """Get coordinates of CA atoms in pdb file.

    Args:
        model (pdb2sql): pdb2sql object
        chains (list): list of chain IDs
        residues (list): list of residue IDs

    Returns:
       list : list of coordinates
    """
    #optionally implement speedup with dictionary instead of repeated receptor and ligand lookup
    if model in LOCAL_COORDINATE_DICTIONARY:
        pdb = LOCAL_COORDINATE_DICTIONARY[model]
    else:
        pdb = pdb2sql(model)
        LOCAL_COORDINATE_DICTIONARY[model] = pdb
    xyz = pdb.get('x,y,z', chainID = chains, resSeq = residues, name = ['CA'])
    return xyz


def post_filter(ft_dict, rot_dict, restraints, receptor, ligand):
    """Apply restraints to receptor and rotated + translated ligand.

    Keep rotations that satisfy restraints.
    
    Args:
        ft_dict (dict): ft file with piper results
        rot_dict (dict): prm file with rotations
        restraints (json): distance restraints
        receptor (str): path to receptor pdb file
        ligand (str): path to ligand pdb file
    
    Returns:
        list: list of indices to keep
    """
    filtered_ft_file = []
    for key, (i, translation) in ft_dict.items():
        rot_mat = rot_dict[key]
        if check_restraints(restraints, rot_mat, translation, receptor, ligand):
            filtered_ft_file.append(i)#index of lines to keep.
    return filtered_ft_file


def check_restraints(restraints, rot_mat, translation, receptor, ligand):
    """Check if restraints are satisfied.

    Args:
        restraints (json): distance restraints
        rot_mat (list): rotation matrix
        translation (tuple): translation
        receptor (str): path to receptor pdb file
        ligand (str): path to ligand pdb file

    Returns:
        bool: True if restraints are satisfied, False otherwise
    """
    groups = restraints['groups']
    total_required = restraints['required']
    total = 0
    for group in groups:
        required_group = group['required']
        valid = 0
        group_res = group['restraints']
        for group_r in group_res:
            #get the restraints
            rec_chain = group_r['rec_chain']
            rec_resid = group_r['rec_resid']
            lig_chain = group_r['lig_chain']
            lig_resid = group_r['lig_resid']
            dmax = group_r['dmax']
            dmin = group_r['dmin']
            #get the coordinates
            coords_receptor = get_pdb_coords(receptor, [rec_chain], [rec_resid])
            coords_ligand = get_pdb_coords(ligand, [lig_chain], [lig_resid])
            #check if restraints are satisfied
            res_bool = validate(coords_receptor, coords_ligand, rot_mat, translation, dmin, dmax)
            if res_bool:
                valid +=1
        if valid >= required_group:
              total += 1
    if total >= total_required:
        return True


def validate(xyz1, xyz2, r1, translation, dmin, dmax):
    """Validate if restraints are satisfied
    
    Args:
        xyz1 (list): coordinates of receptor
        xyz2 (list): coordinates of ligand
        r1 (list): rotation matrix
        translation (tuple): translation
        dmin (float): minimum distance
        dmax (float): maximum distance
    
    returns:
        bool: True if restraints are satisfied, False otherwise
    """
    #take mean of coordinates
    if xyz1 == []:
        return True
    if xyz2 == []:
        return True
    #mean of coordinates
    mean_xyz1 = np.mean(np.array(xyz1), axis = 0)
    mean_xyz2 = np.mean(np.array(xyz2), axis = 0)
    rotation_matrix = np.array(r1).reshape(3, 3)
    #apply rotation and translation
    r = R.from_matrix(rotation_matrix)
    #transorm the coordinates
    transformed_coords = r.apply(mean_xyz2) + np.array(translation)

    distance = coord_distance(mean_xyz1, transformed_coords)
    if distance >= dmin and distance <= dmax:
        return True
    else:
        return False


def coord_distance(coord1, coord2):
    """Calculate distance between two coordinates

    Args:
        coord1 (tuple): coordinates of first atom
        coord2 (tuple): coordinates of second atom

    Returns:
        float: distance between the two coordinates
    """
    (x1, y1, z1) = coord1
    (x2, y2, z2) = coord2
    # Calculate the distance between the two coordinates
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

