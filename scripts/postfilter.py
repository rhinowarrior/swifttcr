"""Script to filter piper results based on distance between atoms.

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
    start = time.time_ns()
    ft_file = output_dir + ft_file
    rot_file = rot_file
    res_file = res_file
    receptor = receptor
    ligand = ligand
    outfilename = output_dir + outfilename

    ft_dict = parse_ft_file(ft_file)
    rot_dict = parse_rot_file(rot_file)
    restraints = parse_res_file(res_file)
    index_to_keep = post_filter(ft_dict, rot_dict, restraints, receptor, ligand)
    #print("Total lines to keep: ", len(index_to_keep))
    filter_file_by_indices(ft_file, outfilename, index_to_keep)
    stop = time.time_ns()
    print("Time to run: ", stop-start)


def filter_file_by_indices(input_file, output_file, indices):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line_num, line in enumerate(infile, start=1):
                if line_num in indices:
                    outfile.write(line)
    except FileNotFoundError:
        print(f"File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


def parse_ft_file(ft_file):
    """ Parse ft file to get rotation index and translation (x y z)

    Returns: ft dict with key: index, value: translation
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

    Returns: rot_dict with key: index, value: rot_matrix
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

    Returns json object
    """
    f = open(res_file, 'r')
    json_line = f.read()
    f.close()
    restraints = json.loads(json_line)

    return restraints


def get_pdb_coords(model, chains, residues):
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

    """
    filtered_ft_file = []
    for key, (i, translation) in ft_dict.items():
        rot_mat = rot_dict[key]
        if check_restraints(restraints, rot_mat, translation, receptor, ligand):
            filtered_ft_file.append(i)#index of lines to keep.
    return filtered_ft_file


def check_restraints(restraints, rot_mat, translation, receptor, ligand):
    
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
            #print("Receptor: ", rec_chain, rec_resid)
            #print("Ligand: ", lig_chain, lig_resid)
            coords_receptor = get_pdb_coords(receptor, [rec_chain], [rec_resid])
            coords_ligand = get_pdb_coords(ligand, [lig_chain], [lig_resid])
            res_bool = validate(coords_receptor, coords_ligand, rot_mat, translation, dmin, dmax)
            if res_bool:
                valid +=1
        if valid >= required_group:
              total += 1
    #print("Satisfied: ", total >= total_required)
    if total >= total_required:
        return True


def validate(xyz1, xyz2, r1, translation, dmin, dmax):
    """Validate if restraints are satisfied

    """
    #take mean of coordinates
    if xyz1 == []:
        return True
    if xyz2 == []:
        return True
    mean_xyz1 = np.mean(np.array(xyz1), axis = 0)
    mean_xyz2 = np.mean(np.array(xyz2), axis = 0)
    rotation_matrix = np.array(r1).reshape(3, 3)
    #rotation_matrix = str_to_rotmat(r1)
    r = R.from_matrix(rotation_matrix)
    transformed_coords = r.apply(mean_xyz2) + np.array(translation)

    distance = coord_distance(mean_xyz1, transformed_coords)
    if distance >= dmin and distance <= dmax:
        return True
    else:
        return False


def coord_distance(coord1, coord2):
    (x1, y1, z1) = coord1
    (x2, y2, z2) = coord2
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

