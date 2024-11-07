"""
Name: apply_results.py
Function: This script is used to apply the results of the ftmap calculations to the pdb files. The script reads the results from the ftmap calculations and applies the results to the pdb files. The output is the pdb files with the results applied.
Date: 25-09-2024
Author: Yannick Aarts
"""

from argparse import FileType
import numpy as np
from itertools import islice
import linecache
from multiprocessing import Pool, Manager
from Bio.PDB import PDBParser, PDBIO, Structure

FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3)), ('E', 'f8')])

def read_rotations(file_or_handle):
    """Reads 3x3 rotation matrices from a file.

    The file may either be a text file with or without the index as
    the first column, or a numpy file. In case it is a numpy file, it
    is assumed to have the correct shape.
    
    Args:
        file_or_handle: The file or file handle to read from.
    
    Returns:
        rotations: The rotation matrices as a numpy array.
    """
    
    rotations = np.loadtxt(file_or_handle)
    if rotations.shape[-1] == 10:
        rotations = rotations[:, 1:]
    return rotations.reshape(-1, 3, 3)

def read_ftresults(filepath, limit=None):
    """Reads ftresults from a file.

    See read_ftresults_stream for details.
    
    Args:
        filepath: The path to the file to read from.
        limit: The maximum number of results to read.
    
    Returns:
        ftresults: The ftresults as a numpy array.
    """
    with open(filepath, "r") as f:
        return read_ftresults_stream(f, limit)


def read_ftresults_stream(f, limit=None):
    """Read ftresults from a stream.

    Ftresults are assumed to be in a text file with at least 5
    columns.  The first column will be the rotation index. The next
    three columns are the translation vector, and the last column is
    the total weighted energy.
    
    Args:
        f: The stream to read from.
        limit: The maximum number of results to read.
    
    Returns:
        ftresults: The ftresults as a numpy array.
    """
    f = iter(f)
    return np.loadtxt(islice(f, 0, limit), dtype=FTRESULT_DTYPE, usecols=(0, 1, 2, 3, 4))


def get_ftresult(filepath, index):
    """Get ftresult at index from file.

    index should be zero offset.
    
    Args:
        filepath: The path to the file to read from.
        index: The index of the ftresult to read.
    
    Returns:
        ftresult: The ftresult as a numpy record.
    """
    line = linecache.getline(filepath, index + 1)
    if not line:
        return None
    ss = line.strip().split()
    return np.array((int(ss[0]), [float(c) for c in ss[1:4]], float(ss[4])), dtype=FTRESULT_DTYPE)


def apply_ftresult(coords, ftresult, rotations, center=None):
    """Apply the ftresult to coords.

    `coords` and `out` cannot point to the same numpy array.
    
    Args:
        coords: The coordinates to apply the ftresult to.
        ftresult: The ftresult to apply.
        rotations: The rotation matrices.
        center: The center of the rotation.
        out: The array to write the output to.
    
    Returns:
        out: The coordinates with the ftresult applied.
    """
    if center is None:
        center = np.mean(coords, axis=0)
    coords_centered = coords - center
    rotated_coords = np.dot(coords_centered, rotations[ftresult['roti']].T)
    transformed_coords = rotated_coords + ftresult['tv'] + center
    return transformed_coords


def apply_results_worker(args):
    """Worker function to apply results in parallel.
    
    Args:
        args: A tuple with the index, ftfile, rotations, pdb_file, and out_prefix
    """
    index, ftfile, rotations, pdb_file, out_prefix, ftresults_shared = args
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure("structure", pdb_file)

    # Get the first model and its coordinates
    model = structure[0]
    coords = np.array([atom.coord for atom in model.get_atoms()])
    center = np.mean(coords, axis=0)
    
    ftresult = get_ftresult(ftfile, index)
    new_coords = apply_ftresult(coords, ftresult, rotations, center)

    # Set new coordinates
    for i, atom in enumerate(model.get_atoms()):
        atom.coord = new_coords[i]

    # Write PDB output
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"{out_prefix}.{index}.pdb")

def apply_results_main(limit, index, rotation, out_prefix, ftfile, rotations_path, pdb_file, cores):
    """ Apply the results of the ftmap calculations to the pdb files.
    
    Args:
        limit (int): The maximum number of results to read.
        index (int): The index of the ftresult to read.
        rotation (list): The rotation matrix.
        out_prefix (str): The prefix of the output file.
        ftfile (str): The path to the file to read from.
        rotations (str): The path to the file to read the rotations from.
        pdb_file (str): The path to the pdb file to read from.
    """
    
    manager = Manager()
    ftresults_shared = manager.list()  # Shared list between processes
    rotations = read_rotations(rotations_path)

    if index is not None:
        ftresults_shared.append(get_ftresult(ftfile, index))    
    elif rotation is not None:
        parser = PDBParser(QUIET=True)
        structure: Structure = parser.get_structure("structure", pdb_file)
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_prefix + ".pdb")
        return
    else:
        ftresults = read_ftresults(ftfile, limit=limit)
        for i in range(len(ftresults)):
            ftresults_shared.append(ftresults[i])
        
    tasks = [(i, ftfile, rotations, pdb_file, out_prefix, ftresults_shared) for i in range(len(ftresults_shared))]
    with Pool(cores) as pool:
        pool.map(apply_results_worker, tasks) 