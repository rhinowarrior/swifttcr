"""
Name: apply_results.py
Function: This script is used to apply the results of the ftmap calculations to the pdb files. The script reads the results from the ftmap calculations and applies the results to the pdb files. The output is the pdb files with the results applied.
Date: 25-09-2024
Author: Yannick Aarts
"""
from argparse import FileType
from prody import parsePDBStream, writePDB
import numpy as np
from itertools import islice
import linecache

FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3)), ('E', 'f8')])
from signal import signal, SIGPIPE, SIG_DFL

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

    return np.loadtxt(
        islice(f, 0, limit),
        dtype=FTRESULT_DTYPE,
        usecols=(0, 1, 2, 3, 4))


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
    return np.array(
        (int(ss[0]), [float(c) for c in ss[1:4]], float(ss[4])),
        dtype=FTRESULT_DTYPE)


def apply_ftresult(coords, ftresult, rotations, center=None, out=None):
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

    if out is None:
        out = np.empty_like(coords)

    # Apply the rotation
    out = np.dot(coords - center, rotations[ftresult['roti']].T)
    np.add(out, ftresult['tv'] + center, out)

    return out  


def apply_ftresults_atom_group(ag,
                               ftresults,
                               rotations,
                               center=None,
                               out=None):
    """Apply ftresult(s) to an atomgroup, returning a new atomgroup.

    The new atomgroup will have one coordinate set for each ftresult passed.
    ftresult can either be a single ftresult object, or an array of ftresult
    objects.
    
    Args:
        ag: The atomgroup to apply the ftresults to.
        ftresults: The ftresults to apply.
        rotations: The rotation matrices.
        center: The center of the rotation.
        out: The atomgroup to write the output to.
    
    Raises:
        ValueError: If ftresults is not an ndarray or void.
        
    Returns:
        out: The atomgroup with the ftresults applied.
    """
    orig_coords = ag.getCoords()  # This returns a copy so we can mutate it
    if center is None:
        center = np.mean(orig_coords, axis=0)
    np.subtract(orig_coords, center, orig_coords)

    # Ensure ftresults is a 1D array
    try:
        if len(ftresults.shape) == 0:
            ftresults = np.expand_dims(ftresults, 0)
    except:
        raise ValueError("ftresults does not seem to be an ndarray or void")

    if out is None:
        out = ag.copy()

    # Apply the rotations
    new_coords = np.dot(rotations[ftresults['roti']],
                        orig_coords.T).transpose(0, 2, 1)
    np.add(new_coords, np.expand_dims(ftresults['tv'] + center, 1), new_coords)
    out._setCoords(new_coords, overwrite=True)

    return out


signal(SIGPIPE, SIG_DFL)

def apply_results_main(limit, index, rotation, out_prefix, ftfile, rotations, pdb_file):
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

    pdb_file = FileType('r')(pdb_file)
    pdb = parsePDBStream(pdb_file)
    rotations = read_rotations(rotations)

    coords = pdb.getCoords()
    center = np.mean(coords, axis=0)
    np.subtract(coords, center, coords)
    if index is not None:
        ftresult = get_ftresult(ftfile, index)

        coords = pdb.getCoords()
        pdb.setCoords(apply_ftresult(coords, ftresult, rotations))

        writePDB(out_prefix+".pdb", pdb)
    elif rotation is not None:
        ftreults = read_ftresults(ftfile)

        writePDB(out_prefix+".pdb", pdb)
    else:
        ftresults = read_ftresults(ftfile, limit=limit)

        new_ag = apply_ftresults_atom_group(pdb, ftresults, rotations)

        # Write the new atomgroup to a PDB file
        for i in range(new_ag.numCoordsets()):
            new_ag.setACSIndex(i)
            writePDB("{}.{}.pdb".format(out_prefix, i), new_ag, csets=i)