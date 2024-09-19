
from argparse import ArgumentParser, FileType
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
    """
    """
    

    try:
        rotations = np.load(file_or_handle)
        return rotations
    except IOError:
        pass
    """
    rotations = np.loadtxt(file_or_handle)

    if rotations.shape[-1] == 10:
        rotations = rotations[:, 1:]

    return rotations.reshape(-1, 3, 3)


def read_ftresults(filepath, limit=None):
    """Reads ftresults from a file.

    See read_ftresults_stream for details."""
    with open(filepath, "r") as f:
        return read_ftresults_stream(f, limit)


def read_ftresults_stream(f, limit=None):
    """Read ftresults from a stream.

    Ftresults are assumed to be in a text file with at least 5
    columns.  The first column will be the rotation index. The next
    three columns are the translation vector, and the last column is
    the total weighted energy.
    """
    f = iter(f)

    return np.loadtxt(
        islice(f, 0, limit),
        dtype=FTRESULT_DTYPE,
        usecols=(0, 1, 2, 3, 4))


def get_ftresult(filepath, index):
    """Get ftresult at index from file.

    index should be zero offset.
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
    """
    if center is None:
        center = np.mean(coords, axis=0)

    if out is None:
        out = np.empty_like(coords)

    out = np.dot(coords - center, rotations[ftresult['roti']].T)
    np.add(out, ftresult['tv'] + center, out)
    #np.add(out, ftresult['tv'], out)

    return out  


def apply_ftresults_atom_group(ag,
                               ftresults,
                               rotations,
                               center=None,
                               out=None):
    """Apply ftresult(s) to an atomgroup, returning a new atomgroup.

    The new atomgroup will have one coordinate set for each ftresult passed.
    ftresult can either be a single ftresult object, or an array of ftresult
    objects."""
    orig_coords = ag.getCoords()  # This returns a copy so we can mutate it
    if center is None:
        center = np.mean(orig_coords, axis=0)
    np.subtract(orig_coords, center, orig_coords)

    try:
        if len(ftresults.shape) == 0:
            ftresults = np.expand_dims(ftresults, 0)
    except:
        raise ValueError("ftresults does not seem to be an ndarray or void")

    if out is None:
        out = ag.copy()

    new_coords = np.dot(rotations[ftresults['roti']],
                        orig_coords.T).transpose(0, 2, 1)
    np.add(new_coords, np.expand_dims(ftresults['tv'] + center, 1), new_coords)
    out._setCoords(new_coords, overwrite=True)

    return out


signal(SIGPIPE, SIG_DFL)

def apply_results_main(limit, index, rotation, out_prefix, ftfile, rotations, pdb_file):

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

        for i in range(new_ag.numCoordsets()):
            new_ag.setACSIndex(i)
            writePDB("{}.{}.pdb".format(out_prefix, i), new_ag, csets=i)