#!/usr/bin/env python

"""
Usage: Can be used to renumber a bound TCR-p-MHC complex to match the numbering of the p-MHC or can be use to renumber a p-MHC to the numbering of a bound TCR-p-MHC.

"""

from typing import Dict, List
import sys
import os
from uuid import uuid4
from tempfile import gettempdir

from Bio.PDB.Residue import Residue
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB as PDB
from Bio.PDB import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionException
from pymol import cmd as pymol_cmd

# PDB parser from Biopython
pdb_parser = PDBParser()

# Directory for storing temporary files
tmp_dir = gettempdir()

def align(mobile_path: str, reference_path: str) -> str:
    """
    Calls PyMOL to align the input structure to the reference structure.

    Args:
        mobile_path: Path to the mobile PDB file.
        reference_path: Path to the reference PDB file (bound PMHC).
    Returns:
        Path to the alignment file in clustal format.
    """

    # Where to store the alignment
    alignment_name = uuid4()
    alignment_path = os.path.join(tmp_dir, f"{alignment_name}.aln")

    # Init PyMOL
    pymol_cmd.reinitialize()

    # Load mobile and reference structures
    pymol_cmd.load(mobile_path, 'mobile')
    pymol_cmd.load(reference_path, 'reference')

    # Align the structures
    # Might want to use super for this not suer if it is better
    r = pymol_cmd.align("mobile", "reference", object="alignment")
    if r[1] == 0:
        raise ValueError("No residues aligned")

    # Save the alignment
    pymol_cmd.save(alignment_path, selection="alignment", format="aln")

    # Clean up PyMOL
    pymol_cmd.remove("all")

    # Return the alignment file path (string)
    return alignment_path


def parse_clustal(path: str) -> Dict[str, str]:
    """
    Parses clustal alignment format file.

    Args:
        path: Path to alignment file.
    Returns:
        Alignment dictionary: {id 1: aligned sequence 1, id 2: aligned sequence 2, ...}
    """

    alignment = {}
    with open(path, 'rt') as f:
        for line in f:
            # Line format: key and sequence separated by whitespaces
            i = line.find(' ')
            k = line[:i].strip()
            v = line[i:].strip()

            if len(k) > 0 and len(v) > 0:
                alignment[k] = alignment.get(k, '') + v

    return alignment


def list_aligned_residues(alignment_path: str, pdb_path: str) -> List[Residue]:
    """
    List all residues that were aligned to the reference structure based on PyMOL's raw alignment.

    Args:
        alignment_path: Path to the alignment file in clustal format.
        pdb_path: Path to the mobile PDB file.
    Returns:
        List of aligned Residue objects from the mobile structure.
    """

    # Get residues from pdb structure
    structure = pdb_parser.get_structure('mobile', pdb_path)
    model = list(structure.get_models())[0]
    residues = []
    for chain in sorted(model.get_chains(), key=lambda c: c.get_id()):
        residues += list(chain.get_residues())

    # Get reference and mobile aligned sequences
    alignment = parse_clustal(alignment_path)
    aligned_reference = alignment['reference']
    aligned_mobile = alignment['mobile']

    # Find all aligned residues
    aligned_indices = []
    for i in range(len(aligned_reference)):
        if aligned_reference[i] != '-':
            aligned_indices.append(len(aligned_mobile[:i].replace('-', '')))

    # Isolate aligned residues
    aligned_residues = [residues[i] for i in aligned_indices]  # Fix this line

    return aligned_residues


def renumber_residues(pdb_path: str, aligned_residues: List[PDB.Residue.Residue], output_path: str):
    """
    Renumber residues in the PDB structure and save to a new file.
    Only renumber residues corresponding to ATOM records (not HETATM).

    Args:
        pdb_path: Path to the original PDB file.
        aligned_residues: List of aligned Residue objects to be renumbered.
        output_path: Path to save the renumbered PDB file.
    """

    # Load the structure to renumber
    structure = pdb_parser.get_structure('mobile', pdb_path)

    # Loop over the aligned residues and renumber them
    for new_index, residue in enumerate(aligned_residues, start=1):
        # Ensure the residue only contains ATOM records, not HETATM
        if all(atom.get_fullname().strip() == "ATOM" for atom in residue.get_atoms()):
            # Get the parent chain of the residue
            chain = residue.get_parent()
            residue_id = residue.get_id()

            # Save the old residue id
            old_residue_id = residue_id

            # Remove the residue from the chain's dictionary with the old id
            del chain.child_dict[old_residue_id]

            # Update the residue number (create new ID with updated index)
            new_residue_id = (residue_id[0], new_index, residue_id[2])

            # Assign the new ID to the residue
            residue.id = new_residue_id

            # Add the residue back to the chain with the new id
            chain.child_dict[new_residue_id] = residue

    # Save the modified structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)
    
    
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: script.py <mobile_pdb_path> <reference_pmhc_path> <output_pdb_path>")
        sys.exit(1)

    mobile_pdb_path = sys.argv[1]
    reference_pmhc_path = sys.argv[2]
    output_pdb_path = sys.argv[3]

    alignment_path = align(mobile_pdb_path, reference_pmhc_path)
    residues = list_aligned_residues(alignment_path, mobile_pdb_path)

    # Renumber residues and save to the output PDB file
    renumber_residues(mobile_pdb_path, residues, output_pdb_path)

    # Print confirmation
    print(f"Renumbered {mobile_pdb_path} -> saved as {output_pdb_path}")
