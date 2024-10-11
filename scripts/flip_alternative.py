"""
Title: flip_alternative.py
Function: This script reorders the residues in a PDB file so that alternative residues come before the normal residues. This is done so that the residues are shown correctly in PyMOL.
Date: 10-10-2024
Author: Nils Smit
"""

from Bio import PDB

def reorder_residues_in_structure(pdb_file, output_file):
    """ Reorder the residues in a PDB file so that alternative residues come before the normal residues.
    
    Args:
        pdb_file (str): Path to the input PDB file.
        output_file (str): Path to the output PDB file with reordered residues.
    """
    # Initialize the PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', pdb_file)
    
    # Create a new structure that will hold the reordered chains
    io = PDB.PDBIO()
    
    # Loop over all models, chains, and residues in the structure
    for model in structure:
        for chain in model:
            residues = list(chain)
            reordered_residues = []

            # Dictionary to hold the base residues and their alternatives
            residue_dict = {}
            for residue in residues:
                # Tuple (hetfield, resseq, icode)
                res_id = residue.get_id()
                base_residue = str(res_id[1])
                # If there is a alternative residue
                if res_id[2] != ' ': 
                    base_residue += res_id[2]

                # Store residue based on whether it's an alternative or normal
                if base_residue.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ') not in residue_dict:
                    residue_dict[base_residue.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ')] = []

                residue_dict[base_residue.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ')].append(residue)

            # Reorder the residues so that alternatives come first
            for base_residue, res_list in residue_dict.items():
                # Sort by alternative residues coming first (sorted by insertion code)
                sorted_residues = sorted(res_list, key=lambda res: res.get_id()[2], reverse=True)
                reordered_residues.extend(sorted_residues)

            # Clear the original chain and add the reordered residues
            chain.child_list = reordered_residues

    # Write the reordered structure to a new PDB file
    io.set_structure(structure)
    io.save(output_file)