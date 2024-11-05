"""
Title: flip_alternative.py
Function: This script reorders the residues in a PDB file so that alternative residues come before the normal residues. This is done so that the residues are shown correctly in PyMOL.
Date: 10-10-2024
Author: Nils Smit
"""

from Bio import PDB

def reorder_residues_in_structure(pdb_file, output_file):
    """ Reorder the residues in a PDB file according to specified rules for normal and alternative residues.
    
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

            # Dictionary to hold residues by their sequence numbers
            residue_dict = {}
            for residue in residues:
                res_id = residue.get_id()
                resseq = res_id[1]  # Residue sequence number
                insertion_code = res_id[2]  # Insertion code, if any

                # Initialize a list in the dictionary if not already present
                if resseq not in residue_dict:
                    residue_dict[resseq] = []

                # Append the residue to the list for this residue number
                residue_dict[resseq].append(residue)

            # Now we will construct the ordered list of residues
            resseq_keys = sorted(residue_dict.keys(), key=int)

            for idx, resseq in enumerate(resseq_keys):
                res_list = residue_dict[resseq]

                # Separate primary residues (no insertion code) from alternatives
                primary_residues = [res for res in res_list if res.get_id()[2] == ' ']
                alternative_residues = [res for res in res_list if res.get_id()[2] != ' ']

                # Sort alternative residues in descending order by insertion code
                alternative_residues.sort(key=lambda res: res.get_id()[2], reverse=True)

                # If there are alternatives for the next residue
                if idx + 1 < len(resseq_keys):
                    next_resseq = resseq_keys[idx + 1]
                    next_res_list = residue_dict[next_resseq]

                    # Check if the next residue has alternatives
                    next_alternatives = [res for res in next_res_list if res.get_id()[2] != ' ']

                    # If the next residue has alternatives, keep current normal first
                    if next_alternatives:
                        # Add primary residues first, then alternatives
                        reordered_residues.extend(primary_residues)
                        reordered_residues.extend(alternative_residues)
                    else:
                        # If the next residue has no alternatives, add alternatives first
                        reordered_residues.extend(alternative_residues)
                        reordered_residues.extend(primary_residues)
                else:
                    # Last residue case: just add normal and alternative residues
                    reordered_residues.extend(primary_residues)
                    reordered_residues.extend(alternative_residues)

            # Set the reordered residues back to the chain
            chain.child_list = reordered_residues

    # Write the reordered structure to a new PDB file
    io.set_structure(structure)
    io.save(output_file)