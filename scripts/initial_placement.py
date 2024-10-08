#!/usr/bin/env python
# Li Xue
# 16-Mar-2023 17:22
"""
Name: initial_placement.py
Function: This script is used to superimpose the target TCR and p-MHC to the reference TCR and p-MHC. The script uses PyMOL to superimpose the structures and then finds the closest residue to the reference residues in the target structures. The script then changes the chain of the closest residues to the chain of the reference residues. The output is the superimposed structures with the chains renamed to the reference chains.
date: 25-09-2024
Author: Nils Smit, Li Xue
"""

"""
In the future it is beter to renumber the pMHC to IMGT numbering so that we can look at only a specific residue in both the reference and target. This will make it faster to find which chain is superimposed to which chain in the reference.

Todo: rewrite script summary
"""

# Have to import cmd because that stops the warning from PyMOL because now i think it can overwrite the cmd module and otherwise pymol will give a warning
import cmd
from pymol import cmd as pymol_cmd
from pathlib import Path

 
def superpose_and_change_chain_IDs(reference, target, output):
    """ Superpose the target structure to the reference structure and change the chain IDs based on the alignment.
    
    Args:
        reference (str): path to the reference PDB file
        target (str): path to the target PDB file
        output (str): path to save the superposed target PDB file
    """
    # Initialize PyMOL and clear the workspace
    pymol_cmd.reinitialize()
    
    # Load the reference and target PDB files
    pymol_cmd.load(reference, "ref")
    pymol_cmd.load(target, "target")
    
    # Use super instead of align
    pymol_cmd.super("target and name CA", "ref and name CA", object="alignment")

    # Get the raw alignment
    alignment = pymol_cmd.get_raw_alignment("alignment")
    
    if alignment:
        ref_atoms = pymol_cmd.get_model("ref").atom 
        target_atoms = pymol_cmd.get_model("target").atom
        
        chain_mapping = {}
        
        for ref_idx, target_idx in alignment:
            ref_atom = ref_atoms[ref_idx[1] - 1]  # ref_idx is a tuple (object index, atom index)
            target_atom = target_atoms[target_idx[1] - 1]

            ref_chain = ref_atom.chain
            target_chain = target_atom.chain

            # Create a chain mapping based on the alignment
            if target_chain not in chain_mapping:
                chain_mapping[target_chain] = ref_chain
                print(f"Mapping target chain {target_chain} to reference chain {ref_chain}")
        
        # Separate chains into different objects, rename, and then merge
        for target_chain, ref_chain in chain_mapping.items():
            # Select the target chain and create a new object
            pymol_cmd.create(f"target_chain_{target_chain}", f"target and chain {target_chain}")
            pymol_cmd.alter(f"target_chain_{target_chain}", f"chain='{ref_chain}'")
        
        # Combine the renamed chains back into a single object
        pymol_cmd.create("renamed_target", "target_chain_*")  # Use wildcard to combine all separate chains
        
        # Save the newly created object
        pymol_cmd.save(output, "renamed_target")
        pymol_cmd.remove("all")
    
    else:
        print("No alignment was produced. Please check your input files.")


def initial_placement_main(receptor, ligand, outputdir, reference_receptor, reference_ligand):
    """ Superimposes the target stuctures to reference structures and renames the chains to the reference chains.

    Args:
        receptor (str): The path to the target p-MHC structure
        ligand (str): The path to the target TCR structure
        outputdir (str): The path to the output directory
        reference_receptor (str): The path to the reference p-MHC structure
        reference_ligand (str): The path to the reference TCR structure
    """
    receptor = Path(receptor)
    ligand = Path(ligand)
    outputdir = Path(outputdir)
    reference_receptor = Path(reference_receptor)
    reference_ligand = Path(reference_ligand)
    
    output_receptor_path = Path(outputdir, receptor.name)
    output_ligand_path = Path(outputdir, ligand.name)

    # Superpose the target p-MHC to the reference p-MHC
    superpose_and_change_chain_IDs(reference_receptor, receptor, output_receptor_path)
    
    # Superpose the target TCR to the reference TCR
    superpose_and_change_chain_IDs(reference_ligand, ligand, output_ligand_path)