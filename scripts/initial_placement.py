#!/usr/bin/env python
# Li Xue
# 16-Mar-2023 17:22
"""
Name: initial_placement.py
Function: This script is used to superimpose the target TCR and p-MHC to the reference TCR and p-MHC. The script uses PyMOL to superimpose the structures and gets the alignment from pymol. The script then changes the chainID of the target to the chainID of the reference. The output is the superimposed structures with the chains renamed to that of the reference chains.
date: 25-09-2024
Author: Nils Smit, Li Xue
"""

"""
In the future it is beter to also renumber p-MHC so that this procces can become even faster. 
Todo: 
Add a more secure way to change the chainIDs.
"""

# Have to import cmd because that stops the warning from PyMOL because now i think it can overwrite the cmd module and otherwise pymol will give a warning
import subprocess
import os
import cmd
from pymol import cmd as pymol_cmd
from pathlib import Path


def run_pdb_tools(input_pdb, chain_target, ref_chain):
    """
    Use PDB Tools to select chain and rename chain/segment based on reference chain.

    Args:
        input_pdb (str): The input PDB file path.
        chain_target (str): The chain from the target PDB file that we want to modify.
        ref_chain (str): The reference chain ID that we want to apply to the target PDB file.

    Returns:
        str: The path to the temporary modified PDB file.
    
    Raises:
        RuntimeError: If the PDB Tools command fails.
    """
    temp_output = f"temp_{chain_target}.pdb"
    
    pdb_command = f"pdb_selchain -{chain_target} {input_pdb} | pdb_chain -{ref_chain} | pdb_seg -{ref_chain} | grep '^ATOM' > {temp_output}"
    
    try:
        subprocess.run(pdb_command, shell=True, check=True)

        if os.path.exists(temp_output):
            with open(temp_output, 'r') as temp_file:
                content = temp_file.read()
        else:
            raise RuntimeError(f"Temporary file {temp_output} was not created.")
    
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error occurred while running PDB Tools for chain {chain_target}: {e}")
    
    return temp_output


def merge_pdb_files(temp_files, output):
    """
    Merge the temporary PDB files using pdb_merge and save the final output.

    Args:
        temp_files (list): List of temporary PDB files to merge.
        output (str): Output file path for the final merged PDB file.
    
    Raises:
        RuntimeError: If the pdb_merge command fails.
    """
    temp_files_str = ' '.join(temp_files)
    merge_command = f"pdb_merge {temp_files_str} > {output}"
    
    try:
        subprocess.run(merge_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error occurred while merging PDB files: {e}")


def superpose_and_change_chain_IDs(reference, target, output):
    """
    Superpose the target structure to the reference structure, change the chain IDs using PDB Tools,
    and only superpose the final merged structure.

    Args:
        reference (str): path to the reference PDB file
        target (str): path to the target PDB file
        output (str): path to save the final superposed and chain-modified target PDB file
    
    Raises:
        RuntimeError: If any subprocess or PyMOL operation fails.
    """
    # Initialize PyMOL and clear the workspace
    pymol_cmd.reinitialize()
    
    # Load the reference and target PDB files for initial superposition
    pymol_cmd.load(reference, "ref")
    pymol_cmd.load(target, "target")
    
    # Step 1: Initial Superposition (before chain ID changes)
    pymol_cmd.super("target and name CA", "ref and name CA", object="alignment")

    # Get the raw alignment between reference and target
    alignment = pymol_cmd.get_raw_alignment("alignment")
    
    if alignment:
        ref_atoms = pymol_cmd.get_model("ref").atom
        target_atoms = pymol_cmd.get_model("target").atom
        
        chain_mapping = {}
        temp_files = []

        # Build chain mapping from target to reference based on alignment
        for ref_idx, target_idx in alignment:
            ref_atom = ref_atoms[ref_idx[1] - 1]
            target_atom = target_atoms[target_idx[1] - 1]
            
            ref_chain = ref_atom.chain
            target_chain = target_atom.chain

            if target_chain not in chain_mapping:
                chain_mapping[target_chain] = ref_chain
                print(f"Mapping target chain {target_chain} to reference chain {ref_chain}")

        # Apply PDB Tools to rename chains and save temporary files
        for target_chain, ref_chain in chain_mapping.items():
            try:
                temp_pdb_file = run_pdb_tools(target, target_chain, ref_chain)
                temp_files.append(temp_pdb_file)
            except RuntimeError as e:
                print(e)
                for temp_file in temp_files:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                return

        # Merge the chain-modified temporary files using pdb_merge
        try:
            merge_pdb_files(temp_files, output)
        except RuntimeError as e:
            print(e)
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            return

        # Load and Superpose the final merged structure to reference
        print(f"Loading and superposing the final merged structure: {output}")
        pymol_cmd.load(output, "final_target")
        pymol_cmd.super("final_target and name CA", "ref and name CA")

        # Save the final superposed structure
        pymol_cmd.save(output, "final_target")
        print(f"Final superposed structure saved as: {output}")

        # Clean up temporary files after merging
        for temp_file in temp_files:
            os.remove(temp_file)

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