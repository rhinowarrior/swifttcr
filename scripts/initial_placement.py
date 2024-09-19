#!/usr/bin/env python
# Li Xue
# 16-Mar-2023 17:22
#REWRITE using 2 reference structures and do superposition take 2bnr_l_u, 2bnr_r_u from input_new/
"""Script to rotate the ligand and receptor based on alpha helices and conserved cys residues.

https://pdb2sql.readthedocs.io/en/latest/_modules/pdb2sql/superpose.html#superpose

Usage: python3 pca.py <receptor> <ligand> <output_dir>
"""
# Have to import cmd because that stops the warning from PyMOL because now i think it can overwrite the cmd module and otherwise pymol will give a warning
import cmd
from pymol import cmd
from pathlib import Path

def initial_placement_main(receptor, ligand, outputdir, reference_receptor, reference_ligand):
    receptor = Path(receptor)
    ligand = Path(ligand)
    outputdir = Path(outputdir)
    reference_receptor = Path(reference_receptor)
    reference_ligand = Path(reference_ligand)

    # Load receptor, ligand, and reference structures into PyMOL
    cmd.load(str(receptor), 'receptor')
    cmd.load(str(ligand), 'ligand')
    cmd.load(str(reference_receptor), 'ref_receptor')
    cmd.load(str(reference_ligand), 'ref_ligand')
    
    # Superpose ligand to reference ligand (align entire structures)
    results_ligand = cmd.super('receptor', 'ref_receptor')
    final_rmsd_receptor = results_ligand[0]
    print(f"Final RMSD receptor: {final_rmsd_receptor}")
    
    chains_ligand = cmd.get_chains('ligand')
    chains_ref_ligand = cmd.get_chains('ref_ligand')
    
    best_chain_mapping = {}

    # Loop through each chain in pdb1 and try to superpose it on each chain in pdb2
    for chain_pdb1 in chains_ligand:
        best_rmsd = float('inf')
        best_chain_pdb2 = None
        
        for chain_pdb2 in chains_ref_ligand:
            # Superpose each chain individually
            selection1 = f"ligand and chain {chain_pdb1}"
            selection2 = f"ref_ligand and chain {chain_pdb2}"
            try:
                rmsd = cmd.super(selection1, selection2)[0]
                # print(f"RMSD between ligand chain {chain_pdb1} and ref_ligand chain {chain_pdb2}: {rmsd:.2f}")
                
                # If this superposition gives a better (lower) RMSD, update the best chain match
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_chain_pdb2 = chain_pdb2
            except:
                print(f"Failed to superpose ligand chain {chain_pdb1} with ref_ligand chain {chain_pdb2}")

        # Store the best matching chain
        if best_chain_pdb2:
            #print(f"Best match for ligand chain {chain_pdb1} is ref_ligand chain {best_chain_pdb2} with RMSD {best_rmsd:.2f}")
            best_chain_mapping[chain_pdb1] = best_chain_pdb2

    # Now that we have the best chain mappings, rename chains in ligand to match the reference ligand
    for chain_pdb1, chain_pdb2 in best_chain_mapping.items():
        if chain_pdb1 != chain_pdb2:
            #print(f"Renaming chain {chain_pdb1} in ligand to {chain_pdb2}")
            cmd.alter(f"(ligand and chain {chain_pdb1})", f"chain='{chain_pdb2}'")

    # Perform superposition based on the selected regions
    results = cmd.super('ligand', 'ref_ligand')
    final_rmsd_ligand = results[0]
    print(f"Final RMSD ligand: {final_rmsd_ligand}")

    # Check if the final RMSD exceeds the threshold and issue a warning
    if final_rmsd_receptor > 2:
        print(f"Warning: Final RMSD ({final_rmsd_receptor:.2f}) exceeds the threshold of 2")
    if final_rmsd_ligand > 2:
        print(f"Warning: Final RMSD ligand ({final_rmsd_ligand:.2f}) exceeds the threshold of 2")

    # Save the updated structures to the output directory
    output_receptor_path = Path(outputdir, receptor.name)
    output_ligand_path = Path(outputdir, ligand.name)

    cmd.save(str(output_receptor_path), 'receptor')
    cmd.save(str(output_ligand_path), 'ligand')
    
    # Clean up and quit PyMOL
    cmd.delete('all')