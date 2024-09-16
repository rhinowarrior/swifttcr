#!/usr/bin/env python
# Li Xue
# 16-Mar-2023 17:22
#REWRITE using 2 reference structures and do superposition take 2bnr_l_u, 2bnr_r_u from input_new/
"""Script to rotate the ligand and receptor based on alpha helices and conserved cys residues.

https://pdb2sql.readthedocs.io/en/latest/_modules/pdb2sql/superpose.html#superpose

Usage: python3 pca.py <receptor> <ligand> <output_dir>
"""
from pymol import cmd
from sys import argv
import sys
from sklearn.decomposition import PCA
import numpy as np
# from pdb2sql import pdb2sql
# from pdb2sql import transform
# from pdb2sql import superpose
import os
# from pdb2sql import align
from pathlib import Path

def initial_placement_main(receptor, ligand, outputdir, reference_receptor, reference_ligand):
    receptor = Path(receptor)
    ligand = Path(ligand)
    outputdir = Path(outputdir)
    reference_receptor = Path(reference_receptor)
    reference_ligand = Path(reference_ligand)
    
    cmd.load(receptor, "receptor")
    cmd.load(reference_receptor, "ref_receptor")
    
    cmd.super("receptor", "ref_receptor")
    cmd.save(os.path.join(outputdir, receptor.name), "receptor")
    cmd.delete("all")

    cmd.load(ligand, "ligand")
    cmd.load(reference_ligand, "ref_ligand")
    
    cmd.super("ligand", "ref_ligand")
    cmd.save(os.path.join(outputdir, ligand.name), "ligand")
    cmd.delete("all")

# def initial_placement_main(receptor, ligand, outputdir, reference_receptor, reference_ligand, chains, variable_domain):
#     receptor = Path(receptor)
#     ligand = Path(ligand)
#     outputdir = Path(outputdir)
#     reference_receptor = Path(reference_receptor)
#     reference_ligand = Path(reference_ligand)
#     chains = chains
#     start, end = map(int, variable_domain.split('-'))
#     print("Receptor: ", receptor)
#     print("Ligand: ", ligand)

#     db_mhc = pdb2sql(str(receptor))
#     db_tcr = pdb2sql(str(ligand))
#     db_ref_mhc = pdb2sql(str(reference_receptor))
#     db_ref_tcr = pdb2sql(str(reference_ligand))

#     updated_mhc = superpose(db_mhc, db_ref_mhc, export = False)
#     updated_tcr = superpose(db_tcr, db_ref_tcr, export = False, chainID = [chains[3], chains[4]], resSeq = list(range(start,end)))


    # updated_mhc.exportpdb(Path(outputdir, receptor.name))
    # updated_tcr.exportpdb(Path(outputdir, ligand.name))