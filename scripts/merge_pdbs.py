"""
Name: merge_pdbs.py
Function: Script to merge pdb files and rename chains to A for receptor and D for ligand. The receptor (pMHC) has 3 chains A,B and C and the ligand (TCR) has 2 chains D and E.
Date: 2021-07-07
Author: Jan Aarts
"""

from pathlib import Path
import subprocess
import os


def merge_pdbs_main(receptor, ligand, output_dir):
    """Merge pdb files and rename chains to A for receptor and D for ligand. The receptor (pMHC) has 3 chains A,B and C and the ligand (TCR) has 2 chains D and E.

    Args:
        receptor (str): Path to the receptor pdb file
        ligand (str): Path to the ligand pdb file
        output_dir (str): Path to the output directory
    """
    p = Path(ligand)
    p_rec = Path(receptor)
    p_out = Path(output_dir)
    print("Ligand: ", p)
    print("Outputdir: ", p_out)
    print("Receptor: ", p_rec)
    os.chdir(p_rec.parent)
    receptor_name = "{}".format(p_rec.stem + "_rename.pdb")
    
    #rename receptor
    command = "pdb_tidy {} | pdb_selchain -A,B,C | pdb_chain -A | pdb_reres -1 > {}".format(str(p_rec), receptor_name)
    os.system(command)
    
    if p.is_dir():
        for f in p.iterdir():
            #rename ligand
            #merge receptor and ligand
            #clean up
            if f.suffix ==".pdb":
                print(f.name)
                ligand_name = "{}".format(f.stem + "_rename.pdb")
                output_E = "{}".format(f.stem + "_Eshift.pdb")
                output_D = "{}".format(f.stem + "_Dshift.pdb")
                
                #had to change this because some results became weird if i didn't shift the residues with 2000 inplace of 1000
                command_shift = "pdb_tidy {} | pdb_selchain -E | pdb_shiftres -2000 | pdb_chain -D > {}".format(str(f), output_E)
               
                #select chain D
                command_lig = "pdb_tidy {} | pdb_selchain -D > {}".format(str(f), output_D)
    
                #merge chain D and E
                command_DE = "pdb_merge {} {} | pdb_sort > {}".format(output_D, output_E, ligand_name)
                
                #merge receptor and ligand
                merged_name = "merged_" + str(f).split(".")[1] + ".pdb"
                merge_command = "pdb_merge {} {} | pdb_tidy -strict > {}".format(receptor_name, ligand_name, str(Path(p_out, merged_name)))

                os.system(command_shift)

                os.system(command_lig)

                os.system(command_DE)
                os.system(merge_command)
                
                #remove renamed file
                os.system("rm {}".format(output_D))
                os.system("rm {}".format(output_E))
                os.system("rm {}".format(ligand_name))

    else:
        pass
