"""
Name: merge_pdbs.py
Function: Script to merge pdb files and rename chains to A for receptor and D for ligand. The receptor (pMHC) has 3 chains A,B and C and the ligand (TCR) has 2 chains D and E.
Date: 2021-07-07
Author: Jan Aarts
"""

from pathlib import Path
import subprocess
import os

def run_command(command):
    """
    Helper function to run a command with subprocess.run and handle errors.
    
    Args:
        command (str): The shell command to be executed.
        
    Returns:
        str: The standard output of the command if successful.
    
    Raises:
        subprocess.CalledProcessError: If the command fails.
    """
    try:
        # Run the command and wait for it to complete
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        
        # Return the standard output
        return result.stdout

    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {command}")
        print(f"Error message: {e.stderr}")
        raise  # Re-raise the exception after logging it
    

def merge_pdbs_main(receptor, ligand, output_dir):
    """Merge pdb files and rename chains to A for receptor and D for ligand. The receptor (pMHC) has 3 chains A, B, and C, and the ligand (TCR) has 2 chains D and E.

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

    # This command makes it so the peptide chain is seperated from chain A in the receptor because the peptide chain with the current code made part of the mhc and it tries to connect to it in pymol.
    # This command doesn't work because in pairwise_rmsd.py it tries to pad the chains which crashes the programm so this will work once we replace gradpose with our own code.
    # command = (
    # "pdb_tidy {} | "  # Clean the PDB
    # "pdb_selchain -A,B | pdb_chain -A | pdb_reres -1000 > mhc.pdb; "  # Select chains A and B, rename to A, renumber residues starting from 1000, and save as mhc.pdb
    # "pdb_tidy {} | "  # Clean the PDB again
    # "pdb_selchain -C | pdb_chain -A | pdb_reres -1 > pep.pdb; "  # Select chain C, rename to A, renumber starting from 1, and save as pep.pdb
    # "pdb_merge pep.pdb mhc.pdb | pdb_tidy > {}; "  # Properly concatenate mhc.pdb and pep.pdb into pMHC.pdb, appending to the output
    # "rm mhc.pdb pep.pdb"  # Remove intermediate files mhc.pdb and pep.pdb
    # ).format(str(p_rec), str(p_rec), receptor_name)
    
    # Command for receptor
    command = f"pdb_tidy {p_rec} | pdb_selchain -A,B,C | pdb_chain -A | pdb_reres -1  > {receptor_name}"
    
    # Run receptor command
    run_command(command)

    if p.is_dir():
        for f in p.iterdir():
            if f.suffix == ".pdb":
                ligand_name = "{}".format(f.stem + "_rename.pdb")
                output_E = "{}".format(f.stem + "_Eshift.pdb")
                output_D = "{}".format(f.stem + "_Dshift.pdb")

                # Shift chain E by 2000 residues and rename to chain D
                command_shift = f"pdb_tidy {f} | pdb_selchain -E | pdb_shiftres -2000 | pdb_chain -D > {output_E}"
                
                # Select chain D
                command_lig = f"pdb_tidy {f} | pdb_selchain -D > {output_D}"
    
                # Merge chain D and E
                command_DE = f"pdb_merge {output_D} {output_E} | pdb_sort > {ligand_name}"
                
                # Merge receptor and ligand
                merged_name = "merged_" + str(f).split(".")[1] + ".pdb"
                merge_command = f"pdb_merge {receptor_name} {ligand_name} | pdb_tidy -strict > {Path(p_out, merged_name)}"

                # Run all commands in sequence
                run_command(command_shift)
                run_command(command_lig)
                run_command(command_DE)
                run_command(merge_command)
                
                # Remove intermediate files
                os.remove(output_D)
                os.remove(output_E)
                os.remove(ligand_name)

    else:
        print("Ligand path is not a directory.")    
