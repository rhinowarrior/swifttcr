"""
Name: merge_pdbs.py
Function: Script to merge pdb files and rename chains to A for receptor and D for ligand. The receptor (pMHC) has 3 chains A,B and C and the ligand (TCR) has 2 chains D and E.
Date: 2021-07-07
Author: Yannick Aarts
"""

"""
For now had to remove pdb_tidy and the last pdb_merge because these caused the alternative residues to be placed after the normal residues. For now this works but I think it is better to run the pdb-tools scripts as functions so that we don't have to use the commandline and when we make it runnable via python functions we can also change the sorting algorithm so that they are sorted correctly with alternative residues first.
"""

from pathlib import Path
import subprocess
import os
import flip_alternative
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

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
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {command}")
        print(f"Error message: {e.stderr}")
        raise
    
def process_ligand(file, receptor_name, p_out):
    """Process each ligand PDB file to create the merged output.
    
    Args:
        file (str): Path to the ligand file
        receptor_name (str): Name of the modified receptor file
        p_out (str): Path to the output directory
    """
    with NamedTemporaryFile(suffix=".pdb", delete=True) as temp_E, \
        NamedTemporaryFile(suffix=".pdb", delete=True) as temp_D, \
        NamedTemporaryFile(suffix=".pdb", delete=True) as temp_ligand:  # Temporary file for ligand

        # Chain E shift command
        command_shift = f"pdb_selchain -E {file} | pdb_shiftres -2000 | pdb_chain -D > {temp_E.name}"
        # Chain D command
        command_lig = f"pdb_selchain -D {file} > {temp_D.name}"
        # Merge E and D into the temporary ligand file
        command_DE = f"pdb_merge {temp_D.name} {temp_E.name} > {temp_ligand.name}"

        # Run commands
        run_command(command_shift)
        run_command(command_lig)
        run_command(command_DE)

        # Generate merged output name based on number after the last period
        file_number = file.stem.split('.')[-1]  # Extract the number portion
        merged_name = f"merged_{file_number}.pdb"  # Construct the new file name

        # Merge receptor and ligand into a single file
        merge_command = f"cat {receptor_name} {temp_ligand.name} | grep '^ATOM ' > {Path(p_out, merged_name)}"
        run_command(merge_command)
        
def merge_pdbs_main(receptor, ligand, output_dir, num_cores):
    """
    Merge pdb files and rename chains to A for receptor and D for ligand.
    
    Args:
        receptor (str): Path to the receptor pdb file
        ligand (str): Path to the ligand pdb file or directory
        output_dir (str): Path to the output directory
        num_cores (int): Number of cores to use for parallel processing
    """
    p = Path(ligand)
    p_rec = Path(receptor)
    p_out = Path(output_dir)
    os.chdir(p_rec.parent)

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
    receptor_name = f"{p_rec.stem}_rename.pdb"
    command = f"pdb_tidy {p_rec} | pdb_selchain -A,B,C | pdb_chain -A | pdb_reres -1 > {receptor_name}"

    # Process receptor only once
    run_command(command)
    flip_alternative.reorder_residues_in_structure(receptor_name, receptor_name)

    if p.is_dir():
        ligand_files = [f for f in p.iterdir() if f.suffix == ".pdb"]

        # Use multiprocessing pool with optimized tasks
        with Pool(num_cores, maxtasksperchild=10) as pool:
            pool.starmap(process_ligand, [(f, receptor_name, p_out) for f in ligand_files])
    else:
        print("Ligand path is not a directory.")