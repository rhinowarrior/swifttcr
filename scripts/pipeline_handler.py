"""
Name: pipeline_handler.py
Function: This script is used to handle the arguments from the user. It checks if the files exist, if the file extensions are correct and if the amount of chains in the pdb files are as expected. The output is the arguments from the user.
Date: 25-09-2024
Author: Nils Smit
"""
import os
from argparse import ArgumentParser


def get_arguments():
    """Gets the arguments from the user using the command line

    Returns:
        args: The arguments from the user
    """
    parser = ArgumentParser(description="SwiftTCR")
    parser.add_argument("--pmhc", "-r", required=True, help="Path to mhc pdb file")
    parser.add_argument("--tcr", "-l", required=True, help="Path to tcr pdb file")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--outprefix", "-op", required=True, help="Name of the output file")
    args = parser.parse_args()
    return args


def check_files(receptor, ligand):
    """Checks if the files exist

    Args:
        receptor (str): Path to receptor pdb file
        ligand (str): Path to ligand pdb file
    """
    if not os.path.exists(receptor):
        print(f"Receptor file {receptor} does not exist")
        exit(1)
    if not os.path.exists(ligand):
        print(f"Ligand file {ligand} does not exist")
        exit(1)


def check_file_extensions(receptor, ligand):
    """Checks if the file extensions are correct
    
    Args: 
        receptor (str): Path to receptor pdb file
        ligand (str): Path to ligand pdb file
    """
    if not receptor.endswith(".pdb"):
        print("Receptor file must be a pdb file")
        exit(1)
    if not ligand.endswith(".pdb"):
        print("Ligand file must be a pdb file")
        exit(1)

def check_amount_of_chains_pdb(receptor, ligand):
    """Checks if the amount of chains are as excepected in the pdb files
    
    Args:
        receptor (str): Path to receptor pdb file
        ligand (str): Path to ligand pdb file
    """
    with open(receptor, "r") as f:
        # count the amount of unique chains in the receptor file
        count_chains = set()
        for line in f:
            if line.startswith("ATOM"):
                count_chains.add(line[21])
        if len(count_chains) != 3:
            print("peptide-MHC file must have 3 chains")
            exit(1)
    with open(ligand, "r") as f:
        # count the amount of unique chains in the ligand file
        count_chains = set()
        for line in f:
            if line.startswith("ATOM"):
                count_chains.add(line[21])
        if len(count_chains) != 2:
            print("TCR file must have 2 chains")
            exit(1)