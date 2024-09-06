import os
from argparse import ArgumentParser

def get_arguments():
    """Gets the arguments from the user using the command line

    Returns:
        args: The arguments from the user
    """
    parser = ArgumentParser(description="SwiftTCR")
    parser.add_argument("--receptor", "-r", required=True, help="Path to receptor pdb file")
    parser.add_argument("--ligand", "-l", required=True, help="Path to ligand pdb file")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--restraints", "-rs", required=True, help="Path to restraints file")
    #parser.add_argument("--rotations", "-ro", required=True, help="Path to rotations file", default="swifttcr/example/input/filtered_cr_in_60.prm")
    parser.add_argument("--outprefix", "-op", required=True, help="Name of the output file")
    
    args = parser.parse_args()
    return args

def check_files(receptor, ligand, output, restraints):
    """Checks if the files exist

    Args:
        receptor (str): Path to receptor pdb file
        ligand (str): Path to ligand pdb file
        output (str): Output directory
        restraints (str): Path to restraints file
    """
    if not os.path.exists(receptor):
        print(f"Receptor file {receptor} does not exist")
        exit(1)
    if not os.path.exists(ligand):
        print(f"Ligand file {ligand} does not exist")
        exit(1)
    if not os.path.exists(restraints):
        print(f"Restraints file {restraints} does not exist")
        exit(1)
    if not os.path.exists(output):
        print(f"Output directory {output} does not exist")
        exit(1)

def check_file_extensions(receptor, ligand, restraints, rotations):
    """Checks if the file extensions are correct
    
    Args: 
        receptor (str): Path to receptor pdb file
        ligand (str): Path to ligand pdb file
        restraints (str): Path to restraints file
        rotations (str): Path to rotations file
    """
    if not receptor.endswith(".pdb"):
        print("Receptor file must be a pdb file")
        exit(1)
    if not ligand.endswith(".pdb"):
        print("Ligand file must be a pdb file")
        exit(1)
    if not restraints.endswith(".json"):
        print("Restraints file must be a json file")
        exit(1)
    if not rotations.endswith(".prm"):
        print("Rotations file must be a prm file")
        exit(1)
