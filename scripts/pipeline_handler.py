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
    # parser.add_argument("--restraints", "-rs", required=True, help="Path to restraints file")
    # parser.add_argument("--rotations", "-ro", required=True, help="Path to rotations file", default="example/input/filtered_cr_in_60.prm")
    parser.add_argument("--outprefix", "-op", required=True, help="Name of the output file")
    # parser.add_argument("--chains", "-c", required=False, help="Chains to use in the pdb files, the chains should be added in the following order 1: alpha p-MHC chain, 2: beta alpha p-MHC chain, 3: Peptide chain, 4: alpha TCR chain, 5: beta TCR chain", default=["A", "B", "C", "D", "E"], nargs="*")
    # parser.add_argument("--variabledomain", "-vd", required=True, help="Variable domain of the TCR. fill in like 1-128", type=str)
    # parser.add_argument(
    # '--attractive_res', '-ar',
    # type=str,
    # required=True,
    # help="JSON string representing attractive regions. Example: '{\"D\": {\"start\": [26, 55, 104], \"end\": [39, 66, 118]}, \"E\": {\"start\": [26, 55, 104], \"end\": [39, 66, 118]}, \"C\": {\"start\": [-1], \"end\": [1000]}}'")
    
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


def check_amount_of_chains_pdb(receptor, ligand):
    """Checks if the amount of chains are as excepected in the pdb files
    
    Args:
        receptor (str): Path to receptor pdb file
        ligand (str): Path to ligand pdb file
    """
    with open(receptor, "r") as f:
        count_chains = set()
        for line in f:
            if line.startswith("ATOM"):
                count_chains.add(line[21])
        if len(count_chains) != 3:
            print("peptide-MHC file must have 3 chains")
            exit(1)
    with open(ligand, "r") as f:
        count_chains = set()
        for line in f:
            if line.startswith("ATOM"):
                count_chains.add(line[21])
        if len(count_chains) != 2:
            print("TCR file must have 2 chains")
            exit(1)