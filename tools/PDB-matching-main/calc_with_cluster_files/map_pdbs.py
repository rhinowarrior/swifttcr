"""
Name: map_pdbs.py
Function: Map the residues between reference and decoy PDBs, then save them.
Date: 16-10-2024
Author: Jan Aarts, Farzaneh Meimandi Parizi, Nils Smit
"""

"""
Input:
    - pipeline_dir: Directory where decoys are stored (The directory where the target structures are stored)
    - reference_dir: Directory where reference models are stored
    - output_dir: Output directory to save the mapped files
    - cluster_input: Name of the cluster file example: "clustering.txt"
    
Example: python map_pdbs.py pipeline_dir reference_dir output_dir cluster_input
"""

from Bio.PDB import PDBParser, PDBIO, Select
from pathlib import Path
from Bio import pairwise2
from Bio.SeqUtils import seq1
import sys
import multiprocessing as mp

def main():
    pipeline_dir = sys.argv[1]  # Directory where decoys are stored
    reference_dir = sys.argv[2] # Directory where reference models are stored
    output_dir = sys.argv[3]    # Output directory to save the mapped files
    cluster_input = sys.argv[4] # Name of the cluster file example: "clustering.txt"
    num_threads = int(sys.argv[5])  # Number of threads (cores) to use
    
    max_threads = max(1, num_threads)
    
    if not Path(output_dir).exists():
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Create a dictionary with the paths to all the decoy models that are contained in the clustering file
    model_dict = create_model_clus_dict(pipeline_dir, cluster_input)

    # Prepare tasks for multiprocessing
    tasks = []
    for model_identifier, decoys in model_dict.items():
        # Create the path to the reference model
        reference = str(Path(reference_dir, model_identifier + "_renumb.pdb"))

        for decoy in decoys:
            # Create the path to the decoy model
            val_path = str(Path(pipeline_dir, model_identifier, "merged", decoy))
            # Append task (function arguments) to the tasks list
            tasks.append((reference, val_path, output_dir, model_identifier))

    # Use multiprocessing with Pool
    with mp.Pool(processes=max_threads) as pool:
        # Map the tasks to the pool of processes
        pool.starmap(map_and_save_pdbs, tasks)
                
            
def create_model_clus_dict(pipeline_dir, cluster_input="clustering.txt"):
    """Parse the clustering file and return a dictionary of models.
    
    Args: 
        p_dir (str): Path to the directory containing folders that contain the clustering files.
        cluster_input (str): Name of the cluster file standard is: "clustering.txt"
    
    Returns:
        model_dict (dict): Dictionary containing the paths to the models
    """
    model_dict = {}
    pipeline = Path(pipeline_dir)
    for model_dir in pipeline.iterdir():
        if  model_dir.name:
            models = parse_clustering_file(str(Path(model_dir, cluster_input)))
            # Get the first four characters of the model directory name which is the model key
            model_identifier = model_dir.name[:4]
            model_dict[model_identifier] = models
    return model_dict

def parse_clustering_file(clustering_path):
    """Parse a clustering file and return a list of model names.
    
    Args:
        clustering_path (str): Path to the clustering file
    
    Returns:
        models (list): List of model names
    """
    models = []
    with open(clustering_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("Cluster center: "):
            model_name = line.split("Cluster center: ")[-1]
            model_name = model_name.split("with")[0].strip()
            models.append(model_name)
    return models

def map_and_save_pdbs(reference, decoy, output_dir, model_identifier):
    """Map the reference and decoy PDBs, then save them.
    
    Args:
        reference (str): Path to the reference PDB file
        decoy (str): Path to the decoy PDB file
        output_dir (str): Path to the output directory
        model_identifier (str): Key to identify the model
    
    Returns:
        pdb_ref (Bio.PDB.Structure.Structure): Mapped reference PDB structure
        pdb_decoy (Bio.PDB.Structure.Structure): Mapped decoy PDB structure
    """
    pdb_ref = PDBParser().get_structure('reference', reference)
    pdb_decoy = PDBParser().get_structure('target', decoy)

    # Extract the number from the merged filename (assuming format is like merged_X.pdb)
    merged_number = Path(decoy).stem.split('_')[-1]

    # Map the reference and decoy PDBs based on the residues and remove non-mapped residues
    pdb_ref, pdb_decoy = map_PDBs(pdb_ref, pdb_decoy)

    # Renumber the residues in both PDBs so they have the same numbering
    pdb_ref = reres(pdb_ref)
    pdb_decoy = reres(pdb_decoy)

    # Save the mapped PDBs with the unique naming format
    ref_outfile = str(Path(output_dir, f"{model_identifier}_reference_{merged_number}_mapped.pdb"))
    decoy_outfile = str(Path(output_dir, f"{model_identifier}_merged_{merged_number}_mapped.pdb"))
    
    write_PDB(pdb_ref, ref_outfile)
    write_PDB(pdb_decoy, decoy_outfile)
    
    return pdb_ref, pdb_decoy


def map_PDBs(pdb_ref, pdb_decoy, remove_non_mapped=True, custom_map={}):
    """Map the residues between reference and decoy PDBs.
    
    Args:
        pdb_ref (Bio.PDB.Structure.Structure): Reference PDB structure
        pdb_decoy (Bio.PDB.Structure.Structure): Decoy PDB structure
        remove_non_mapped (bool): Remove non-mapped residues (default=True)
        custom_map (dict): Custom residue mapping (default={})
    
    Returns:
        pdb_ref (Bio.PDB.Structure.Structure): Mapped reference PDB structure
        pdb_decoy (Bio.PDB.Structure.Structure): Mapped decoy PDB structure
    """
   
    # Extract the sequences from the PDBs and sort them
    ref_sequences = [[chain.id, seq1(''.join([res.resname for res in chain]), custom_map=custom_map)] for chain in pdb_ref.get_chains()]

    ref_sequences.sort()
    
    decoy_sequences = [[chain.id, seq1(''.join([res.resname for res in chain]), custom_map=custom_map)] for chain in pdb_decoy.get_chains()]
                            
    decoy_sequences.sort()
    
    # Check if the sequences contain the same amount of chains if not raise an error
    assert(len(ref_sequences) == len(decoy_sequences))
    
    # align the sequences per chain against each other and update the sequences
    for chain in range(len(ref_sequences)):
            pair = pairwise2.align.globalxx(ref_sequences[chain][1], decoy_sequences[chain][1])[0]
            ref_sequences[chain][1]   = pair.seqA
            decoy_sequences[chain][1] = pair.seqB

    # Extract the residue numbers from the sequences
    ref_sequences = [[seq[0],[i+1 for i,res in enumerate(seq[1]) if res != '-']] for seq in ref_sequences]
    decoy_sequences = [[seq[0],[i+1 for i,res in enumerate(seq[1]) if res != '-']] for seq in decoy_sequences]

    # Extract the gapped residues
    ref_gapped_residues = [[ref_sequences[j][0],list(set(ref_sequences[j][1]) -set(decoy_sequences[j][1])) ] for j in range(len(ref_sequences))]  
    decoy_gapped_residues = [[decoy_sequences[j][0],list(set(decoy_sequences[j][1]) -set(ref_sequences[j][1])) ] for j in range(len(decoy_sequences))]  

    # Create dictionaries from the lists
    ref_sequences = {chain:res for chain, res in ref_sequences}
    decoy_sequences = {chain:res for chain, res in decoy_sequences}
    
    ref_gapped_residues = {chain:res for chain, res in ref_gapped_residues}
    decoy_gapped_residues = {chain:res for chain, res in decoy_gapped_residues}

    # Assign the sequences to the PDBs and remove the non-mapped residues if not specified
    if remove_non_mapped:
        pdb_ref = assign(pdb_ref, ref_sequences,ref_gapped_residues )
        pdb_decoy = assign(pdb_decoy, decoy_sequences,  decoy_gapped_residues)
    else:
        pdb_ref = assign(pdb_ref, ref_sequences )
        pdb_decoy = assign(pdb_decoy, decoy_sequences)
    return pdb_ref, pdb_decoy

def assign(pdb, pdb_sequences, residues_to_remove=None):
    """Renumber the residues based on alignment.
    
    Args:
        pdb (Bio.PDB.Structure.Structure): PDB structure
        pdb_sequences (dict): Dictionary of sequences
        residues_to_remove (dict): Dictionary of residues to remove (default=None)
    
    Returns:
        pdb (Bio.PDB.Structure.Structure): PDB structure with renumbered residues
    """
    for chain in pdb.get_chains():
        seq = pdb_sequences.get(chain.id)
        for ind, res in enumerate(chain):
            res.id = ('X', seq[ind], res.id[2])

    for chain in pdb.get_chains():
        for res in chain:
            res.id = (' ', res.id[1], ' ')

    if residues_to_remove:
        # Remove the selected residues from the PDB
        pdb = remove_selected_residues(pdb, residues_to_remove)

    return pdb

def remove_selected_residues(pdb, residues_to_remove):
    """Remove specific residues from a PDB structure.

    Args:
        pdb (Bio.PDB.Structure.Structure): PDB structure
        residues_to_remove (dict): Dictionary of residues to remove
    
    Returns:
        pdb (Bio.PDB.Structure.Structure): PDB structure with removed residues
    """
    # Create a list of residues to remove from the corresponding chains
    chains_to_remove = [k for k in residues_to_remove]

    for chain in pdb.get_chains():
        if chain.id in chains_to_remove:
            for residue in chain.get_residues():
                if residue.id[1] in residues_to_remove.get(chain.id):
                    residue.id = ('X', residue.id[1], residue.id[2])

            residues = reversed(list(enumerate(chain)))
            for indx, res in residues:
                if 'X' == res.id[0]:
                    # Remove the residue from the chain
                    chain.detach_child(('X', res.id[1], res.id[2]))
    return pdb

def reres(pdb):
    """Renumber residues in a PDB structure.
    
    Args:
        pdb (Bio.PDB.Structure.Structure): PDB structure
    
    Returns:
        pdb (Bio.PDB.Structure.Structure): PDB structure with renumbered residues
    """
    # Create a dictionary with the new residue numbers
    new_residueNum = [[chain.id, [i + 1 for i, res in enumerate(chain.get_residues())]] for chain in pdb.get_chains()]
    new_residueNum = {chain: res for chain, res in new_residueNum}
    # Assign the new residue numbers to the PDB
    assign(pdb, new_residueNum)
    return pdb

def write_PDB(pdb, file_name):
    """Write a PDB structure to a file.
    
    Args: pdb (Bio.PDB.Structure.Structure): PDB structure
            file_name (str): Name of the output file
    """
    io = PDBIO()
    io.set_structure(pdb)
    io.save(file_name, select=NotDisordered())

class NotDisordered(Select):
    """Class to handle disordered atoms.
    
    Args:
        Select (class): Select class
    
    Returns:
        bool: True if the atom is not disordered, False otherwise
    """
    def accept_atom(self, atom):
        keepAltID = 'A'
        if not atom.is_disordered() or atom.get_altloc() == keepAltID:
            atom.set_altloc(" ")  # Eliminate alt location ID before output.
            return True
        return False

if __name__ == "__main__":
    main()
