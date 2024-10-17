"""
Name: map_all_pdbs.py
Function: Map the residues between reference and all the decoy PDBs, then saves them.
Date: 17-10-2024
Author: Jan Aarts, Farzaneh Meimandi Parizi, Nils Smit
"""

"""
Input:
    - pipeline_dir: Directory where decoys are stored
    - reference_dir: Directory where reference models are stored
    - output_dir: Output directory to save the mapped files
    - num_threads: Number of threads (cores) to use
    
Example: python map_all_pdbs.py pipeline_dir reference_dir output_dir num_threads
"""

from Bio.PDB import PDBParser, PDBIO, Select
from pathlib import Path
from Bio import pairwise2
from Bio.SeqUtils import seq1
import sys
import multiprocessing as mp
import gc

def main():
    pipeline_dir = sys.argv[1]  # Directory where decoys are stored
    reference_dir = sys.argv[2] # Directory where reference models are stored
    output_dir = sys.argv[3]    # Output directory to save the mapped files
    num_threads = int(sys.argv[4])  # Number of threads (cores) to use

    max_threads = max(1, num_threads)

    # Create the output directory if it doesn't exist
    if not Path(output_dir).exists():
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Create a dictionary with the paths to all the decoy models in the pipeline directory
    model_dict = create_model_dict(pipeline_dir)

    # Prepare tasks for multiprocessing
    tasks = []
    for model_identifier, decoys in model_dict.items():
        # Create the path to the reference model
        reference = str(Path(reference_dir, model_identifier + "_renumb.pdb"))

        for decoy in decoys:
            # Append task (function arguments) to the tasks list
            tasks.append((reference, decoy, output_dir, model_identifier))

    # Use multiprocessing with Pool
    with mp.Pool(processes=max_threads) as pool:
        # Map the tasks to the pool of processes with error handling
        results = pool.starmap(map_and_save_pdbs, tasks)

def create_model_dict(pipeline_dir):
    """Return a dictionary of models based on the structure of the directory.
    
    Args: 
        p_dir (str): Path to the directory containing folders that contain the decoy models.
    
    Returns:
        model_dict (dict): Dictionary containing the paths to the models
    """
    model_dict = {}
    pipeline = Path(pipeline_dir)

    # Iterate over all subdirectories in the pipeline_dir
    for model_dir in pipeline.iterdir():
        if model_dir.is_dir():  # Ensure it's a directory
            # Get all 'merged_X.pdb' files in the "merged" subdirectory
            merged_dir = Path(model_dir, "merged")
            if merged_dir.exists():
                merged_files = list(merged_dir.glob("merged_*.pdb"))
                if merged_files:
                    # The model_identifier is the first 4 characters of the directory name (as per the original logic)
                    model_identifier = model_dir.name[:4]
                    model_dict[model_identifier] = merged_files
    return model_dict

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
    try:
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

    except Exception as e:
        print(f"Error processing {decoy}: {e}")
    finally:
        # Clean up large objects to free memory
        del pdb_ref, pdb_decoy
        gc.collect()  # Trigger garbage collection

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
   
    ref_sequences = [[chain.id, seq1(''.join([res.resname for res in chain]), custom_map=custom_map)] for chain in pdb_ref.get_chains()]
    ref_sequences.sort()
    
    decoy_sequences = [[chain.id, seq1(''.join([res.resname for res in chain]), custom_map=custom_map)] for chain in pdb_decoy.get_chains()]
    decoy_sequences.sort()
    
    assert(len(ref_sequences) == len(decoy_sequences))
    
    for chain in range(len(ref_sequences)):
        pair = pairwise2.align.globalxx(ref_sequences[chain][1], decoy_sequences[chain][1])[0]
        ref_sequences[chain][1]   = pair.seqA
        decoy_sequences[chain][1] = pair.seqB

    ref_sequences = [[seq[0],[i+1 for i,res in enumerate(seq[1]) if res != '-']] for seq in ref_sequences]
    decoy_sequences = [[seq[0],[i+1 for i,res in enumerate(seq[1]) if res != '-']] for seq in decoy_sequences]

    ref_gapped_residues = [[ref_sequences[j][0],list(set(ref_sequences[j][1]) -set(decoy_sequences[j][1])) ] for j in range(len(ref_sequences))]  
    decoy_gapped_residues = [[decoy_sequences[j][0],list(set(decoy_sequences[j][1]) -set(ref_sequences[j][1])) ] for j in range(len(decoy_sequences))]  

    ref_sequences = {chain:res for chain, res in ref_sequences}
    decoy_sequences = {chain:res for chain, res in decoy_sequences}
    
    ref_gapped_residues = {chain:res for chain, res in ref_gapped_residues}
    decoy_gapped_residues = {chain:res for chain, res in decoy_gapped_residues}
 
    if remove_non_mapped:
        pdb_ref = assign(pdb_ref, ref_sequences, ref_gapped_residues)
        pdb_decoy = assign(pdb_decoy, decoy_sequences,  decoy_gapped_residues)
    else:
        pdb_ref = assign(pdb_ref, ref_sequences)
        pdb_decoy = assign(pdb_decoy, decoy_sequences)

    return pdb_ref, pdb_decoy

def assign(pdb, pdb_sequences, residues_to_remove=None):
    for chain in pdb.get_chains():
        seq = pdb_sequences.get(chain.id)
        for ind, res in enumerate(chain):
            res.id = ('X', seq[ind], res.id[2])

    for chain in pdb.get_chains():
        for res in chain:
            res.id = (' ', res.id[1], ' ')

    if residues_to_remove:
        pdb = remove_selected_residues(pdb, residues_to_remove)

    return pdb

def remove_selected_residues(pdb, residues_to_remove):
    chains_to_remove = [k for k in residues_to_remove]

    for chain in pdb.get_chains():
        if chain.id in chains_to_remove:
            for residue in chain.get_residues():
                if residue.id[1] in residues_to_remove.get(chain.id):
                    residue.id = ('X', residue.id[1], residue.id[2])

            residues = reversed(list(enumerate(chain)))
            for indx, res in residues:
                if 'X' == res.id[0]:
                    chain.detach_child(res.id)

    return pdb

def reres(pdb):
    """Renumber residues in the PDB file.
    
    Args:
        pdb (Bio.PDB.Structure.Structure): PDB structure
    
    Returns:
        pdb (Bio.PDB.Structure.Structure): Renumbered PDB structure
    """
    for model in pdb.get_models():
        for chain in model.get_chains():
            for i, residue in enumerate(chain.get_residues()):
                residue.id = (' ', i + 1, ' ')  # Renumber residues
    return pdb

def write_PDB(pdb, filename):
    """Write the PDB structure to a file.
    
    Args:
        pdb (Bio.PDB.Structure.Structure): PDB structure
        filename (str): Filename to save the structure
    """
    io = PDBIO()
    io.set_structure(pdb)
    io.save(filename)

if __name__ == "__main__":
    main()
