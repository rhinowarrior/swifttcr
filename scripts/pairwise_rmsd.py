"""
Name: pairwise_rmsd.py
Function: This script calculates the pairwise RMSD between PDB files. The script can calculate the RMSD between all residues of two chains or only the interface residues of two chains. The output is a list of tuples with the RMSD values and a list of the file paths to the PDB files.
Date: 25-09-2024
Author: Kevin van Geemen
"""
import os
import time
import tempfile
import itertools
import glob
from itertools import repeat
import torch
import torch.multiprocessing as mp
import gradpose
from tqdm import tqdm
import re

# Setting these to stop issues on the cluster
mp.set_sharing_strategy('file_system')
mp.set_start_method("spawn", force=True)

# Helper Functions
def load_ligand_ca_xyzr(pdb_path, chain):
    """Loads CA atoms of a chain from a PDB.
    
    Args:
        pdb_path (str): Path to the PDB file.
        chain (str): Chain to load.
    
    Returns:
        torch.Tensor: Tensor with the CA atoms and residue numbers.
    """
    with open(pdb_path, encoding='utf-8') as pdb_file:
        return torch.tensor(
            [
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                    int(line[22:26])
                ]
                for line in pdb_file.read().split('\n') \
                if line.startswith('ATOM ') and line[13:15] == 'CA' \
                    and line[21] == chain
            ]
        )

def load_ligand_ca_xyz(pdb_path, chain, residues, residues_dict):
    """Loads specified residues and chain's CA atoms.
    
    Args:
        pdb_path (str): Path to the PDB file.
        chain (str): Chain to load.
        residues (list): List of residues to load.
        residues_dict (dict): Dictionary of residue numbers to indices.
    
    Returns:
        torch.Tensor: Tensor with XYZ coordinates of the CA atoms.
    """
    with open(pdb_path, 'r', encoding='utf-8') as pdb_file:
        pxyz = torch.tensor(
            [
                [
                    residues_dict[int(line[22:26])],
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54])
                ]
                for line in pdb_file.read().split('\n') \
                if line.startswith('ATOM ') and line[13:15] == 'CA' \
                    and line[21] == chain \
                    and int(line[22:26]) in residues
            ]
        )
    xyz = torch.zeros((len(residues), 3))
    xyz[pxyz[:, 0].long()] = pxyz[:, 1:]
    return xyz

def calc_selection_rmsd(vars):
    """Calculates the RMSD between two selections of atoms.

    Args:
        vars (tuple): Tuple with the following variables:
            template_index (int): Index of the template PDB.
            xyz_all (torch.Tensor): Tensor with the XYZ coordinates of all PDBs.
            del_mask (torch.Tensor): Tensor with the deletion mask of all PDBs.

    Returns:
        list: List of tuples with the RMSD values.
    """
    return _calc_selection_rmsd(*vars)

def _calc_selection_rmsd(template_index, xyz_all, del_mask):
    """Calculate RMSDs compared to one template.
    Only looks at PDBs after the template index to avoid duplicates.
    
    Args:
        template_index (int): Index of the template PDB.
        xyz_all (torch.Tensor): Tensor with the XYZ coordinates of all PDBs.
        del_mask (torch.Tensor): Tensor with the deletion mask of all PDBs.
    
    Returns:
        list: List of tuples with the RMSD values.
    """
    rmsd_list = []
    xyz_all_current = xyz_all[template_index:]
    del_mask_current = del_mask[template_index:]

    pdb_range = torch.arange(start=1, end=len(xyz_all_current))
    del_mask_current = del_mask_current[:1] * del_mask_current[1:]

    rmsds = ((xyz_all_current[0] - xyz_all_current[pdb_range]) * del_mask_current).pow(2).reshape(len(xyz_all_current)-1, -1).sum(1).div(del_mask_current.sum(1).squeeze()).sqrt()
    for i in range(1, len(xyz_all_current)):
        rmsd_list.append((template_index, i + template_index, rmsds[i-1].item()))
    return rmsd_list

def _load_ligand_ca_xyzr(vars):
    """Wrapper function for loading ligand CA atoms with residue numbers.

    Args:
        vars (tuple): Tuple with the following variables:
            pdb_path (str): Path to the PDB file.
            chain (str): Chain to load.
    
    Returns:
        torch.Tensor: Tensor with the CA atoms and residue numbers.
    """
    return load_ligand_ca_xyzr(*vars)

def _load_ligand_ca_xyz(vars):
    """Wrapper function for loading ligand CA atoms with residue numbers.

    Args:
        vars (tuple): Tuple with the following variables:
            pdb_path (str): Path to the PDB file.
            chain (str): Chain to load.
            residues (list): List of residues to load.
            residues_dict (dict): Dictionary of residue numbers to indices.

    Returns:
        torch.Tensor: Tensor with XYZ coordinates of the CA atoms.
    """    
    return load_ligand_ca_xyz(*vars)

def _pad_chain(vars):
    """Wrapper function for padding missing residues in a residue range.

    Args:
        vars (tuple): Tuple with the following variables:
            xyzr (torch.Tensor): Tensor with the CA atoms and residue numbers.
            low (int): Lowest residue number.
            high (int): Highest residue number.
            
    Returns:
        torch.Tensor: Tensor with the CA atoms and residue numbers, padded with empty atoms.
    """
    return pad_chain(*vars)

def pad_chain(xyzr, low, high):
    """Fill out missing residues in a residue range with empty atoms.
    
    Args:
        xyzr (torch.Tensor): Tensor with the CA atoms and residue numbers.
        low (int): Lowest residue number.
        high (int): Highest residue number.
    
    Returns:
        torch.Tensor: Tensor with the CA atoms and residue numbers, padded with empty atoms.
    """
    # Get all unique residue numbers from the xyzr tensor
    existing_residues = set(xyzr[:, -1].tolist())
    
    # Identify missing residues in the range
    missing_residues = [res for res in range(low, high + 1) if float(res) not in existing_residues]
    
    # Create a tensor for missing residues with zero coordinates
    if missing_residues:
        padding = torch.zeros((len(missing_residues), 4))
        padding[:, -1] = torch.tensor(missing_residues, dtype=torch.float)
        
        # Concatenate existing xyzr with the padding and sort only once
        xyzr = torch.cat((xyzr, padding), dim=0)
    
    # Sort by residue number in the last column
    sorted_indices = xyzr[:, -1].argsort()
    return xyzr[sorted_indices]

def load_stack_xyz_all(file_names, chain, residues, n_cores):
    """Function that pools loading xyz for multiple cores.
    
    Args:
        file_names (list): List of PDB file paths.
        chain (str): Chain to load.
        residues (list): List of residues to load.
        n_cores (int): Number of CPU cores to use.
    
    Returns:
        torch.Tensor: Tensor with the XYZ coordinates of the CA atoms.
        torch.Tensor: Tensor with the deletion mask.
    """
    residues_dict = {i:residues.index(i) for i in residues}
    with mp.Pool(n_cores) as pool:
        xyz_list = list(tqdm(
            pool.imap(
                func=_load_ligand_ca_xyz,
                iterable=zip(file_names, repeat(chain), repeat(residues), repeat(residues_dict)),
                chunksize=max(2, round(len(file_names) / (n_cores * 16)))
            ),
            total=len(file_names),
            desc="Loading atoms"
        ))
    xyz = torch.stack(xyz_list)
    del_mask = (xyz.abs().sum(2) != 0).unsqueeze(2).float()

    return xyz, del_mask

def load_stack_xyzr_all(file_names, chain, n_cores):
    """This function loads the xyz and residue from PDBs and pads any missing residues.
    This way many PDBs can be loaded in a matrix even if they have deletions.
    
    Args:
        file_names (list): List of PDB file paths.
        chain (str): Chain to load.
        n_cores (int): Number of CPU cores to use.
    
    Returns:
        torch.Tensor: Tensor with the CA atoms and residue numbers.
        torch.Tensor: Tensor with the deletion mask.
        int: Lowest residue number.
        int: Highest residue number.
    """
    with mp.Pool(n_cores) as pool:
        xyzr_list = list(tqdm(
            pool.imap(
                func=_load_ligand_ca_xyzr,
                iterable=zip(file_names, repeat(chain)),
                chunksize=max(2, round(len(file_names) / (n_cores * 16)))
            ),
            total=len(file_names),
            desc="Loading atoms"
        ))

    low, high = float('inf'), 0
    for xyzr in xyzr_list:
        c_low = xyzr[:, -1].min().int().item()
        c_high = xyzr[: -1].max().int().item()
        low, high = min(low, c_low), max(high, c_high)

    with mp.Pool(n_cores) as pool:
        fixed_list = list(tqdm(
            pool.imap(
                func=_pad_chain,
                iterable=zip(xyzr_list, repeat(low), repeat(high)),
                chunksize=max(2, round(len(file_names) / (n_cores * 16)))
            ),
            total=len(xyzr_list),
            desc="Padding chains"
        ))

    xyz = torch.stack(fixed_list)
    del_mask = (xyz.abs().sum(2) != 0).unsqueeze(2).float()

    return xyz, del_mask, low, high


def get_interface_residues(chain1_xyzr, chain2_xyzr, chain1_del_mask, chain2_del_mask, cutoff=8.5):
    """Returns the interface residues of two CA tensors.
    
    Args:
        chain1_xyzr (torch.Tensor): Tensor with the CA atoms and residue numbers of chain 1.
        chain2_xyzr (torch.Tensor): Tensor with the CA atoms and residue numbers of chain 2.
        chain1_del_mask (torch.Tensor): Deletion mask of chain 1.
        chain2_del_mask (torch.Tensor): Deletion mask of chain 2.
        cutoff (float, optional): Cutoff for interface detection. Defaults to 8.5.
    
    Returns:
        list: List of interface residues of chain 1.
        list: List of interface residues of chain 2.
    """
    chain_1_xyz, chain_2_xyz = chain1_xyzr[:, :-1], chain2_xyzr[:, :-1]
    valid_chain_1 = chain_1_xyz[chain1_del_mask.flatten() != 0]
    valid_chain_2 = chain_2_xyz[chain2_del_mask.flatten() != 0]
    distances = torch.cdist(valid_chain_1, valid_chain_2)
    interface_pairs = distances <= cutoff
    chain_1_indices = torch.nonzero(interface_pairs.sum(1)).flatten()
    chain_2_indices = torch.nonzero(interface_pairs.sum(0)).flatten()
    chain1_residues = [int(chain1_xyzr[i, -1].item()) for i in chain_1_indices]
    chain2_residues = [int(chain2_xyzr[i, -1].item()) for i in chain_2_indices]

    return chain1_residues, chain2_residues


def _get_interface_residues(vars):
    """Wrapper function for getting the interface residues.

    Args:
        vars (tuple): Tuple with the following variables:
            chain1_xyzr (torch.Tensor): Tensor with the CA atoms and residue numbers of chain 1.
            chain2_xyzr (torch.Tensor): Tensor with the CA atoms and residue numbers of chain 2.
            chain1_del_mask (torch.Tensor): Deletion mask of chain 1.
            chain2_del_mask (torch.Tensor): Deletion mask of chain 2.
            cutoff (float): Cutoff for interface detection.
    
    Returns:
        list: List of interface residues of chain 1.
        list: List of interface residues of chain 2.
    """
    return get_interface_residues(*vars)

def calc_rmsd(models_path, rmsd_path, chain_1, chain_2, interface_cutoff=10., n_cores=mp.cpu_count(), type="interface"):
    """
    Type 'interface':
    Calculates the interface residues of two chains,
    aligns all PDBs on the union of interface residues for the first chain.
    Then calculates the pairwise RMSD with the union of interface residues for the second chain.

    Type 'ligand':
    Aligns all PDBs on the first chain.
    Then calculates the pairwise RMSD with all residues from the second chain.

    Args:
        models_path (str): Path to a directory with PDBs.
        rmsd_path (str): Output csv file path. (If None: do not write to file.)
        chain_1 (str): Receptor chain.
        chain_2 (str): Ligand chain.
        interface_cutoff (float, optional): Cutoff to use for interface detection in Angstrom. Defaults to 10.0.
        n_cores (int, optional): Number of CPU cores to use for multiprocessing. Defaults to all CPU cores.
        type (str, optional): RMSD type. Use 'ligand' or 'interface'. Defaults to "interface".

    Raises:
        Exception: "RMSD Type not 'ligand' or 'interface'."

    Returns:
        list: List of tuples with RMSDs. (file_index_1, file_index_2, rmsd)
            File indices refer to the number in the files list.
        list: List of file paths to PDBs. The RMSD tuples refer to an index in this list.
    """    
    if type not in ['ligand', 'interface']:
        raise Exception("RMSD Type not 'ligand' or 'interface'.")

    # set the number of threads to 1 to avoid issues where the amount of cores is higher than specified
    os.environ["OMP_NUM_THREADS"] = str(n_cores)
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    torch.set_num_threads(1)

    print("RMSD TYPE:", type)

    # ensure the context is spawn and avoids tne issue of using more cores than allocated
    ctx = mp.get_context("spawn")

    # Read all PDB file names from the names file.
    file_names = glob.glob(os.path.join(models_path, "*.pdb"))
    # Makes sure the files are sorted by the number in the file name
    file_names = sorted(file_names, key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[1]))
    print(f"{len(file_names)=}")

    start_time = time.perf_counter()
    # Load all residues from both chains.
    if type == 'interface':
        chain1_xyzr, chain1_del_mask, *_ = load_stack_xyzr_all(file_names, chain_1, n_cores)
        torch.cuda.empty_cache()
    chain2_xyzr, chain2_del_mask, chain2_first, chain2_last = load_stack_xyzr_all(file_names, chain_2, n_cores)
    torch.cuda.empty_cache()
    
    if type == 'interface':
        # Determine the interface residues.
        with ctx.Pool(n_cores) as pool:
            result = list(tqdm(
                pool.imap(
                    func=_get_interface_residues,
                    iterable=zip(
                        [xyzr for _, xyzr in enumerate(chain1_xyzr)],
                        [xyzr for _, xyzr in enumerate(chain2_xyzr)],
                        [del_mask for _, del_mask in enumerate(chain1_del_mask)],
                        [del_mask for _, del_mask in enumerate(chain2_del_mask)],
                        repeat(interface_cutoff)),
                    chunksize=max(2, round(len(file_names) / (n_cores * 16)))
                ),
                total=len(file_names),
                desc="Determining interface residues"
            ))
        interface_res_sets_1, interface_res_sets_2 = zip(*result)
        interface_res_1 = sorted(list(set().union(*interface_res_sets_1)))
        interface_res_2 = sorted(list(set().union(*interface_res_sets_2)))
        torch.cuda.empty_cache()
    
        print(f"{interface_res_1=}")
        print(f"{interface_res_2=}")

    with tempfile.TemporaryDirectory() as tmp_folder:
        # Align all PDBs on the interface residues of chain 1.
        if type == "ligand":
            superposition_residues = None
        else:
            superposition_residues = interface_res_1
        # runs gradpose.superpose
        gradpose.superpose(file_names, file_names[0], output=tmp_folder, residues=superposition_residues, chain=chain_1, cores=n_cores, gpu=False)
        torch.cuda.empty_cache()

        aligned_file_names = [os.path.join(tmp_folder, os.path.basename(file_name)) for file_name in file_names]

        # Here is where the types differ.
        # Ligand-RMSD: Load all residues from chain 2.
        # Interface-RMSD: Load the interface residues from chain 2.
        if type == "ligand":
            residues = list(range(chain2_first, chain2_last+1))
            xyzr_all, del_mask = load_stack_xyz_all(aligned_file_names, chain_2, residues, n_cores)
        
        if type == "interface":
            xyzr_all, del_mask = load_stack_xyz_all(aligned_file_names, chain_2, interface_res_2, n_cores)

    # Calc pairwise RMSDs of the current selection of atoms.
    with mp.get_context("spawn").Pool(n_cores) as pool:
        rmsd_lists = list(tqdm(
            pool.imap(
                func=calc_selection_rmsd,
                iterable=list(zip(reversed(range(len(xyzr_all) - 1)), repeat(xyzr_all), repeat(del_mask))),
                chunksize=max(2, round(len(file_names) / (n_cores * 16)))
            ),
            total=len(file_names) - 1,
            desc="Calculating RMSDs"
        ))
    torch.cuda.empty_cache()

    # Reverse the order of the list.
    print("Processing RMSDs...")
    rmsd_list = list(itertools.chain.from_iterable(reversed(rmsd_lists)))

    if rmsd_path:
        with open(rmsd_path, "w", encoding="utf-8") as rmsd_file:
            for rmsd in tqdm(rmsd_list, desc="Saving RMSDs to file"):
                rmsd_file.write(",".join([os.path.basename(file_names[rmsd[0]]), os.path.basename(file_names[rmsd[1]]), str(rmsd[2])]) + "\n")

    # Finished
    end_time = time.perf_counter()
    print(f"Done in {end_time - start_time:.1f}s")
    return rmsd_list, file_names