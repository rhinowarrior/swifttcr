######################################################################
# Pairwise RMSD calculation. Ligand RMSD and special interface RMSD. #
# Author: Kevin van Geemen                                           #
######################################################################
import os
import time
import sys
import tempfile
import itertools
import glob
from itertools import repeat
import torch
import torch.multiprocessing as mp
import gradpose
from tqdm import tqdm

# Setting these to stop issues on the cluster
mp.set_sharing_strategy('file_system')
mp.set_start_method("spawn", force=True)

def load_ligand_ca_xyzr(pdb_path, chain):
    """Loads CA atoms of a chain from a PDB.
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
    return _calc_selection_rmsd(*vars)

def _calc_selection_rmsd(template_index, xyz_all, del_mask):
    """Calculate RMSDs compared to one template.
    Only looks at PDBs after the template index to avoid duplicates.
    """
    rmsd_list = []
    xyz_all_current = xyz_all[template_index:]
    del_mask_current = del_mask[template_index:]

    pdb_range = torch.arange(start=1, end=len(xyz_all_current))
    del_mask_current = del_mask_current[:1] * del_mask_current[1:]
    # torch.save(xyz_all_current, "/home/kevinvg/superpose/xyz_all_current.pt")
    # torch.save(del_mask_current, "/home/kevinvg/superpose/del_mask_current.pt")

    rmsds = ((xyz_all_current[0] - xyz_all_current[pdb_range]) * del_mask_current).pow(2).reshape(len(xyz_all_current)-1, -1).sum(1).div(del_mask_current.sum(1).squeeze()).sqrt()
    for i in range(1, len(xyz_all_current)):
        # cur_rmsd = rmsds[i-1].item()
        # sql_rmsd = pdb2sql.StructureSimilarity.get_rmsd(xyz_all_current[0].numpy(), xyz_all_current[i].numpy())
        # print(f"COMPARE\t{cur_rmsd}\t{sql_rmsd}")
        rmsd_list.append((template_index, i + template_index, rmsds[i-1].item()))
    return rmsd_list

def _load_ligand_ca_xyzr(vars):
    return load_ligand_ca_xyzr(*vars)

def _load_ligand_ca_xyz(vars):
    return load_ligand_ca_xyz(*vars)

def _pad_chain(vars):
    return pad_chain(*vars)

def pad_chain(xyzr, low, high):
    """Fill out missing residues in a residue range with empty atoms.
    """
    for cur_res in range(low, high+1):
        if float(cur_res) not in xyzr[:, -1]:
            tmp = torch.zeros((1,4))
            tmp[:, -1] = float(cur_res)
            xyzr = torch.cat((xyzr, tmp), dim=0)
    ind = xyzr[:, -1].argsort(dim=0)
    xyzr = xyzr[ind]
    return xyzr

def load_stack_xyz_all(file_names, chain, residues, n_cores):
    """Function that pools loading xyz for multiple cores.
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
    #xyzr_list = [load_ligand_ca_xyzr(file, chain) for file in file_names]

    low = 1e10
    high = 0
    for xyzr in xyzr_list:
        c_low = xyzr[:, -1].min().int().item()
        c_high = xyzr[: -1].max().int().item()
        low = c_low if c_low < low else low
        high = c_high if c_high > high else high

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
    """
    chain_1_xyz = chain1_xyzr[:, :-1]
    chain_2_xyz = chain2_xyzr[:, :-1]

    distances = torch.cdist(chain_1_xyz, chain_2_xyz)

    inside_radius = distances <= cutoff
    inside_radius[chain1_del_mask.flatten() == 0, :] = False
    inside_radius[:, chain2_del_mask.flatten() == 0] = False
    chain_1_indices, chain_2_indices = inside_radius.nonzero(as_tuple=True)
    return [chain1_xyzr[index][-1].int().item() for _, index in enumerate(chain_1_indices)], [chain2_xyzr[index][-1].int().item() for _, index in enumerate(chain_2_indices)]


def _get_interface_residues(vars):
    return get_interface_residues(*vars)


def calc_rmsd(models_path, rmsd_path, chain_1, chain_2, interface_cutoff=10., n_cores=mp.cpu_count(), type="interface"):
    """Type 'interface':
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

    print("RMSD TYPE:", type)

    # Read all PDB file names from the names file.
    file_names = glob.glob(os.path.join(models_path, "*.pdb"))
    print(f"{len(file_names)=}")

    start_time = time.perf_counter()
    # Load all residues from both chains.
    if type == 'interface':
        chain1_xyzr, chain1_del_mask, *_ = load_stack_xyzr_all(file_names, chain_1, n_cores)
    chain2_xyzr, chain2_del_mask, chain2_first, chain2_last = load_stack_xyzr_all(file_names, chain_2, n_cores)

    if type == 'interface':
        # Determine the interface residues.
        with mp.Pool(n_cores) as pool:
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
    
        print(f"{interface_res_1=}")
        print(f"{interface_res_2=}")

    with tempfile.TemporaryDirectory() as tmp_folder:
        # Align all PDBs on the interface residues of chain 1.
        if type == "ligand":
            superposition_residues = None
        else:
            superposition_residues = interface_res_1
        gradpose.superpose(file_names, file_names[0], output=tmp_folder, residues=superposition_residues, chain=chain_1, cores=n_cores, gpu=torch.cuda.is_available())

        aligned_file_names = [os.path.join(tmp_folder, os.path.basename(file_name)) for file_name in file_names]

        # Here is where the types differ.
        # Ligand-RMSD: Load all residues from chain 2.
        # Interface-RMSD: Load the interface residues from chain 2.
        if type == "ligand":
            residues = list(range(chain2_first, chain2_last+1))
            xyzr_all, del_mask = load_stack_xyz_all(aligned_file_names, chain_2, residues, n_cores)
        
        if type == "interface":
            xyzr_all, del_mask = load_stack_xyz_all(aligned_file_names, chain_2, interface_res_2, n_cores)
            # interface_res_1 = list(zip(repeat(chain_1), interface_res_1))
            # interface_res_2 = list(zip(repeat(chain_2), interface_res_2))

            # chain_res_pairs = interface_res_1 + interface_res_2
            # xyzr_all, del_mask = load_union_xyz(aligned_file_names, chain_res_pairs)


    # print("After loading data", flush=True)
    # print(xyzr_all.shape, flush=True)
    # print(del_mask.shape, flush=True)
    # print("===")

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

    # Reverse the order of the list.
    print("Processing RMSDs...")
    rmsd_list = list(itertools.chain.from_iterable(reversed(rmsd_lists)))

    # Write RMSDs to file.
    if rmsd_path:
        with open(rmsd_path, "w", encoding="utf-8") as rmsd_file:
            for rmsd in tqdm(rmsd_list, desc="Saving RMSDs to file"):
                # rmsd_file.write(",".join([str(num) for num in rmsd]) + "\n")
                rmsd_file.write(",".join([os.path.basename(file_names[rmsd[0]]), os.path.basename(file_names[rmsd[1]]), str(rmsd[2])]) + "\n")

    # Finished
    end_time = time.perf_counter()
    print(f"Done in {end_time - start_time:.1f}s")
    return rmsd_list, file_names


def pairwise_rmsd_main(models_path, rmsd_file, chain_1, chain_2, rmsd_type, n_cores):
    """Arguments: path/to/pdbs output_file.csv ChainID1 ChainID2 RMSDtype num_cpu_cores
    RMSDtype can be ligand or interface.
    num_cpu_cores is optional and will default to all cores.
    Interface cutoff cannot be set like this at the moment.
    (This whole section should be rewritten with argparse,
    or import this script and call calc_rmsd() with the arguments.)
    """
    
    # Now uses calc_rmsd() directly.

    print(f"{models_path=}")
    print(f"{rmsd_file=}")
    print(f"chains: {chain_1}, {chain_2}")
    print(f"{rmsd_type=}")
    print(f"{n_cores=}")
    calc_rmsd(models_path, rmsd_file, chain_1, chain_2, interface_cutoff=10., n_cores=n_cores, type=rmsd_type)

# if __name__ == "__main__":
#     pairwise_rmsd_main()