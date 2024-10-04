"""Script to parse clustering data for all cases.

Example usage: python3 python_programs/get_results_pipeline.py pipeline/ rmsds_benchmark_7.txt ExpandedBenchmark
/imgt/ clustering_7.txt

"""

from sys import argv
from pathlib import Path
from Bio.PDB import Select
from Bio.PDB import PDBParser, PDBIO, Selection
from Bio import pairwise2
from Bio import BiopythonWarning
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner

import os
import sys
import subprocess
import glob
import warnings

import pdb2sql
from pdb2sql.StructureSimilarity import StructureSimilarity
warnings.simplefilter('ignore', BiopythonWarning)
PATH_extension = ""

def main():
    
    pipeline_dir = argv[1]#pipeline/
    outfile = argv[2]
    reference_dir = argv[3]#path to directory with reference models, name should match case id
    cluster_input = argv[4]#name of cluster file example: "clustering_8.txt"
    model_dict = create_model_clus_dict(pipeline_dir, cluster_input)
    p = Path(outfile)
    f = open(str(add_suffix(p, "lrmsd_")), 'w')
    g = open(str(add_suffix(p, "irmsd_")), 'w')
    h = open(str(add_suffix(p, "fnat_")), 'w')
    
    for key, value in model_dict.items():
        reference = str(Path(reference_dir, key + "_merged.pdb"))
        print(reference)
        val_paths = []
        for val in value:
             val_path = str(Path(pipeline_dir, key + PATH_extension, "merged", val))
             val_paths.append(val_path)


        print(reference, val_paths)
        lrmsds, irmsds, fnats = calc_LRMSD_decoys(reference, val_paths)
        outline = key + "".join(["\t" + str(lrmsd) for lrmsd in lrmsds]) + "\n"
        i_outline = key + "".join(["\t" + str(irmsd) for irmsd in irmsds]) + "\n"
        f_outline = key + "".join(["\t" + str(fnat) for fnat in fnats]) + "\n"
        f.write(outline)
        g.write(i_outline)
        h.write(f_outline)
    print("Wrote results to ", outfile)
    f.close()
    g.close()
    h.close()
    
def add_suffix(p, suffix):
    before = p.parent
    after = p.stem
    return Path(before, suffix + after)    


def parse_clustering_file(clustering_path):
    models = []
    f = open(clustering_path)
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.startswith("Cluster center: "):
            model_name = line.split("Cluster center: ")[-1]
            model_name = model_name.split("with")[0].strip()
            models.append(model_name)
    return models

def create_model_clus_dict(p_dir, cluster_input = "clustering.txt"):
    model_dict = {}
    p = Path(p_dir)
    for model_dir in p.iterdir():
        if PATH_extension in model_dir.name:
            models = parse_clustering_file(str(Path(model_dir, cluster_input)))
            key = model_dir.name[:4]#-3
            model_dict[key] = models
    return model_dict


def calc_LRMSD_decoys(ref_file, decoys):
    """
    Calculates the LRMSD between a reference PDB structure and a set of decoy PDB structures.
    Args:
        in_DIR: str, input directory containing the PDB structures.
        target_id: str, ID of the target to calculate LRMSD for.
        ref_name: str, name of the reference PDB structure (default='ref.pdb').
        reg_expr: str, regular expression to match the decoy PDB structures (default='*.pdb').
    Returns:
        None.
    """
    final_scores = {}

    pdb_ref = PDBParser().get_structure('pdb1', ref_file)    


    # iterating over all files
    #print(decoy_files)
    lrmsds = []
    irmsds = []
    fnats = []
    for decoy_file in decoys:
        pdb_decoy = PDBParser().get_structure('pdb2', decoy_file)

        
        # # remove c-like domain and keep only g domain
        # pdb_decoy = remove_C_like_domain(pdb_decoy)
        # pdb_ref = remove_C_like_domain(pdb_ref)


        pdb_ref , pdb_decoy= map_PDBs(pdb_ref,pdb_decoy)

        pdb_ref = reres(pdb_ref)
        pdb_decoy = reres(pdb_decoy)

        write_PDB(pdb_ref, ref_file[:-4]+'_mapped.pdb')
        write_PDB(pdb_decoy, decoy_file[:-4]+'_mapped.pdb')
        lrmsd_pdb2sql = 1000
        irmsd_pdb2sql = 1000
        fnat_pdb2sql = 1000
        # print(f"Reference file: {ref_file[:-4]+'_mapped.pdb'}")
        # print(f"Decoy file: {decoy_file[:-4]+'_mapped.pdb'}")
        
        try:
            sim = StructureSimilarity(ref_file[:-4]+'_mapped.pdb',decoy_file[:-4]+'_mapped.pdb', enforce_residue_matching=False)

            lrmsd_pdb2sql = sim.compute_lrmsd_fast(lzone='ref.lzone')
            irmsd_pdb2sql = sim.compute_irmsd_fast()
            fnat_pdb2sql = sim.compute_fnat_fast()
            print(os.path.basename(decoy_file), lrmsd_pdb2sql)

        except Exception as e:
            print("Error in")
            print(f"Error in {os.path.basename(decoy_file)}: {str(e)}")
            
        try:
            os.remove(decoy_file[:-4]+'_mapped.pdb')
            os.remove(ref_file[:-4]+'_mapped.pdb')
        except:
            print('error deleting mapped files')
        lrmsds.append(lrmsd_pdb2sql)
        irmsds.append(irmsd_pdb2sql)
        fnats.append(fnat_pdb2sql)

    return lrmsds, irmsds, fnats

def write_PDBs(pdb_ref, pdb_decoy, output_dir, target_id):
    decoy_path = '%s/%s_decoy_mapped.pdb' % (output_dir, target_id)

    io = PDBIO()
    io.set_structure(pdb_decoy)
    io.save(decoy_path, select=NotDisordered())

    ref_path = '%s/%s_ref_mapped.pdb' % (output_dir, target_id)

    io = PDBIO()
    io.set_structure(pdb_ref)
    io.save(ref_path, select=NotDisordered())
    return True

def check_atom_mapping( pdb1, pdb2):
    """
    Checks if the atoms in the input PDB structures are correctly mapped.
    
    Parameters:
    pdb1 (Bio.PDB): The reference PDB structure
    pdb2 (Bio.PDB): The PDB structure to be compared with the reference
    """

    # Get lists of chains in each structure
    chains1 = pdb1.get_chains()
    chains2 = pdb2.get_chains()

    # Iterate over chains and residues in each structure
    for chain1, chain2 in zip(chains1, chains2):
       residues1 = chain1.get_residues()
       residues2 = chain2.get_residues()
   
       for residue1, residue2 in zip(residues1, residues2):
           atoms1 = residue1.get_atoms()
           atoms2 = residue2.get_atoms()

           atom_names1= []
           atom_names2= []
           # Check if corresponding atoms have the same name
           for atom1 in atoms1:
               atom_names1.append(atom1.id)
           for atom2 in atoms2:
               atom_names2.append(atom2.id)

           mismatched_atoms = list(set(atom_names2).symmetric_difference(set(atom_names1)))


           backbone_mismatches =  any(elem in mismatched_atoms for elem in ['C', 'CA', 'N', 'O'])
           if backbone_mismatches:
               print('Backbone mismatch',backbone_mismatches)
        
           if (mismatched_atoms):
                print(mismatched_atoms)

    
def remove_selected_residues( pdb, residues_to_remove):
    """
    Removes selected residues from a PDB structure.
    
    Parameters:
    pdb (Bio.PDB): The PDB structure to modify
    residues_to_remove (dict): A dictionary with chain IDs as keys and lists of residue IDs as values
    
    Returns:
    Bio.PDB: The modified PDB structure
    """
    chains_to_remove = [k for k in residues_to_remove]

    for chain in pdb.get_chains():
        if chain.id in chains_to_remove:
            for residue in chain.get_residues():
                if residue.id[1] in residues_to_remove.get(chain.id):
                    #print(residue.id ,' -> ',residue.id[1])
                    residue.id = ('X', residue.id[1], residue.id[2])
                    #print(residue.id ,' <- ')
            
            residues= reversed(list(enumerate(chain)))
            for indx, res in residues: 
                if 'X' == res.id[0]:
                    #print('Changed: ',res.id, indx)
                    chain.detach_child(('X', res.id[1], res.id[2]))
    return pdb

def remove_to_last_residue(pdb, chain_id, begin_range):
        chain = pdb[0][chain_id]
        # get the last residue object in the chain
        last_residue = chain.get_list()[-1]
        last_residue_number = last_residue.get_id()[1]

        residue_ids_to_remove = {chain_id: range(begin_range , last_residue_number)}
        if residue_ids_to_remove:
            pdb = remove_selected_residues( pdb, residue_ids_to_remove)
        return pdb


def assign(pdb, pdb_sequences ,  residues_to_remove=[]):
        ''' Renumbers the pdb using aligned sequences with residue numberings in pdb_sequences. 
        Args:
            pdb_ref:   Bio.PDB object
            pdb_decoy: Bio.PDB object
            
        Returns: Bio.PDB objects with renumbered residues
        '''
        
        for chain in pdb.get_chains():

            seq = pdb_sequences.get(chain.id)
            #print(seq)      
            for ind, res in enumerate(chain):
                        #print(res.id , ' ->',('X', seq[ind], res.id[2]))
                        res.id = ('X', seq[ind], res.id[2])

        for chain in pdb.get_chains():
            for res in chain:
                res.id = (' ', res.id[1], ' ')
        
        if residues_to_remove:
            pdb = remove_selected_residues( pdb, residues_to_remove)

        return pdb


def map_PDBs(pdb_ref,pdb_decoy,remove_non_mapped=True,custom_map={}):
    
    ref_sequences = [[chain.id, seq1(''.join([res.resname for res in chain]), custom_map=custom_map)] for chain in pdb_ref.get_chains()]

    ref_sequences.sort()
    decoy_sequences = [[chain.id, seq1(''.join([res.resname for res in chain]), custom_map=custom_map)] for chain in pdb_decoy.get_chains()]
                            
    decoy_sequences.sort()
    
    assert(len(ref_sequences) == len(decoy_sequences))
    
    for ind in range(len(ref_sequences)):
            pair = pairwise2.align.globalxx(ref_sequences[ind][1], decoy_sequences[ind][1])[0]
            ref_sequences[ind][1]   = pair.seqA
            decoy_sequences[ind][1] = pair.seqB


    ref_sequences = [[seq[0],[i+1 for i,res in enumerate(seq[1]) if res != '-']] for seq in ref_sequences]
    decoy_sequences = [[seq[0],[i+1 for i,res in enumerate(seq[1]) if res != '-']] for seq in decoy_sequences]

    ref_gapped_residues = [[ref_sequences[j][0],list(set(ref_sequences[j][1]) -set(decoy_sequences[j][1])) ] for j in range(len(ref_sequences))]  
    decoy_gapped_residues = [[decoy_sequences[j][0],list(set(decoy_sequences[j][1]) -set(ref_sequences[j][1])) ] for j in range(len(decoy_sequences))]  


    ref_sequences = {key:val for key, val in ref_sequences}
    decoy_sequences = {key:val for key, val in decoy_sequences}
    
    ref_gapped_residues = {key:val for key, val in ref_gapped_residues}
    decoy_gapped_residues = {key:val for key, val in decoy_gapped_residues}
 
    if not all(x==[] for x in ref_gapped_residues.values()):
            print('After mapping: Residues removed from reference PDB:')
            print(ref_gapped_residues)
    if not all(x==[] for x in decoy_gapped_residues.values()):
            print('After mapping: Residues removed from decoy     PDB:')
            print(decoy_gapped_residues)

    if remove_non_mapped:
        pdb_ref = assign(pdb_ref, ref_sequences,ref_gapped_residues )
        pdb_decoy = assign(pdb_decoy, decoy_sequences,  decoy_gapped_residues)
    else:
        pdb_ref = assign(pdb_ref, ref_sequences )
        pdb_decoy = assign(pdb_decoy, decoy_sequences)
    return pdb_ref, pdb_decoy

def write_PDB(pdb,  file_name):
    io = PDBIO()
    io.set_structure(pdb)
    io.save(file_name, select=NotDisordered())
    


def remove_C_like_domain(pdb):
    """ Removes the C-like domain from a MHC struture and keeps only the G domain
    Args:
        pdb: (Bio.PDB): Bio.PDB object with chains names M (N for MHCII) and P
        need_to_be_removed (list, optional):list of atoms to remove from M chain. Defaults to None.
    Returns: (Bio.PDB): Bio.PDB object without the C-like domain
    """
    chain_ids = [chain.id for chain in pdb.get_chains()]

    # If MHCII, remove the C-like domain from the M-chain (res 80 and higher) and the N-chain (res 90 and higher)
    if 'N' in chain_ids:
        pdb = remove_to_last_residue(pdb, 'N' , 90)
        if 'M' in chain_ids:
            pdb = remove_to_last_residue(pdb, 'M' , 80)


    # If MHCI, remove the C-like domain, which is from residue 180+
    if 'N' not in chain_ids and 'M' in chain_ids:

        pdb = remove_to_last_residue(pdb, 'M' , 180)
        if 'B'  in chain_ids:
            pdb = pdb[0].detach_child('B')
            #pdb = remove_to_last_residue(pdb, 'B' , 1)
                
    return pdb


class NotDisordered(Select):  # Inherit methods from Select class
    '''
    Keep one Alternative location for the given atom
    '''
    def accept_atom(self, atom):
        keepAltID = 'A'
        if (not atom.is_disordered()) or atom.get_altloc() == keepAltID:
            atom.set_altloc(" ")  # Eliminate alt location ID before output.
            return True
        else:  # Alt location was not one to be output.
            return False

def reres(pdb):
    """
    Renumbers the residues in a PDB structure.
    """

    new_residueNum = [[chain.id,[i+1 for i,res in enumerate(chain.get_residues())]] for chain in pdb.get_chains()]    
    new_residueNum = {key:val for key, val in new_residueNum}
    assign(pdb, new_residueNum)
    return pdb

if __name__ == "__main__":
    main()