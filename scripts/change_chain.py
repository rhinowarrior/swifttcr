# import sys
# import os
# import subprocess

# project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# sys.path.append(project_dir)
# sys.path.append('tools/PDB_matching')

# from tools.PDB_matching import pdb_match_chn, pdb_rename_chain

# def get_chain_ids(pdb_file):
#     with open(pdb_file, 'r') as f:
#         lines = f.readlines()
#     chains = []
#     for line in lines:
#         if line.startswith('ATOM'):
#             if line[21] not in chains:
#                 chains.append(line[21])
#     return chains

# def get_shortest_chain(chains, pdb_file):
#     shortest_chain = None
#     shortest_chain_length = None
#     for chain in chains:
#         chain_length = 0
#         with open(pdb_file, 'r') as f:
#             lines = f.readlines()
#         for line in lines:
#             if line.startswith('ATOM'):
#                 if line[21] == chain:
#                     chain_length += 1
#         if shortest_chain_length == None or chain_length < shortest_chain_length:
#             shortest_chain = chain
#             shortest_chain_length = chain_length
#     return shortest_chain

# def compare_to_ref(ref_pdb, other_pdb, peptide = None):
#     # Define the output directory and file
#     output_dir = "example/input/match"
#     output_file = os.path.join(output_dir, "match.txt")
    
#     # Ensure the directory exists
#     os.makedirs(output_dir, exist_ok=True)
    
#     # Open the file and write the output
#     with open(output_file, "w") as f:
#         # Save the current sys.stdout (so you can restore it later)
#         original_stdout = sys.stdout
        
#         # Redirect sys.stdout to the file
#         sys.stdout = f
        
#         try:
#             # Call the function that normally prints to stdout
#             pdb_match_chn.pdb_match_chn_main(ref_pdb, other_pdb)
#         finally:
#             # Restore the original stdout so print statements go back to the console
#             sys.stdout = original_stdout
    
#     # Now, let's check the content of the file and replace "None" with the peptide variable
#     with open(output_file, "r") as f:
#         content = f.read()  # Read the entire file content
    
#     # Replace "None" with the peptide value
#     if peptide:
#         updated_content = content.replace("None", peptide)
#     else:
#         updated_content = content
    
#     # Write the updated content back to the file
#     with open(output_file, "w") as f:
#         f.write(updated_content)
        
# def compare_to_ref_tcr(ref_pdb, other_pdb):
#     # Define the output directory and file
#     output_dir = "example/input/match"
#     output_file = os.path.join(output_dir, "match.txt")
    
#     # Ensure the directory exists
#     os.makedirs(output_dir, exist_ok=True)
    
#     # Open the file and write the output
#     with open(output_file, "w") as f:
#         # Save the current sys.stdout (so you can restore it later)
#         original_stdout = sys.stdout
        
#         # Redirect sys.stdout to the file
#         sys.stdout = f
        
#         try:
#             # Call the function that normally prints to stdout
#             pdb_match_chn.pdb_match_chn_main(ref_pdb, other_pdb)
#         finally:
#             # Restore the original stdout so print statements go back to the console
#             sys.stdout = original_stdout
    
# def change_chain(pdb_file):
#     # Define the output directory and file
#     outdir = "example/input/changed_chain"
#     os.makedirs(outdir, exist_ok=True)
#     output_file = os.path.join(outdir, os.path.basename(pdb_file).replace(".pdb", "_new.pdb"))
#     match_file = "example/input/match/match.txt"
#     pdb_rename_chain.rename_chain(pdb_file, match_file, output_file)
#     return output_file

# def are_chain_ids_equal(pdb_file1, pdb_file2):
#     # Get chain IDs for both files
#     chains1 = get_chain_ids(pdb_file1)
#     chains2 = get_chain_ids(pdb_file2)
    
#     # Convert lists to sets for easy comparison
#     set_chains1 = set(chains1)
#     set_chains2 = set(chains2)
    
#     # Check if both sets are equal
#     return set_chains1 == set_chains2

# def change_chain_main(ref_tcr, ref_pmhc, tcr, pmhc):
#     # Get the arguments from the user
#     chains = get_chain_ids(pmhc)
#     peptide = get_shortest_chain(chains, pmhc)
#     if not are_chain_ids_equal(ref_pmhc, pmhc):
#         compare_to_ref(ref_pmhc, pmhc, peptide)
#         pmhc_path = change_chain(pmhc)
#     else:
#         pmhc_path = pmhc
    
#     if not are_chain_ids_equal(ref_tcr, tcr):
#         compare_to_ref_tcr(ref_tcr, tcr)
#         tcr_path = change_chain(tcr)
#     else:
#         tcr_path = tcr
#     return pmhc_path, tcr_path
    
# if __name__ == '__main__':
#     ref_pmhc = "/home/nils/swifttcr/ref/2bnr_l_u.pdb"
#     ref_tcr = "/home/nils/swifttcr/ref/2bnr_r_u.pdb"
#     pmhc = "/home/nils/swifttcr/example/input/5c0c/5c0c_pmhc_pnon.pdb"
#     tcr = "/home/nils/swifttcr/example/input/5c0c/5c0c_tcr_pnon.pdb"
#     pmhc_path, tcr_path = change_chain_main(ref_tcr, ref_pmhc, tcr, pmhc)
#     print(pmhc_path, tcr_path)