"""
Name: calculate_rmsd_all.py
Function: Calculate LRMSD, IRMSD, and FNAT for all models in the mapped directory. Data is read from first file to the last of each model. So, the files wil be read from 0 to n for each model.
Date: 17-10-2024
Author: Jan Aarts, Farzaneh Meimandi Parizi, Nils Smit
"""

"""
Input:
    - mapped_dir: Directory containing the mapped files
    - outfile: Base output file name
    - num_threads: Number of threads (cores) to use for multiprocessing

Example: python calculate_rmsd_all.py mapped_dir outfile num_threads
"""

from Bio.PDB import PDBParser
from pathlib import Path
import os
from sys import argv
from pdb2sql.StructureSimilarity import StructureSimilarity
from multiprocessing import Pool
import gc

def add_suffix(path, suffix):
    """Adds a suffix to the filename and returns the new Path object.
    
    Args:
        path (Path): The original Path object.
        suffix (str): The suffix to add to the filename.
    
    Returns:
        Path: The new Path object with the added suffix.
    """
    before = path.parent
    after = path.stem
    return Path(before, f"{after}_{suffix}{path.suffix}")

def calc_LRMSD_decoys(ref_file, decoy_files, thread_num):
    """Calculates LRMSD, IRMSD, and FNAT between reference and decoy files.
    
    Args:
        ref_file (str): Path to the reference file.
        decoy_files (list): List of paths to the decoy files.
        thread_num (int): The thread number for creating a unique Lzone file.
        
    Returns:
        tuple: Tuple containing lists of LRMSD, IRMSD, and FNAT values.
    """
    lrmsds, irmsds, fnats = [], [], []

    # Create a thread-specific Lzone file name
    lzone_file = f'ref_lzone_{thread_num}.lzone'
    
    # Delete any existing Lzone file before processing
    if os.path.exists(lzone_file):
        try:
            os.remove(lzone_file)
        except Exception as e:
            print(f"Error deleting pre-existing Lzone file {lzone_file}: {str(e)}")

    try:
        for decoy_file in decoy_files:
            try:
                # Create a StructureSimilarity object
                sim = StructureSimilarity(decoy_file, ref_file, enforce_residue_matching=False)

                # Calculate LRMSD, IRMSD, and FNAT values
                lrmsd_pdb2sql = sim.compute_lrmsd_fast(lzone=lzone_file)
                irmsd_pdb2sql = sim.compute_irmsd_fast()
                fnat_pdb2sql = sim.compute_fnat_fast()

                lrmsds.append(lrmsd_pdb2sql)
                irmsds.append(irmsd_pdb2sql)
                fnats.append(fnat_pdb2sql)

                print(f"Thread {thread_num}: {os.path.basename(decoy_file)} - LRMSD: {lrmsd_pdb2sql}, IRMSD: {irmsd_pdb2sql}, FNAT: {fnat_pdb2sql}")

            except Exception as e:
                print(f"Error in {os.path.basename(decoy_file)}: {str(e)}")
            
            finally:
                # Clean up the StructureSimilarity object
                del sim
                gc.collect()

    finally:
        # Ensure the Lzone file is deleted after processing
        if os.path.exists(lzone_file):
            try:
                os.remove(lzone_file)
            except Exception as e:
                print(f"Error deleting Lzone file {lzone_file}: {str(e)}")

    return lrmsds, irmsds, fnats

def process_model(index, model_identifier, file_number, mapped_dir, thread_num):
    """Processes a single model and calculates RMSD values.
    
    Args:
        index: Index of the model
        model_identifier: Model identifier
        file_number: File number
        mapped_dir: Path to the mapped directory
        thread_num: Thread number for unique Lzone file
    
    Returns:
        Tuple with the results of the RMSD calculations.
    """
    # Use only the first four characters for the model identifier
    short_model_identifier = model_identifier[:4]
    
    # Create paths to the reference and decoy files
    ref_file = mapped_dir / f"{short_model_identifier}_reference_{file_number}_mapped.pdb"
    decoy_file = mapped_dir / f"{short_model_identifier}_merged_{file_number}_mapped.pdb"

    if ref_file.exists() and decoy_file.exists():
        print(f"Thread {thread_num}: Calculating RMSD for {ref_file.name} and {decoy_file.name}")
        lrmsd_values, irmsd_values, fnat_values = calc_LRMSD_decoys(ref_file, [decoy_file], thread_num)
        return index, short_model_identifier, file_number, lrmsd_values, irmsd_values, fnat_values
    else:
        print(f"Files do not exist: {ref_file} or {decoy_file}")
        return index, short_model_identifier, file_number, [], [], []  # Empty results for missing files

def set_thread_limit(num_threads):
    """
    Sets thread limit for multi-threaded libraries like NumPy, MKL, OpenMP.
    
    Args:
        num_threads (int): Number of threads to use.
    """
    os.environ["OMP_NUM_THREADS"] = str(num_threads)
    os.environ["MKL_NUM_THREADS"] = str(num_threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(num_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(num_threads)

def main():
    mapped_dir = Path(argv[1])        # Directory containing the mapped files
    outfile = Path(argv[2])           # Base output file name
    num_threads = int(argv[3])        # Number of threads (cores) specified by the user

    # Prepare a list of model identifiers and their corresponding file numbers
    model_dict = {}
    
    set_thread_limit(1)  # Set thread limit for multi-threaded libraries
    
    for file in mapped_dir.glob("*_mapped.pdb"):
        parts = file.stem.split("_")
        # Check if the file name has at least 4 parts
        if len(parts) >= 4:
            # Extract the model ID
            model_identifier = "_".join(parts[:-2])
            # Extract the file number from the second-to-last part
            file_number = parts[-2]
            short_model_identifier = model_identifier[:4]
            if short_model_identifier not in model_dict:
                model_dict[short_model_identifier] = []
            if file_number not in model_dict[short_model_identifier]:
                model_dict[short_model_identifier].append(file_number)

    results = {}

    # Sort the model_dict by model identifier and then by file number
    for model_identifier in model_dict:
        model_dict[model_identifier].sort(key=int)  # Sort numerically

    with Pool(num_threads) as pool:
        futures = []
        index = 0
        
        # Submit jobs to the pool and sort the results
        for model_identifier, file_numbers in sorted(model_dict.items()):  
            for file_number in file_numbers:
                thread_num = (index % num_threads) + 1
                futures.append(pool.apply_async(process_model, (index, model_identifier, file_number, mapped_dir, thread_num)))
                index += 1
        
        # As each thread completes, gather results
        for future in futures:
            result = future.get()  # Blocks until the result is available
            model_id = result[1]  # Extract model identifier
            if model_id not in results:
                results[model_id] = [[], [], []]  # Initialize lists for LRMSD, IRMSD, FNAT
            # Extend the lists with the computed values
            results[model_id][0].extend(result[3])  # LRMSD values
            results[model_id][1].extend(result[4])  # IRMSD values
            results[model_id][2].extend(result[5])  # FNAT values

    # Open files for writing the separate results
    lrmsd_outfile = add_suffix(outfile, "lrmsd")
    irmsd_outfile = add_suffix(outfile, "irmsd")
    fnat_outfile = add_suffix(outfile, "fnat")
    combined_outfile = add_suffix(outfile, "combined")  # Combined output file

    with open(lrmsd_outfile, 'w') as f_lrmsd, \
         open(irmsd_outfile, 'w') as f_irmsd, \
         open(fnat_outfile, 'w') as f_fnat, \
         open(combined_outfile, 'w') as f_combined:

        for model_identifier, (lrmsd_values, irmsd_values, fnat_values) in results.items():
            # Write separate files for LRMSD, IRMSD, FNAT
            f_lrmsd.write(f"{model_identifier}\t" + "\t".join(f"{val}" for val in lrmsd_values) + "\n")
            f_irmsd.write(f"{model_identifier}\t" + "\t".join(f"{val}" for val in irmsd_values) + "\n")
            f_fnat.write(f"{model_identifier}\t" + "\t".join(f"{val}" for val in fnat_values) + "\n")

            # Combine results into a formatted output for the combined file
            combined_line = f"{model_identifier}\t" + "\t".join(f"({lrmsd}, {irmsd}, {fnat})" for lrmsd, irmsd, fnat in zip(lrmsd_values, irmsd_values, fnat_values))
            combined_line += "\n"  # Add a newline at the end
            f_combined.write(combined_line)

    print("Results written to", lrmsd_outfile, irmsd_outfile, fnat_outfile, combined_outfile)


if __name__ == "__main__":
    main()
