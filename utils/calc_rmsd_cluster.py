"""
Name: calc_rmsd_all.py
Function: Calculates the LRMSD, IRMSD and Fnat based on a reference and target structure. The RMSD is calculated by getting all the subdirectories in a given directory and getting from those subdirectories the merged directory which contains the files. The program know which reference to compare to because the reference should have the same name as the subdirectory only has the reference _renumberd added to it. The RMSD calculation is done by the tool DockQ.
Date: 25-11-2024
Author: Nils Smit, Yannick Aarts, Farzaneh Meimandi Parizi
"""

from pathlib import Path
import sys
import multiprocessing as mp
import gc
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
import time


def main():
    pipeline_dir = sys.argv[1]  # Directory where decoys are stored
    reference_dir = sys.argv[2]  # Directory where reference models are stored
    cluster_input = sys.argv[3]  # Name of the cluster file example: "clustering.txt"
    outfile_base = sys.argv[4]  # Output file path base (user-specified)
    num_threads = int(sys.argv[5])  # Number of threads (cores) to use
    batch_size = 15 # Size of each batch (can be adjusted)

    max_threads = max(1, num_threads)
    time_start = time.time()

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
            tasks.append((reference, val_path, pipeline_dir, model_identifier))

    # Process tasks in batches using Pool with map
    process_batches(tasks, max_threads, batch_size, outfile_base)
    time_end = time.time()
    print(f"DockQ calculations completed in {time_end - time_start:.2f} seconds.")


def process_batches(tasks, max_threads, batch_size, outfile_base):
    """Process tasks in batches using multiprocessing Pool with map.
    
    Args:
        tasks (list): List of tasks to process.
        max_threads (int): Maximum number of threads (cores) to use.
        batch_size (int): Size of each batch.
        outfile_base (str): Output file path base.
    """
    num_batches = (len(tasks) // batch_size) + (1 if len(tasks) % batch_size != 0 else 0)

    # Initialize dictionaries to store all results for each model_identifier
    results_dict_fnat = {}
    results_dict_lrmsd = {}
    results_dict_irmsd = {}

    for i in range(num_batches):
        batch = tasks[i * batch_size: (i + 1) * batch_size]
        
        print(f"Processing batch {i + 1}/{num_batches} with {len(batch)} tasks...")

        # Use Pool.map to process the batch sequentially in parallel workers
        with mp.Pool(processes=max_threads) as pool:
            results = pool.map(map_and_run_dockq, batch)

        # Accumulate the results for each model
        for result in results:
            model_identifier, fnat, lrmsd, irmsd = result

            # Store the results in the respective dictionaries
            if model_identifier not in results_dict_fnat:
                results_dict_fnat[model_identifier] = []
                results_dict_lrmsd[model_identifier] = []
                results_dict_irmsd[model_identifier] = []
                
            results_dict_fnat[model_identifier].append(fnat)
            results_dict_lrmsd[model_identifier].append(lrmsd)
            results_dict_irmsd[model_identifier].append(irmsd)

            # Print the result for monitoring progress
            print(f"Processed: {model_identifier} | fnat: {fnat} | LRMSD: {lrmsd} | iRMSD: {irmsd}")

        # Explicitly release memory after each batch is processed
        gc.collect()
        print(f"Batch {i + 1} completed. Memory cleaned up.")

    # Write all results after processing all batches
    write_all_results_to_files(results_dict_fnat, results_dict_lrmsd, results_dict_irmsd, outfile_base)


def write_all_results_to_files(fnat_dict, lrmsd_dict, irmsd_dict, outfile_base):
    """Write all accumulated results to separate output files in the desired format.
    
    Args:
        fnat_dict (dict): Dictionary containing FNAT results for each model.
        lrmsd_dict (dict): Dictionary containing LRMSD results for each model.
        irmsd_dict (dict): Dictionary containing IRMSD results for each model.
        outfile_base (str): Output file path base.
    """
    # Define output file paths
    fnat_outfile = Path(outfile_base + "_fnat.txt")
    lrmsd_outfile = Path(outfile_base + "_lrmsd.txt")
    irmsd_outfile = Path(outfile_base + "_irmsd.txt")

    # Write FNAT results to file
    with open(fnat_outfile, 'w') as f_fnat:
        for model_identifier, fnat_values in fnat_dict.items():
            # Join the model identifier with its corresponding fnat values
            result_line = f"{model_identifier}\t" + "\t".join(f"{round(value, 6)}" for value in fnat_values) + "\n"
            f_fnat.write(result_line)

    # Write LRMSD results to file
    with open(lrmsd_outfile, 'w') as f_lrmsd:
        for model_identifier, lrmsd_values in lrmsd_dict.items():
            # Join the model identifier with its corresponding lrmsd values
            result_line = f"{model_identifier}\t" + "\t".join(f"{round(value, 6)}" for value in lrmsd_values) + "\n"
            f_lrmsd.write(result_line)

    # Write IRMSD results to file
    with open(irmsd_outfile, 'w') as f_irmsd:
        for model_identifier, irmsd_values in irmsd_dict.items():
            # Join the model identifier with its corresponding irmsd values
            result_line = f"{model_identifier}\t" + "\t".join(f"{round(value, 6)}" for value in irmsd_values) + "\n"
            f_irmsd.write(result_line)

    print(f"Results written to {fnat_outfile}, {lrmsd_outfile}, {irmsd_outfile}")


def create_model_clus_dict(pipeline_dir, cluster_input="clustering.txt"):
    """Parse the clustering file and return a dictionary of models.
    
    Args:
        pipeline_dir (str): Path to the pipeline directory.
        cluster_input (str): Name of the clustering file.
    
    Returns:
        dict: Dictionary containing model identifiers and their respective models.
    """
    model_dict = {}
    pipeline = Path(pipeline_dir)
    for model_dir in pipeline.iterdir():
        if model_dir.is_dir():
            models = parse_clustering_file(str(Path(model_dir, cluster_input)))
            model_identifier = model_dir.name[:4]
            model_dict[model_identifier] = models
    return model_dict


def parse_clustering_file(clustering_path):
    """Parse a clustering file and return a list of model names.
    
    Args:
        clustering_path (str): Path to the clustering file.
    
    Returns:
        list: List of model names.
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


def map_and_run_dockq(args):
    """Map residues and run DockQ on the reference and decoy PDBs.
    
    Args:
        args (tuple): Tuple containing reference, decoy, pipeline_dir, model_identifier.
        
    Returns:
        tuple: Tuple containing model_identifier, fnat, lrmsd, irmsd.
    """
    reference, decoy, pipeline_dir, model_identifier = args

    # Load PDB files for decoy and reference
    target = load_PDB(decoy)
    native = load_PDB(reference)

    try:
        # Define the chain mapping
        chain_map = {"A": "A", "D": "D"}

        # Run DockQ on the reference and decoy PDBs using all native interfaces
        dockq = run_on_all_native_interfaces(target, native, chain_map=chain_map)

        # Extract the dictionary and DockQ score
        dockq_dict, dockq_score = dockq

        # Extract values of interest
        interface_key = list(dockq_dict.keys())[0]  # 'AD'
        interface_data = dockq_dict[interface_key]  # Get data for the interface

        fnat = round(interface_data['fnat'], 6)
        lrmsd = round(interface_data['LRMSD'], 6)
        irmsd = round(interface_data['iRMSD'], 6)

    finally:
        # Explicitly free memory after processing each PDB
        del target, native, dockq, dockq_dict, interface_data
        gc.collect()  # Force garbage collection after memory is freed

    # Return the result which will be handled in the callback
    return (model_identifier, fnat, lrmsd, irmsd)


if __name__ == "__main__":
    main()
