import os.path
import os
import sys
import pdb2ms
import initial_placement
import subprocess
import postfilter
import apply_results
import merge_pdbs
import pairwise_rmsd
import clustering
import io
from contextlib import redirect_stdout
from argparse import ArgumentParser
import pipeline_handler
import json
import warnings

# Add the project directory to the path so we can import the modules
project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_dir)

from tools.protein_prep import prepare
from tools.ANARCI_master.Example_scripts_and_sequences import ImmunoPDB

os.environ["OMP_NUM_THREADS"] = "6"

if __name__ == "__main__":
    # Get the arguments from the user
    args = pipeline_handler.get_arguments()
    
    # Gets the reference files
    reference_ligand = os.path.realpath("ref/2bnr_r_u.pdb")
    reference_receptor = os.path.realpath("ref/2bnr_l_u.pdb")
    rotations = os.path.realpath("rotations_and_restraints/filtered_cr_in_60.prm")
    restraint_path =  os.path.realpath("rotations_and_restraints/restraintsDE.json")
    
    # checks if the output directory exists, if not it creates it
    if not os.path.exists(os.path.realpath(args.output)):
        os.mkdir(os.path.realpath(args.output))
    
    # make sure the paths are absolute and creates the variables
    receptor_path = os.path.realpath(args.pmhc)
    ligand_path = os.path.realpath(args.tcr)
    output_path = os.path.realpath(args.output) + "/"
    # restraint_path = os.path.realpath(args.restraints)
    # rotations = os.path.realpath(args.rotations)
    piper_path = os.path.realpath("tools/piper")
    #chains = args.chains
    # variable_domain = args.variabledomain
    # attractive_res = json.loads(args.attractive_res)
    
    out_ligand = "renumbered_"+ os.path.basename(ligand_path)
    out_ligand_path = ligand_path.replace(os.path.basename(ligand_path), out_ligand)

    # checks if the files exist and if the extensions are correct
    pipeline_handler.check_files(receptor_path, ligand_path, output_path, restraint_path)
    pipeline_handler.check_file_extensions(receptor_path, ligand_path, restraint_path, rotations)
    pipeline_handler.check_amount_of_chains_pdb(receptor_path, ligand_path)
    
    # Renumbers the TCR file and ignore the warnings that bio.pdb throws
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", module='Bio.PDB')
        ImmunoPDB.immunopdb_main(ligand_path, out_ligand_path)
    
    # prepares the pdb files
    receptor_pnon = prepare.prepare_main(receptor_path)
    ligand_pnon = prepare.prepare_main(out_ligand_path)
        
    # gets extra path names for later use
    receptor = os.path.basename(receptor_pnon)
    receptor_ms = os.path.basename(receptor_pnon).replace(".pdb", ".ms")

    ligand = os.path.basename(ligand_pnon)
    ligand_ms = os.path.basename(ligand_pnon).replace(".pdb", ".ms")
    
    # runs initial placement
    initial_placement.initial_placement_main(receptor_pnon, ligand_pnon, output_path, reference_receptor, reference_ligand)

    # runs pdb2ms
    pdb2ms.pdb2ms_main(output_path)

    os.chdir(output_path)

    # runs piper
    subprocess.run([
        piper_path + "/piper_attr",
        "-k1",
        "--msur_k=1.0",
        "--maskr=1.0",
        "-T", "FFTW_EXHAUSTIVE",
        "-p", piper_path + "/atoms04.prm",
        "-f", piper_path + "/coeffs04.prm",
        "-r", rotations,
        os.path.join(output_path, receptor_ms),
        os.path.join(output_path, ligand_ms)
    ]
    )

    # runs postfilter
    postfilter.post_filter_main(output_path, "ft.000.00", rotations, restraint_path, receptor, ligand, str(args.outprefix))

    # checks if the rotated directory exists, if not it creates it
    if not os.path.exists("rotated"):
        os.mkdir("rotated")
    os.chdir("rotated")

    # runs apply_results
    apply_results.apply_results_main(1000, None, None,  args.outprefix, os.path.join(output_path, "ft.000.00"), os.path.join(output_path,rotations), os.path.join(output_path,ligand))

    os.chdir("..")

    # checks if the merged directory exists, if not it creates it
    if not os.path.exists("merged"):
        os.mkdir("merged")

    # runs merge_pdbs
    merge_pdbs.merge_pdbs_main(receptor,"rotated", "merged")

    # runs pairwise_rmsd
    pairwise_rmsd.calc_rmsd("merged", "irmsd.csv", "A", "D", 10,  6)

    # runs clusteringS
    output_file = os.path.join(output_path, 'clustering.txt')

    output_buffer = io.StringIO()

    with redirect_stdout(output_buffer):
        clustering.clustering_main(os.path.join(output_path, 'irmsd.csv'))

    captured_output = output_buffer.getvalue()

    with open(output_file, 'w') as file:
        file.write(captured_output)

    print(f"Output written to {output_file}")