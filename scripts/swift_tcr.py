import os.path
import os
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


os.environ["OMP_NUM_THREADS"] = "6"

if __name__ == "__main__":
    # Get the arguments from the user
    args = pipeline_handler.get_arguments()
    
    #Do we want to change these ever or keep them the same?
    reference_ligand = "/home/nils/swifttcr/example/input/2bnr_r_u.pdb"
    reference_receptor = "/home/nils/swifttcr/example/input/2bnr_l_u.pdb"
    rotations = "/home/nils/swifttcr/example/input/filtered_cr_in_60.prm"
    
    #make sure the paths are absolute
    receptor_path = os.path.realpath(args.receptor)
    ligand_path = os.path.realpath(args.ligand)
    output_path = os.path.realpath(args.output) + "/"
    restraint_path = os.path.realpath(args.restraints)

    receptor = os.path.basename(receptor_path)
    receptor_ms = os.path.basename(receptor_path).replace(".pdb", ".ms")

    ligand = os.path.basename(ligand_path)
    ligand_ms = os.path.basename(ligand_path).replace(".pdb", ".ms")
    
    #checks if the files exist and if the extensions are correct
    pipeline_handler.check_files(receptor_path, ligand_path, output_path, restraint_path)
    pipeline_handler.check_file_extensions(receptor_path, ligand_path, restraint_path, rotations)
    
    pipeline_handler.check_chains_pdb(receptor_path, ligand_path)
    
    # #runs initial placement
    # initial_placement.initial_placement_main(receptor_path, ligand_path, output_path, reference_receptor, reference_ligand)

    # # runs pdb2ms
    # pdb2ms.pdb2ms_main(output_path)

    # os.chdir(output_path)

    # #runs piper
    # subprocess.run([
    #     "/home/nils/swifttcr/piper/piper_attr",
    #     "-k1",
    #     "--msur_k=1.0",
    #     "--maskr=1.0",
    #     "-T", "FFTW_EXHAUSTIVE",
    #     "-p", "/home/nils/swifttcr/piper/atoms04.prm",
    #     "-f", "/home/nils/swifttcr/piper/coeffs04.prm",
    #     "-r", rotations,
    #     os.path.join(output_path, receptor_ms),
    #     os.path.join(output_path, ligand_ms)
    # ]
    # )

    # #runs postfilter
    # postfilter.post_filter_main(output_path, "ft.000.00", rotations, restraint_path, receptor, ligand, str(args.outprefix))

    # #checks if the rotated directory exists, if not it creates it
    # if not os.path.exists("rotated"):
    #     os.mkdir("rotated")
    # os.chdir("rotated")

    # #runs apply_results
    # apply_results.apply_results_main(100, None, None,  args.outprefix, os.path.join(output_path, "ft.000.00"), os.path.join(output_path,rotations), os.path.join(output_path,ligand))

    # os.chdir("..")

    # #checks if the merged directory exists, if not it creates it
    # if not os.path.exists("merged"):
    #     os.mkdir("merged")

    # #runs merge_pdbs
    # merge_pdbs.merge_pdbs_main(receptor,"rotated", "merged")

    # #runs pairwise_rmsd
    # pairwise_rmsd.calc_rmsd("merged", "irmsd.csv", "A", "D", 10,  6)

    # #runs clusteringS
    # output_file = os.path.join(output_path, 'clustering.txt')

    # output_buffer = io.StringIO()

    # with redirect_stdout(output_buffer):
    #     clustering.clustering_main(os.path.join(output_path, 'irmsd.csv'))

    # captured_output = output_buffer.getvalue()

    # with open(output_file, 'w') as file:
    #     file.write(captured_output)

    # print(f"Output written to {output_file}")