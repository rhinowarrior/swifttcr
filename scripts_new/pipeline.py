from argparse import ArgumentParser, FileType
import os.path

parser = ArgumentParser()
parser.add_argument("--receptor", "-r" , help="input pathway for the protein-MHC PDB file")
parser.add_argument("--ligand", "-l", help="input pathway for the TCR PDB file")
parser.add_argument("--output", "-o", help="output pathway for the docked complex")
parser.add_argument("--restraints", "-rs", help="input pathway for the restraints file")
parser.add_argument("--rotations", "-ro", help="input pathway for the file with rotations")

args = parser.parse_args()
receptor_dir = args.receptor_dir
ligand_dir = args.ligand_dir
output_dir = args.output_dir
restraints = args.restraints
rotations = args.rotations


print("receptor_dir: ", receptor_dir)
print("ligand_dir: ", ligand_dir)
print("output_dir: ", output_dir)
print("restraints_dir: ", restraints)
print("rotations_dir: ", rotations)




    
    