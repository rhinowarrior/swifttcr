"""Script to merge pdb files and rename chains to A for receptor and D for ligand

Receptor (pMHC) with 3 chains A,B,C renamed to A
Ligand (TCR) with 2 chains D,E renamed to D
"""

from sys import argv
from pathlib import Path
import subprocess
import os

"""

def rename(p_rec, p_lig, p_out):
    os.chdir(p_rec.parent)
    receptor_name = "{}".format(p_rec.stem + "_rename.pdb")
    command = "pdb_selchain -A,B,C {}| pdb_chain -A > {}".format(str(p_rec), receptor_name)
    os.system(command)

"""


def merge_pdbs_main(receptor, ligand, output_dir, chains):
    p = Path(ligand)
    p_rec = Path(receptor)
    p_out = Path(output_dir)
    print("Ligand: ", p)
    print("Outputdir: ", p_out)
    print("Receptor: ", p_rec)
    os.chdir(p_rec.parent)
    receptor_name = "{}".format(p_rec.stem + "_rename.pdb")
    
    command = f"pdb_tidy {str(p_rec)} | pdb_selchain -{','.join(chains[:3])} | pdb_chain -{chains[0]} | pdb_reres -1 > {receptor_name}"
    os.system(command)
    
    if p.is_dir():
        for f in p.iterdir():
            #rename ligand
            #merge receptor and ligand
            #clean up
            if f.suffix ==".pdb":
                print(f.name)
                ligand_name = "{}".format(f.stem + "_rename.pdb")
                #command_lig =  "pdb_selchain -D,E {} | pdb_chain -D | pdb_reres -1 > {}".format(str(f), ligand_name)
                #select chain E and shiftres
                output_E = "{}".format(f.stem + "_Eshift.pdb")
                output_D = "{}".format(f.stem + "_Dshift.pdb")
                
                #add pdb_tidy before

                command_shift = f"pdb_tidy {str(f)} | pdb_selchain -{chains[4]} | pdb_shiftres -1000 | pdb_chain -{chains[3]} > {output_E}"#pdb_selchain -E 1ao7_r_u_pnon.pdb | pdb_shiftres -1000 | pdb_chain -D
                command_lig = f"pdb_tidy {str(f)} | pdb_selchain -{chains[3]} > {output_D}"
                #command_shift = "pdb_selchain -E {} | pdb_shiftres -1000 | pdb_chain -D > {}".format(str(f), output_E)#pdb_selchain -E 1ao7_r_u_pnon.pdb | pdb_shiftres -1000 | pdb_chain -D
                #command_lig = "pdb_selchain -D {} > {}".format(str(f), output_D)
                #merge chain D and E
                command_DE = "pdb_merge {} {} | pdb_sort > {}".format(output_D, output_E, ligand_name)
                #command_lig =  "pdb_sort -C {} | pdb_selchain -D,E | pdb_chain -D | pdb_reres -1 > {}".format(str(f), ligand_name)#use pdb_shift res to add 1000 to
                
                merged_name = "merged_" + str(f).split(".")[1] + ".pdb"
                merge_command = "pdb_merge {} {} | pdb_tidy -strict > {}".format(receptor_name, ligand_name, str(Path(p_out, merged_name)))
                #print(os.getcwd())
                #print("Command renaming ligand: ", command_lig)
                #print("Command merge: ", merge_command)

                os.system(command_shift)

                os.system(command_lig)

                os.system(command_DE)
                os.system(merge_command)
                #remove renamed file
                os.system("rm {}".format(output_D))
                os.system("rm {}".format(output_E))
                os.system("rm {}".format(ligand_name))

    else:
        pass



#i'll wait to remove this not sure if it is neccesary for something else

# if __name__=="__main__":
#     merge_pdbs_main()



#pdb_selchain.py -AB your.pdb | pdb_chain.py -A | pdb_reres.py -1 > new.pdb
#python pdb_merge.py <pdb file> <pdb file>
#python pdb_tidy.py -strict 1CTF.pdb

#np.subtract(orig_coords, center, orig_coords)
#out = ag.copy()
#out._setCoords(orig_coords, overwrite=True)
#writePDB("test_move.pdb", out)


#rename chains receptor

#print(command)
