"""
    Psi extension

    Print the phi backbone angle for each residue in the structure.
    Psi angle is determined by the coordinates of the C(i-1), N(i), CA(i), and
    C(i). atoms.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "11 February 2013"
__author__ = "David Hall"

def usage():
    str  = "Print a PDB file named consistently with libmol pdb format"
    return str

def run_extension(routines, outroot, options):
    """
        Print the list of phi angles

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".mol.pdb"
    file = open(outname, "w")

    # Initialize some variables

    protein = routines.protein

    for atom in protein.getAtoms():

        file.write("%s\n" % atom.getPDBString())

    file.close()
