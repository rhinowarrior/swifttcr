"""
Name: test_pairwise_rmsd.py
Function: This script is used to test the pairwise_rmsd.py script. The script tests the calc_rmsd function.
Date: 15-10-2024
Author: Nils Smit
"""

import os
import unittest
import tempfile
import torch
import multiprocessing as mp
from scripts.pairwise_rmsd import calc_rmsd

class TestCalcRMSD(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory to store PDB files
        self.test_dir = tempfile.TemporaryDirectory()
        self.models_path = self.test_dir.name
        self.rmsd_path = os.path.join(self.test_dir.name, "rmsd.csv")

        # Create a dummy PDB file with multiple chains
        pdb_content = (
            "ATOM      1  N   GLY A   1     -16.725  11.718  16.390  1.00 27.64           N\n"
            "ATOM      2  CA  GLY A   1     -15.242  11.785  16.856  1.00 26.32           C\n"
            "ATOM      3  C   GLY A   1     -14.247  11.660  15.713  1.00 26.13           C\n"
            "ATOM      4  O   GLY A   1     -14.656  11.634  14.538  1.00 26.11           O\n"
            "ATOM      5  H   GLY A   1     -16.911  12.551  15.892  1.00  0.00           H\n"
            "ATOM      6  N   SER A   2     -12.956  11.614  16.045  1.00 24.29           N\n"
            "ATOM      7  CA  SER A   2     -11.919  11.396  15.054  1.00 23.11           C\n"
            "ATOM      8  C   SER A   2     -11.901   9.896  14.614  1.00 21.27           C\n"
            "ATOM      1  N   LYS D   2     -13.586  -2.007 -15.767  1.00 66.58           N\n"
            "ATOM      2  CA  LYS D   2     -13.064  -0.618 -15.822  1.00 64.19           C\n"
            "ATOM      3  C   LYS D   2     -12.903  -0.053 -17.241  1.00 61.50           C\n"
            "ATOM      4  O   LYS D   2     -12.592  -0.763 -18.211  1.00 61.69           O\n"
            "ATOM      5  CB  LYS D   2     -11.732  -0.489 -15.094  1.00 64.01           C\n"
            "ATOM      6  CG  LYS D   2     -11.856  -0.107 -13.636  1.00 67.60           C\n"
            "ATOM      7  CD  LYS D   2     -10.526  -0.334 -12.924  1.00 70.19           C\n"
            "ATOM      8  CE  LYS D   2     -10.772  -0.897 -11.537  1.00 73.45           C\n"
            "ATOM      9  NZ  LYS D   2      -9.621  -1.750 -11.125  1.00 75.26           N1+\n"
            "ATOM     10  H   LYS D   2     -14.340  -2.098 -16.416  1.00  0.00           H\n"
            "ATOM     11  HZ1 LYS D   2      -9.426  -2.413 -11.847  1.00  0.00           H\n"
            "ATOM     12  HZ2 LYS D   2      -9.852  -2.233 -10.281  1.00  0.00           H\n"
            "ATOM     13  HZ3 LYS D   2      -8.819  -1.173 -10.973  1.00  0.00           H\n"
            "ATOM     14  N   GLU D   3     -13.109   1.247 -17.327  1.00 57.83           N\n"
            "ATOM     15  CA  GLU D   3     -13.204   1.924 -18.580  1.00 55.06           C\n"
            "ATOM     16  C   GLU D   3     -12.148   3.045 -18.527  1.00 50.19           C\n"
            "END\n"
        )

        # Write the PDB content to a file
        with open(os.path.join(self.models_path, "model1.pdb"), 'w') as f:
            f.write(pdb_content)

        self.chain_1 = "A"
        self.chain_2 = "D"

    # Only test the Ligand RMSD calculation, because we don't use the interface RMSD calculation in the actual code
    def test_calc_rmsd_ligand(self):
        rmsd_list, file_names = calc_rmsd(
            models_path=self.models_path,
            rmsd_path=self.rmsd_path,
            chain_1=self.chain_1,
            chain_2=self.chain_2,
            interface_cutoff=10.0,
            n_cores=1,
            type="ligand"
        )
        self.assertIsInstance(rmsd_list, list)
        self.assertIsInstance(file_names, list)
        self.assertTrue(os.path.exists(self.rmsd_path))

    # Test the Interface RMSD calculation
    def test_invalid_rmsd_type(self):
        with self.assertRaises(Exception) as context:
            calc_rmsd(
                models_path=self.models_path,
                rmsd_path=self.rmsd_path,
                chain_1=self.chain_1,
                chain_2=self.chain_2,
                interface_cutoff=10.0,
                n_cores=1,
                type="invalid_type"
            )
        self.assertTrue("RMSD Type not 'ligand' or 'interface'." in str(context.exception))

if __name__ == '__main__':
    unittest.main()