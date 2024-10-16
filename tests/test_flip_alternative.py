"""
Name: test_flip_alternative.py
Function: This script is used to test the flip_alternative.py script. The script tests the reorder_residues_in_structure function.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
from Bio import PDB
import os
from scripts.flip_alternative import reorder_residues_in_structure

class TestReorderResiduesInStructure(unittest.TestCase):

    def setUp(self):
        # Create a sample PDB file for testing
        self.test_pdb_file = 'test_input.pdb'
        self.output_pdb_file = 'test_output.pdb'
        with open(self.test_pdb_file, 'w') as f:
            f.write("""\
ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1      12.000  14.000  10.000  1.00 20.00           C  
ATOM      3  C   ALA A   1      13.000  13.000  10.000  1.00 20.00           C  
ATOM      4  O   ALA A   1      14.000  13.000  10.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1A     12.000  15.000  10.000  1.00 20.00           C  
ATOM      6  N   GLY A   2      11.104  13.207  11.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2      12.000  14.000  11.000  1.00 20.00           C  
ATOM      8  C   GLY A   2      13.000  13.000  11.000  1.00 20.00           C  
ATOM      9  O   GLY A   2      14.000  13.000  11.000  1.00 20.00           O  
ATOM     10  N   SER A   3      11.104  13.207  12.000  1.00 20.00           N  
ATOM     11  CA  SER A   3      12.000  14.000  12.000  1.00 20.00           C  
ATOM     12  C   SER A   3      13.000  13.000  12.000  1.00 20.00           C  
ATOM     13  O   SER A   3      14.000  13.000  12.000  1.00 20.00           O  
ATOM     14  CB  SER A   3A     12.000  15.000  12.000  1.00 20.00           C  
""")

    def tearDown(self):
        # Remove the test PDB files after tests
        os.remove(self.test_pdb_file)
        if os.path.exists(self.output_pdb_file):
            os.remove(self.output_pdb_file)

    def test_reorder_residues_in_structure(self):
        # Run the function to reorder residues
        reorder_residues_in_structure(self.test_pdb_file, self.output_pdb_file)

        # Parse the output PDB file
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('PDB_structure', self.output_pdb_file)

        # Check the order of residues in the output structure
        residues = list(structure.get_residues())
        self.assertEqual(residues[0].get_id()[1], 1)
        self.assertEqual(residues[0].get_id()[2], 'A')
        self.assertEqual(residues[1].get_id()[1], 1)
        self.assertEqual(residues[1].get_id()[2], ' ')
        self.assertEqual(residues[2].get_id()[1], 2)
        self.assertEqual(residues[2].get_id()[2], ' ')
        self.assertEqual(residues[3].get_id()[1], 3)
        self.assertEqual(residues[3].get_id()[2], 'A')
        self.assertEqual(residues[4].get_id()[1], 3)
        self.assertEqual(residues[4].get_id()[2], ' ')

if __name__ == '__main__':
    unittest.main()