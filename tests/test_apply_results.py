"""
Name: test_apply_results.py
Function: This script is used to test the apply_results.py script. The script tests the read_rotations, read_ftresults, read_ftresults_stream, get_ftresult, apply_ftresult, and apply_ftresults_atom_group functions.
Date: 15-10-2024
Author: Nils Smit
"""

import unittest
import numpy as np
from prody import AtomGroup
import os

from scripts.apply_results import (
    read_rotations,
    read_ftresults,
    read_ftresults_stream,
    get_ftresult,
    apply_ftresult,
    apply_ftresults_atom_group,
    FTRESULT_DTYPE
)

class TestApplyResults(unittest.TestCase):

    def setUp(self):
        # Create temporary test files with sample data
        self.test_rotations_file = open("test_rotations.txt", "w")
        self.test_rotations_file.write("1 0 0 0 1 0 0 0 1\n1 0 0 0 1 0 0 0 1\n")
        self.test_rotations_file.close()

        self.test_ftresults_file = open("test_ftresults.txt", "w")
        self.test_ftresults_file.write("0 1.0 1.0 1.0 -10.0\n")
        self.test_ftresults_file.close()

    def tearDown(self):
        # Remove temporary test files
        os.remove("test_rotations.txt")
        os.remove("test_ftresults.txt")

    def test_read_rotations(self):
        rotations = read_rotations("test_rotations.txt")
        self.assertEqual(rotations.shape, (2, 3, 3))

    def test_read_ftresults(self):
        ftresults = read_ftresults(self.test_ftresults_file.name)
        
        # Check if ftresults is a structured NumPy array
        self.assertIsInstance(ftresults, np.ndarray, "Expected ftresults to be a NumPy array")
        
        # Since it's a single record, we can check for its scalar nature
        self.assertEqual(ftresults.shape, (), "Expected ftresults to have shape () for a single entry")  # Shape should be ()
        
        # Check if the structured array has the correct dtype and values
        # Check the 'roti' field
        self.assertEqual(ftresults['roti'], 0) 
        # Check the 'tv' field
        np.testing.assert_array_equal(ftresults['tv'], [1.0, 1.0, 1.0]) 
        # Check the 'E' field
        self.assertAlmostEqual(ftresults['E'], -10.0)

    def test_read_ftresults_stream(self):
        with open(self.test_ftresults_file.name, "r") as f:
            print(f"Contents of ftresults stream: {f.read()}")
            f.seek(0)  # Reset file pointer to beginning
            ftresults = read_ftresults_stream(f)
        
        # Check if ftresults from stream is a structured NumPy array
        self.assertIsInstance(ftresults, np.ndarray, "Expected ftresults from stream to be a NumPy array")
        
        # Since it's a single record, we can check for its scalar nature
        self.assertEqual(ftresults.shape, (), "Expected ftresults from stream to have shape () for a single entry")
        
        # Check if the structured array has the correct dtype and values
        # Check the 'roti' field
        self.assertEqual(ftresults['roti'], 0)
        # Check the 'tv' field
        np.testing.assert_array_equal(ftresults['tv'], [1.0, 1.0, 1.0])
        # Check the 'E' field
        self.assertAlmostEqual(ftresults['E'], -10.0)

    def test_get_ftresult(self):
        # Test the get_ftresult function
        ftresult = get_ftresult(self.test_ftresults_file.name, 0)
        self.assertEqual(ftresult['roti'], 0)

    def test_apply_ftresult(self):
        # Test the apply_ftresult function
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        ftresult = np.array((0, [1.0, 1.0, 1.0], -10.0), dtype=FTRESULT_DTYPE)
        rotations = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
        result = apply_ftresult(coords, ftresult, rotations)
        expected = np.array([[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_apply_ftresults_atom_group(self):
        # Test the apply_ftresults_atom_group function
        ag = AtomGroup("test")
        ag.setCoords(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        ftresults = np.array([(0, [1.0, 1.0, 1.0], -10.0)], dtype=FTRESULT_DTYPE)
        rotations = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
        new_ag = apply_ftresults_atom_group(ag, ftresults, rotations)
        expected_coords = np.array([[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]])
        np.testing.assert_array_almost_equal(new_ag.getCoords(), expected_coords)

if __name__ == '__main__':
    unittest.main()
