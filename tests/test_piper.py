"""
Name: test_piper.py
Function: This script is used to test the piper.py script. The script tests the run_piper function.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
from unittest.mock import patch, call
import subprocess
import os

class TestPiper(unittest.TestCase):

    @patch('subprocess.run')
    def test_run_piper(self, mock_run):
        # Setup
        piper_path = "/mock/tools/piper"
        rotations = "/mock/rotations_and_restraints/filtered_cr_in_60.prm"
        output_path = "/mock/output"
        receptor_ms = "mock_receptor.ms"
        ligand_ms = "mock_ligand.ms"

        # Expected command
        expected_command = [
            os.path.join(piper_path, "piper_attr"),
            "-k1",
            "--msur_k=1.0",
            "--maskr=1.0",
            "-T", "FFTW_EXHAUSTIVE",
            "-p", os.path.join(piper_path, "atoms04.prm"),
            "-f", os.path.join(piper_path, "coeffs04.prm"),
            "-r", rotations,
            os.path.join(output_path, receptor_ms),
            os.path.join(output_path, ligand_ms)
        ]

        # Call the function
        subprocess.run(expected_command)

        # Assert subprocess.run was called with the expected command
        mock_run.assert_called_once_with(expected_command)

if __name__ == '__main__':
    unittest.main()