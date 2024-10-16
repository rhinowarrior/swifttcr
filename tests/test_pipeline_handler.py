"""
Name: test_pipeline_handler.py
Function: This script is used to test the pipeline_handler.py script. The script tests the get_arguments, check_files and check_file_extensions functions.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
from unittest.mock import patch, mock_open
import argparse
from scripts.pipeline_handler import get_arguments, check_files, check_file_extensions


class TestPipelineHandler(unittest.TestCase):

    @patch('argparse.ArgumentParser.parse_args')
    def test_get_arguments(self, mock_parse_args):
        # Mock the arguments from the user
        mock_parse_args.return_value = argparse.Namespace(
            pmhc='path/to/mhc.pdb',
            tcr='path/to/tcr.pdb',
            output='output/dir',
            outprefix='output_file',
            cores=4,
            threshold=9
        )
        args = get_arguments()
        # Check if the arguments are correct
        self.assertEqual(args.pmhc, 'path/to/mhc.pdb')
        self.assertEqual(args.tcr, 'path/to/tcr.pdb')
        self.assertEqual(args.output, 'output/dir')
        self.assertEqual(args.outprefix, 'output_file')
        self.assertEqual(args.cores, 4)
        self.assertEqual(args.threshold, 9)

    @patch('os.path.exists')
    def test_check_files_exist(self, mock_exists):
        # Mock the files to exist
        mock_exists.side_effect = lambda x: True
        try:
            check_files('path/to/receptor.pdb', 'path/to/ligand.pdb')
        except SystemExit:
            # Check if the function raised SystemExit unexpectedly
            self.fail("check_files() raised SystemExit unexpectedly!")

    @patch('os.path.exists')
    def test_check_files_not_exist(self, mock_exists):
        # Mock the files to not exist
        mock_exists.side_effect = lambda x: False
        with self.assertRaises(SystemExit):
            # Check if the function raises SystemExit when the files do not exist
            check_files('path/to/receptor.pdb', 'path/to/ligand.pdb')

    def test_check_file_extensions_correct(self):
        try:
            # Check if the function identifies correct file extensions
            check_file_extensions('path/to/receptor.pdb', 'path/to/ligand.pdb')
        except SystemExit:
            # Check if the function raised SystemExit unexpectedly
            self.fail("check_file_extensions() raised SystemExit unexpectedly!")

    def test_check_file_extensions_incorrect(self):
        # Check if the function raises SystemExit when the file extensions are incorrect
        with self.assertRaises(SystemExit):
            check_file_extensions('path/to/receptor.txt', 'path/to/ligand.pdb')
        with self.assertRaises(SystemExit):
            check_file_extensions('path/to/receptor.pdb', 'path/to/ligand.txt')

if __name__ == '__main__':
    unittest.main()