"""
Name: test_postfilter.py
Function: This script is used to test the postfilter.py script. The script tests the post_filter_main function.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
from unittest.mock import patch, mock_open
from scripts.postfilter import post_filter_main

class TestPostFilterMain(unittest.TestCase):
    # Mock the open function
    @patch('scripts.postfilter.parse_ft_file')
    @patch('scripts.postfilter.parse_rot_file')
    @patch('scripts.postfilter.parse_res_file')
    @patch('scripts.postfilter.post_filter')
    @patch('scripts.postfilter.filter_file_by_indices')
    @patch('builtins.open', new_callable=mock_open)
    @patch('os.path.exists')
    def test_post_filter_main(self, mock_exists, mock_open, mock_filter_file_by_indices, mock_post_filter, mock_parse_res_file, mock_parse_rot_file, mock_parse_ft_file):
        # Mock the return values of the parse functions
        mock_parse_ft_file.return_value = {1: (0, (0.0, 0.0, 0.0))}
        mock_parse_rot_file.return_value = {1: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]}
        mock_parse_res_file.return_value = {"groups": [], "required": 0}
        mock_post_filter.return_value = [1]
        mock_exists.return_value = True

        # Define the input arguments
        output_dir = "/mock/output/"
        ft_file = "mock_ft_file"
        rot_file = "mock_rot_file"
        res_file = "mock_res_file"
        receptor = "mock_receptor.pdb"
        ligand = "mock_ligand.pdb"
        outfilename = "mock_output_file"

        # Call the function under test
        post_filter_main(output_dir, ft_file, rot_file, res_file, receptor, ligand, outfilename)

        # Check if the parse functions were called with the correct arguments
        mock_parse_ft_file.assert_called_once_with(output_dir + ft_file)
        mock_parse_rot_file.assert_called_once_with(rot_file)
        mock_parse_res_file.assert_called_once_with(res_file)

        # Check if post_filter was called with the correct arguments
        mock_post_filter.assert_called_once_with(mock_parse_ft_file.return_value, mock_parse_rot_file.return_value, mock_parse_res_file.return_value, receptor, ligand)

        # Check if filter_file_by_indices was called with the correct arguments
        mock_filter_file_by_indices.assert_called_once_with(output_dir + ft_file, output_dir + outfilename, mock_post_filter.return_value)

if __name__ == '__main__':
    unittest.main()
