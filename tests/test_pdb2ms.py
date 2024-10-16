"""
Name: test_pdb2ms.py
Function: This script is used to test the pdb2ms.py script. The script tests the pdb2ms_main, open_file, pdb2ms, add_attraction_tcr, write_ms and isPrepared functions.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
from unittest.mock import patch, mock_open, MagicMock
from pathlib import Path
from scripts.pdb2ms import pdb2ms_main, open_file, pdb2ms, add_attraction_tcr, write_ms, isPrepared

class TestPdb2MsMain(unittest.TestCase):

    @patch('scripts.pdb2ms.open_file')
    @patch('scripts.pdb2ms.pdb2ms')
    @patch('scripts.pdb2ms.add_attraction_tcr')
    @patch('scripts.pdb2ms.write_ms')
    @patch('scripts.pdb2ms.isPrepared')
    def test_pdb2ms_main(self, mock_isPrepared, mock_write_ms, mock_add_attraction_tcr, mock_pdb2ms, mock_open_file):
        # Mock the return values of the functions
        mock_isPrepared.return_value = True
        mock_open_file.return_value = ["ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N"]
        mock_pdb2ms.return_value = ["ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N   0.0"]
        mock_add_attraction_tcr.return_value = ["ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N   1.0"]

        # Define the input files
        file_1 = "/home/nils/swifttcr/scripts/test1.pdb"
        file_2 = "/home/nils/swifttcr/scripts/test2.pdb"

        # Call the function
        with patch('builtins.print') as mock_print:
            pdb2ms_main(file_1, file_2)

            # Check if the functions were called with the correct arguments
            mock_open_file.assert_any_call(Path(file_1))
            mock_open_file.assert_any_call(Path(file_2))
            mock_pdb2ms.assert_called()
            mock_add_attraction_tcr.assert_called()
            mock_write_ms.assert_called()

            # Check if the print statement was called with the correct output file names
            mock_print.assert_any_call("Wrote to: /home/nils/swifttcr/scripts/test1.ms")
            mock_print.assert_any_call("Wrote to: /home/nils/swifttcr/scripts/test2.ms")

    @patch('scripts.pdb2ms.open_file')
    @patch('scripts.pdb2ms.pdb2ms')
    @patch('scripts.pdb2ms.add_attraction_tcr')
    @patch('scripts.pdb2ms.write_ms')
    @patch('scripts.pdb2ms.isPrepared')
    def test_pdb2ms_main_skips_non_pdb(self, mock_isPrepared, mock_write_ms, mock_add_attraction_tcr, mock_pdb2ms, mock_open_file):
        # Mock the return values of the functions
        mock_isPrepared.return_value = True

        # Define the input files
        file_1 = "/home/nils/swifttcr/scripts/test1.txt"
        file_2 = "/home/nils/swifttcr/scripts/test2.pdb"

        # Call the function
        with patch('builtins.print') as mock_print:
            pdb2ms_main(file_1, file_2)

            # Check if the functions were called with the correct arguments
            mock_open_file.assert_called_once_with(Path(file_2))
            mock_pdb2ms.assert_called_once()
            mock_add_attraction_tcr.assert_called_once()
            mock_write_ms.assert_called_once()

            # Check if the print statement was called with the correct output file names
            mock_print.assert_any_call("Skipped test1")
            mock_print.assert_any_call("Wrote to: /home/nils/swifttcr/scripts/test2.ms")

if __name__ == '__main__':
    unittest.main()