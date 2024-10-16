"""
Name: test_merge_pdbs.py
Function: This script is used to test the merge_pdbs.py script. The script tests the run_command function.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
import sys
import os
from unittest.mock import patch, MagicMock
import subprocess

# Add the scripts directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'scripts')))

from scripts.merge_pdbs import run_command

class TestRunCommand(unittest.TestCase):

    @patch('subprocess.run')
    def test_run_command_success(self, mock_run):
        # Mock the subprocess.run to return a successful result
        mock_result = MagicMock()
        mock_result.stdout = "Command executed successfully"
        mock_run.return_value = mock_result

        command = "echo 'Hello, World!'"
        result = run_command(command)
        
        # Assert that the command was executed successfully
        mock_run.assert_called_once_with(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        self.assertEqual(result, "Command executed successfully")

    @patch('subprocess.run')
    def test_run_command_failure(self, mock_run):
        # Mock the subprocess.run to raise a CalledProcessError
        mock_run.side_effect = subprocess.CalledProcessError(returncode=1, cmd="invalid_command", stderr="Error message")

        command = "invalid_command"
        
        with self.assertRaises(subprocess.CalledProcessError):
            run_command(command)
        
        # Assert that the command was attempted to be executed
        mock_run.assert_called_once_with(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

if __name__ == '__main__':
    unittest.main()
