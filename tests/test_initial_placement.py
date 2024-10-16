"""
Name: test_initial_placement.py
Function: This script is used to test the initial_placement.py script. The script tests the superpose_and_change_chain_IDs function and the initial_placement_main function.
Date: 14-10-2024
Author: Nils Smit
"""

import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path

from scripts.initial_placement import superpose_and_change_chain_IDs, initial_placement_main

class TestInitialPlacement(unittest.TestCase):

    @patch('scripts.initial_placement.pymol_cmd')
    def test_superpose_and_change_chain_IDs(self, mock_pymol_cmd):
        # Mock the PyMOL commands
        # Mock get_raw_alignment to return a single alignment
        mock_pymol_cmd.get_raw_alignment.return_value = [
        [('ref', 1), ('target', 2)]
        ]
        
        # Mock get_model to return appropriate atoms for each chain (ref and target)
        mock_pymol_cmd.get_model.side_effect = [
            MagicMock(atom=[MagicMock(chain='A', resi='1'), MagicMock(chain='A', resi='2')]),  # ref atoms
            MagicMock(atom=[MagicMock(chain='B', resi='1'), MagicMock(chain='B', resi='2')]),  # target atoms
            MagicMock(atom=[MagicMock(resi='1'), MagicMock(resi='2')]),  # ref chain A
            MagicMock(atom=[MagicMock(resi='1'), MagicMock(resi='2')]),  # target chain B
        ]
        
        # Mock get_chains to return chains for reference and target
        mock_pymol_cmd.get_chains.side_effect = [
            ['A'],  # Reference chains
            ['B']   # Target chains
        ]
        
        # Call the function under test
        superpose_and_change_chain_IDs('reference.pdb', 'target.pdb', 'output.pdb')
            
        # Assertions to check the sequence of PyMOL commands
        
        # Check if PyMOL was reinitialized and the reference and target structures were loaded
        mock_pymol_cmd.reinitialize.assert_called_once()
        mock_pymol_cmd.load.assert_any_call('reference.pdb', 'ref')
        mock_pymol_cmd.load.assert_any_call('target.pdb', 'target')
        
        # Check if superposition was performed correctly
        mock_pymol_cmd.super.assert_called_once_with("target and name CA", "ref and name CA", object="alignment")
        # Check if get_raw_alignment was called with the correct object
        mock_pymol_cmd.get_raw_alignment.assert_called_once_with("alignment")
        
        # Check if chain mapping logic was triggered correctly
        mock_pymol_cmd.get_model.assert_any_call("ref")
        mock_pymol_cmd.get_model.assert_any_call("target")
        mock_pymol_cmd.get_chains.assert_any_call("ref")
        mock_pymol_cmd.get_chains.assert_any_call("target")
        # Check if the peptide chain was identified correctly
        mock_pymol_cmd.create.assert_any_call("target_chain_B", "target and chain B")
        # Check if the chain was renamed correctly
        mock_pymol_cmd.alter.assert_any_call("target_chain_B", "chain='A'")
        mock_pymol_cmd.create.assert_called_with("renamed_target", "target_chain_*")
        # Check if the structure was saved correctly
        mock_pymol_cmd.save.assert_called_once_with('output.pdb', "renamed_target")
        # Check if the workspace was cleared
        mock_pymol_cmd.remove.assert_called_once_with("all")

    @patch('scripts.initial_placement.superpose_and_change_chain_IDs')
    def test_initial_placement_main(self, mock_superpose_and_change_chain_IDs):
        # Prepare test paths
        receptor = Path('receptor.pdb')
        ligand = Path('ligand.pdb')
        outputdir = Path('outputdir')
        reference_receptor = Path('reference_receptor.pdb')
        reference_ligand = Path('reference_ligand.pdb')
        
        # Call the function under test
        initial_placement_main(receptor, ligand, outputdir, reference_receptor, reference_ligand)
        
        # Verify that superpose_and_change_chain_IDs was called twice (once for receptor and once for ligand)
        self.assertEqual(mock_superpose_and_change_chain_IDs.call_count, 2)
        mock_superpose_and_change_chain_IDs.assert_any_call(
            Path('reference_receptor.pdb'), Path('receptor.pdb'), Path('outputdir/receptor.pdb')
        )
        mock_superpose_and_change_chain_IDs.assert_any_call(
            Path('reference_ligand.pdb'), Path('ligand.pdb'), Path('outputdir/ligand.pdb')
        )



if __name__ == '__main__':
    unittest.main()