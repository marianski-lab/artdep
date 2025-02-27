from unittest import TestCase
from unittest.mock import patch, mock_open
import numpy as np
from mol import *


class TestMol(TestCase):

    @patch("builtins.open", new_callable=mock_open, read_data="""
            # opt freq
            Normal termination
            SCF Done:  E(RB3LYP)= -100.12345 A.U. after 10 cycles
            8
            NAtoms= 3
            Coordinates (Angstroms)
            Standard orientation:
            1   8   0   0.000   0.000   0.000
            2   1   0   0.757   0.586   0.000
            3   1   0   -0.757   0.586   0.000
            Sum of electronic and zero-point Energies= -100.23456
            Sum of electronic and thermal Enthalpies= -100.34567
            Sum of electronic and thermal Free Energies= -100.4567
            Frequencies --  1650.0  3200.0  3500.0
            IR Inten    --  10.5    25.2    30.1
            Deg. of freedom 4
            Normal termination
        """)
    def test_gaussian(self, mock_file):

        # Create a test Mol object
        mol = Mol(path="test_dir")  # Ensure path is set

        # Call the gaussian function
        mol.gaussian()

        # Check extracted energy values
        self.assertAlmostEqual(mol.E, -100.12345)
        self.assertAlmostEqual(mol.Ezpe, -100.23456)
        self.assertAlmostEqual(mol.H, -100.34567)
        self.assertAlmostEqual(mol.F, -100.45678)

        # Check number of atoms
        self.assertEqual(mol.NAtoms, 3)

        # Check extracted geometry
        expected_geom = np.array([
            [0.000, 0.000, 0.000],
            [0.757, 0.586, 0.000],
            [-0.757, 0.586, 0.000]
        ])
        np.testing.assert_array_almost_equal(mol.xyz, expected_geom)

        # Check extracted frequencies
        expected_freq = np.array([1650.0, 3200.0, 3500.0])
        np.testing.assert_array_almost_equal(mol.Freq, expected_freq)

        # Check IR intensities
        expected_ints = np.array([10.5, 25.2, 30.1])
        np.testing.assert_array_almost_equal(mol.Ints, expected_ints)

        # Ensure that `open()` was called with the correct file path
        mock_file.assert_called_with("test_dir/input.log", 'r')

    def test_connectivity_matrix(self):
        # Create a test molecule (water H2O)
        mol = Mol()
        mol.NAtoms = 3
        mol.xyz = np.array([
            [0.000, 0.000, 0.000],  # Oxygen
            [0.757, 0.586, 0.000],  # Hydrogen 1
            [-0.757, 0.586, 0.000]  # Hydrogen 2
        ])
        mol.atoms = ['O', 'H', 'H']

        # Call the connectivity matrix function
        mol.connectivity_matrix()

        # Expected connectivity matrix for H₂O:
        expected_conn_mat = np.array([
            [0, 1, 1],  # Oxygen is connected to both hydrogens
            [1, 0, 0],  # Hydrogen 1 is connected to oxygen
            [1, 0, 0]   # Hydrogen 2 is connected to oxygen
        ])

        # Check if the connectivity matrix matches the expected matrix
        np.testing.assert_array_equal(mol.conn_mat, expected_conn_mat)

        # Check if the molecule is correctly identified as a single connected component
        self.assertEqual(mol.Nmols, 1)

        # Ensure that the generated graph has the correct number of nodes and edges
        self.assertEqual(mol.graph.number_of_nodes(), 3)
        self.assertEqual(mol.graph.number_of_edges(), 2)