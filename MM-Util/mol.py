import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from Cython.Shadow import nonecheck


class Mol():

    """
    A class that can parse log files from different Computational Software.
    It can parse log files from Gaussian 16, xvg files from Gromacs, out files from ORCA, and log files from CP2K.
    """

    def __init__(self, path=None):

        """
        Constructs a molecule object.

        :param path: (path) Path to the directory containing the log file(s).
        """

        # Constructor Attributes
        self.path = path
        self.xyz = None #Nx3 array containg the xyz coords for all atoms
        self.conn_mat = None
        self.atoms = None

        self.energy = None
        self.f_Energy = None
        self.enthalpy = None

        self.coupling = None
        self.rmsd = []
        self.rmsf = []



class reaction():

    """
    A class that organizes several molecules into a reaction
    """

    def __init__(self, mol_list, mol_label):

        """
        Constructs reaction from a list of molecules.

        :param mol_list: (list) List of molecule objects.
        :param mol_label: (list) List of reaction labels to be used (Reactant, Intermediate, Product, Transition State)
        """

        self.energies = []
        self.enthalpies = []
        self.f_energies = []

        self.mol_list = mol_list
        self.mol_label = mol_label

    @staticmethod
    def combiner(mol_list):

        """
        Adds the energies for multiple molecules and returns a new mol object.

        :param mol_list: (list) List of molecule objects.
        :return: mol
        """
        new_mol = Mol()

        for mol in mol_list:

            new_mol.energy += mol.energy
            new_mol.enthalpy += mol.enthalpy
            new_mol.f_Energy += mol.f_Energy

        return new_mol