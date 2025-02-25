import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re
import networkx as nx
from networkx.utils.misc import flatten

from utilities import *


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

    def element_symbol(A):

        """ A dictionary for atomic number and atomic symbol
        :param A: either atomic number or atomic symbol for Hydrogen, Carbon, Nitrogen, Oxygen, Fluorine and Silicon
        :return: the corresponding atomic symbol or atomic number
        """

        periodic_table = {'1': 'H', '2': 'He',
                          '3': 'Li', '4': 'Be', '5': 'B', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '10': 'Ne',
                          '11': 'Na', '12': 'Mg', '13': 'Al', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl', '18': 'Ar',
                          '19': 'K', '20': 'Ca', '35': 'Br'
                          }

    def gaussian(self, path=None):

        """
        Creates a mol object from a Gaussian 16 log file.
        :param path: Path to the directory containing the log file.
        :return: Mol Object
        """

        flags = {'freq_flag': False, 'nmr_flag': False, 'opt_flag': False, 'jcoup_flag': False, 'normal_mode': False,
                 'read_geom': False}
        job_type = None

        # temprorary variables to hold the data
        freq = [];
        ints = [];
        vibs = [];
        geom = [];
        atoms = [];
        nmr = []
        self.NAtoms = None

        for line in open("/".join([self.path, "input.log"]), 'r').readlines():

            if job_type is None and line is not "\n":
                if "opt" in line:
                    if "freq" in line:
                        job_type = 'optfreq'
                    else:
                        job_type = 'opt'
                elif "freq" in line:
                    if "opt" in line:
                        job_type = 'optfreq'
                    else:
                        job_type = 'freq'
                        flags["freq_flag"] = True
                elif "nmr" in line:
                    job_type = 'nmr'
                else:
                    job_type = 'sp'

            if self.NAtoms is None and re.search('^ NAtoms=', line):
                self.NAtoms = int(line.split()[1])

            if job_type == 'optfreq' or job_type == "freq":

                if flags['freq_flag'] == False and re.search('Normal termination', line): flags['freq_flag'] = True
                # We skip the opt part of optfreq job, all info is in the freq part
                if flags['freq_flag'] == True:

                    if re.search('SCF Done', line):
                        self.E = float(line.split()[3])
                    elif re.search('Sum of electronic and zero-point Energies', line):
                        self.Ezpe = float(line.split()[-1])
                    elif re.search('Sum of electronic and thermal Enthalpies', line):
                        self.H = float(line.split()[-1])
                    elif re.search('Sum of electronic and thermal Free Energies', line):
                        self.F = float(line.split()[-1])

                    elif re.search('Coordinates', line) and len(geom) == 0:
                        flags['read_geom'] = True

                    elif flags['read_geom'] == True and re.search('^\s*.\d', line):
                        geom.append([float(x) for x in line.split()[3:6]])
                        atoms.append(Mol.element_symbol(line.split()[1]))
                        if int(line.split()[0]) == self.NAtoms:
                            flags['read_geom'] = False


                    elif re.search('Deg. of freedom', line):
                        self.NVibs = int(line.split()[-1])

                    elif re.search('^ Frequencies', line):
                        freq_line = line.strip()
                        for f in freq_line.split()[2:5]: freq.append(float(f))
                        flags['normal_mode'] = False

                    elif re.search('^ IR Inten', line):
                        ir_line = line.strip()
                        for i in ir_line.split()[3:6]: ints.append(float(i))

                    elif re.search('^  Atom  AN', line):
                        flags['normal_mode'] = True  # locating normal modes of a frequency
                        mode_1 = [];
                        mode_2 = [];
                        mode_3 = [];
                        # continue

                    elif flags['normal_mode'] == True and re.search('^\s*\d*\s*.\d*', line) and len(line.split()) > 3:
                        mode_1.append([float(x) for x in line.split()[2:5]])
                        mode_2.append([float(x) for x in line.split()[5:8]])
                        mode_3.append([float(x) for x in line.split()[8:11]])

                    elif flags['normal_mode'] == True:
                        flags['normal_mode'] = False
                        for m in [mode_1, mode_2, mode_3]: vibs.append(np.array(m))

            elif job_type == 'opt':

                if re.search('SCF Done', line): E = float(line.split()[4])
                if re.search('Optimization completed.', line):
                    self.E = E;
                    flags['opt_flag'] = True
                if flags['opt_flag'] == True:
                    if re.search('Standard orientation:', line):
                        flags['read_geom'] = True

                    elif flags['read_geom'] == True and re.search('^\s*.\d', line):
                        geom.append([float(x) for x in line.split()[3:6]])
                        atoms.append(Mol.element_symbol(line.split()[1]))
                        if int(line.split()[0]) == self.NAtoms:
                            flags['read_geom'] = False

            elif job_type == 'nmr':

                if re.search('SCF Done', line):
                    self.E = float(line.split()[4])
                elif re.search('Coordinates', line) and len(geom) == 0:
                    flags['read_geom'] = True

                elif flags['read_geom'] == True and re.search('^\s*.\d', line):
                    geom.append([float(x) for x in line.split()[3:6]])
                    atoms.append(Mol.element_symbol(line.split()[1]))
                    if int(line.split()[0]) == self.NAtoms:
                        flags['read_geom'] = False

                elif re.search('Total nuclear spin-spin coupling J', line):
                    spin = [[] for i in range(self.NAtoms)]
                    flags['jcoup_flag'] = True

                elif flags['jcoup_flag'] == True and re.search('-?\d\.\d+[Dd][+\-]\d\d?', line):
                    for x in line.split()[1:]:
                        spin[int(line.split()[0]) - 1].append(float(x.replace('D', 'E')))

                elif flags['jcoup_flag'] == True and re.search('End of Minotr F.D. properties file', line):
                    flags['jcoup_flag'] = False

            elif job_type == 'sp':

                if re.search('SCF Done', line):
                    self.E = float(line.split()[4])
                elif re.search('Standard orientation:', line):
                    flags['read_geom'] = True
                elif flags['read_geom'] == True and re.search('^\s*.\d', line):
                    geom.append([float(x) for x in line.split()[3:6]])
                    atoms.append(Mol.element_symbol(line.split()[1]))
                    if int(line.split()[0]) == self.NAtoms:
                        flags['read_geom'] = False

        # postprocessing:
        if job_type == 'freq' or job_type == 'optfreq':
            self.Freq = np.array(freq)
            self.Ints = np.array(ints)
            self.Vibs = np.zeros((self.NVibs, self.NAtoms, 3))
            for i in range(self.NVibs): self.Vibs[i, :, :] = vibs[i]

        if job_type == 'nmr':
            for at in spin:
                while len(at) < self.NAtoms: at.append(0)
            self.nmr = np.tril(spin)

        self.xyz = np.array(geom)
        self.atoms = atoms

    def connectivity_matrix(self, distXX=1.65, distXH=1.15):

        """ Creates a connectivity matrix of the molecule. A connectivity matrix holds the information of which atoms are bonded and to what.

        :param distXX: The max distance between two atoms (not hydrogen) to be considered a bond
        :param distXH: The max distance between any atom and a hydrogen atom to be considered a bond
        """

        Nat = self.NAtoms
        self.conn_mat = np.zeros((Nat, Nat))

        for at1 in range(Nat):
            for at2 in range(Nat):

                dist = get_distance(self.xyz[at1], self.xyz[at2])

                if at1 == at2: pass

                elif (self.atoms[at1] == 'H' or self.atoms[at2] == 'H') and dist < distXH:
                    self.conn_mat[at1,at2] = 1; self.conn_mat[at2,at1] = 1
                elif (self.atoms[at1] != 'H' and self.atoms[at2] != 'H') and dist < distXX:
                    self.conn_mat[at1,at2] = 1; self.conn_mat[at2,at1] = 1

        #Remove bifurcated Hs:
        for at1 in range(Nat):
            if self.atoms[at1] == 'H' and np.sum(self.conn_mat[at1,:]) > 1:

                    at2list = np.where(self.conn_mat[at1,:] == 1)
                    at2list = at2list[0].tolist()

                    at2dist = [ round(get_distance(self.xyz[at1], self.xyz[at2x]), 3) for at2x in at2list]
                    for at,dist in zip(at2list, at2dist):
                        if self.atoms[at] == 'H':
                            at2list.remove(at)
                            at2dist.remove(dist)

                    at2 = at2list[at2dist.index(min(at2dist))]
                    for at2x in at2list:
                        if at2x != at2:
                            print('remove', self._id, at2x, at1, at2)
                            self.conn_mat[at1, at2x] = 0 ; self.conn_mat[at2x, at1] = 0

        graph = nx.graph.Graph(self.conn_mat)
        self.Nmols = nx.number_connected_components(graph)
        self.graph = graph

class Reaction():

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
class Reaction():

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