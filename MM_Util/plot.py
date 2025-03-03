import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
import os
import sys
from .utilities import *
from .mol import *

class Plot():
    """
    A class that creates plots from a numpy array.
    It can make reaction coordinate diagrams, RMSD and RMSF trajectory plots,
    Free energy plots, Scatter Plots, and SNFG Figures.
    """

    def __init__(self, data=None):
        """
        Constructs a plot object.
        :param data: (array) a numpy array containing the data to be plotted.
        """

        # Constructor Attributes
        self.data = data
        self.fig = None
        self.ax = None
        self.path = None

    def trajectory(self, var_name = 'colvar'):
        """ Plots MD trajectory with histogram. Takes in data for CP2K or Gromacs via Mol.

        :param var_name: (list) Name of the collective variable you are plotting on your y-axis.
        """
        
        time = self.data[:, 0]
        colvar = self.data[:, 1]
        
        timestep = np.abs(time[0] - time[1])

        fig, ax = plt.subplots(1,2, figsize=(11,3), gridspec_kw={'width_ratios': [3.5, 1]})

        ax[0].plot(time, colvar, linewidth=0.2)
        ax[0].set_xlabel(f"time (fs); stepsize = {timestep}fs")
        ax[0].set_ylabel(var_name)
        # ax1.set_title(f"file: {xyz_file}", fontsize = 10)

        midpt = int(np.round(len(colvar) / 2))


        ax[1].hist(colvar[0:midpt], bins='rice', fc=(0, 0, 1, 0.3), orientation="horizontal") # First half shown in blue
        ax[1].hist(colvar[midpt:-1], bins='rice', fc=(0, 0, 1, 0.5), orientation="horizontal") # Second half shown in red
        # ax2.axhline(y=np.average(colvar), color='b', linewidth=2)

        ax[1].set_title(f"average = {np.round(np.average(colvar), 3)}", fontsize = 10)

        ax[1].set_xlabel('structures')

        plt.tight_layout()
        
        self.fig = fig, self.ax = ax


    def reaction_profile(self, reaction,type='delta E', linewidth=3, scale=0.32, annotate=True, color='black'):
        """
        Plots a reaction coordinate diagram.
        """
        mol_list= reaction.mol_list
        labels = reaction.mol_label

        energies = []
        # labels = []

        for mol in mol_list:
            energy = mol.E  # Access the 'E' attribute which stores the energy
            energies.append(energy)
        
        # for mol in mol_label:
        #     label = mol.label
        #     labels.append(label)

        fig, ax = plt.subplots(figsize=(8, 6))
        
        relative_energies = [hartree_to_kcal(e - energies[0]) for e in energies]
        print (relative_energies)

        annotation_offset=0.05
        
        for j, energy in enumerate(relative_energies):
            # Draw Horizontal Bars at Each Energy Level
            ax.plot([(j + 1 - scale), (j + 1 + scale)], [energy, energy], 
                    color=color, linewidth=linewidth)   

            # Annotate Energy Values
            if annotate and j != 0:
                ax.text(j + 1,  energy + annotation_offset, f"{energy:.2f}", fontsize=12, ha='center', color='black')

            # Draw Dashed Connecting Lines
            if j < len(relative_energies) - 1:
                ax.plot([(j + 1 + scale), (j + 2 - scale)], 
                        [energy, relative_energies[j + 1]], 
                        linestyle=":", color=color, linewidth=linewidth)



        # Set Y-axis Label
        if type == 'delta E':
            reaction_type = '$\\Delta E$ (kcal $\\cdot$ mol${}^{-1}$)' 
        elif type == 'delta F':
            reaction_type = '$\\Delta F$ (kcal $\\cdot$'

        # Invisible plot for the legend label
        ax.plot([], [], color=color, linewidth=linewidth, label=reaction_type)
        
        # Customize X-axis Ticks
        ax.set_xticks(range(1, len(energies) + 1))
        # ax.set_xticklabels(labels[i] for i in mol_label)
        ax.set_xticklabels(i for i in labels)
        # Add Legend
        ax.legend(loc="lower left", frameon=False, fontsize=14)

        self.fig = fig; self.ax = ax
        # return fig, ax