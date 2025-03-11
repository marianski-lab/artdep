import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
import os
import sys

from contourpy.util import data
from matplotlib import colors

from .utilities import *
from .mol import *

class Plot():
    """
    A class that creates plots from a numpy array.
    It can make reaction coordinate diagrams, RMSD and RMSF trajectory plots,
    Free energy plots, Scatter Plots, and SNFG Figures.
    """

    def __init__(self, data, labels:list = None, desc:list = None, xtick:int = None, ytick:int = None, xrange:list = None, yrange:list = None, colors:list = None):
        """
        Constructs a plot object.
        :param data: (array) a numpy array containing the data to be plotted.
        """

        # Constructor Attributes
        self.data = data
        self.labels = labels
        self.desc = desc
        self.xtick = xtick
        self.ytick = ytick
        self.xrange = xrange
        self.yrange = yrange
        self.colors = colors

    def trajectory(self, molecule, var_name = 'colvar'):
        """ Plots MD trajectory with histogram. Takes in data for CP2K or Gromacs via Mol.
        :param molecule: (Mol) Class Mol. 
        :param var_name: (list) Name of the collective variable you are plotting on your y-axis.
        """
        
        time = molecule.data[:, 0]
        colvar = molecule.data[:, 1]
        
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


    def reaction_profile(self, Reaction, type, linewidth=3, scale=0.32, annotate=True, color='black'):
        """
        Plots a reaction coordinate diagram.
        """
        mol_list = Reaction.mol_list
        labels = Reaction.mol_label
        
        energies = []

        '''
        type = 'delta E' or 'delta F' or 'delta H'
        '''

        for mol in mol_list:
            if type == 'delta E':
                energy = mol.E  
            elif type == 'delta F':
                energy = mol.F
            elif type == 'delta H':
                energy = mol.H
            else:
                print("Unsupported Energy Type")
                return  
            
            energies.append(energy)  
       
        if not energies:
            raise ValueError("No energies found. Check the input data.")

        # Dynamic Figure Size Based on Number of Reaction Steps
        num_steps = len(mol_list)
        fig_width = max(6, num_steps * 2)  # Adjust width based on number of steps
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        

        relative_energies = [hartree_to_kcal(e - energies[0]) for e in energies]

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
            reaction_type = '$\\Delta F$ (kcal $\\cdot$ mol${}^{-1}$)'
        elif type == 'delta H':
            reaction_type = '$\\Delta H$ (kcal $\\cdot$ mol${}^{-1}$)' 

        # Invisible plot for the legend label
        ax.plot([], [], color=color, linewidth=linewidth)

        # Add X-axis Guide Line at the halfway point
        ax.axhline(0, color="black", linestyle=":", linewidth=1.5, zorder=-4)
        
        ax.set_ylabel(f'{reaction_type}', fontsize=16)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: proper_minus(x)))

        # Customize X-axis Ticks
        ax.set_xticks(range(1, len(energies) + 1))
        ax.set_xticklabels(labels) 
            
        # Make bottom and left borders thicker and visible
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_linewidth(1.5)  # Thicker bottom border
        ax.spines['left'].set_linewidth(1.5)    # Thicker left border

        # Final Formatting
        ax.tick_params(labelsize=14)
        ax.legend(loc="lower left", frameon=False, fontsize=14)

        plt.tight_layout()

        self.fig = fig; self.ax = ax
        # return fig, ax
        return fig, ax

    def scatter(self, r2:bool=False):

        """
        Generates a scatter plot from data
        :return fig, ax: Matplotlib figure and axis objects.
        """

        data = self.data
        xtick = self.xtick
        ytick = self.ytick
        xrange = self.xrange
        yrange = self.yrange

        labels = self.labels
        desc = self.desc

        colors = self.colors if self.colors is not None else ['b', 'r', 'g', 'c', 'm', 'y', 'k']
        print(colors)

        data_x = data[:, 0]
        data_ys = []

        fig, ax = plt.subplots(1,1, figsize=(5,5))
        ax.set_title(labels[0] if labels is not None else "Scatter Plot")
        ax.set_xlabel(labels[1] if labels is not None else "")
        ax.set_ylabel(labels[2] if labels is not None else "")


        if xtick is not None: ax.set_xticks(xtick)
        if ytick is not None: ax.set_yticks(ytick)
        if xrange is not None: ax.set_xlim(xrange)
        if yrange is not None: ax.set_ylim(yrange)

        for col in range(1, len(data[0,:])):

            data_y = data[:,col]
            data_ys.append(data_y)

            fit = np.polyfit(data_x, data_y, 1)
            val = np.polyval(fit, data_x)

            if r2:
                corr = np.corrcoef(data_x, data_y)[0,1]
                r2_val = corr**2
            else:
                r2_val = None

            ax.scatter(data_x, data_y, marker='.', label=desc[col-1], color = colors[col-1])
            ax.plot(data_x, val, label=f"R2 = {r2_val:.2f}", color = colors[col-1])

        fig.tight_layout()
        ax.legend(loc='best')
        self.fig = fig
        self.ax = ax
