import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
import os
import sys
from .utilities import *
from .mol import *
from MM_Util import mol

class Plot():
    """
    A class that creates plots from a numpy array.
    It can make reaction coordinate diagrams, RMSD and RMSF trajectory plots,
    Free energy plots, Scatter Plots, and SNFG Figures.
    """

    def __init__(self, data):
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
        
        return fig, ax
    
    
    def reaction_energies(labels, type, energies1, energies2, mult=False):
    # creating the numpy arrays containing the energy data for one or two reactions
        labels = []    # The name of the reaction step (R, Int, TS, P)
        name_1 = []     # The name of the reaction
        name_2 = []   # the name of a second reaction
        energies1 = []
        energies2 = []
        if mult and type=='delta E':
            for item in Mol:
                energies1.append(mol.E(item))
                labels.append(item)
                energies2.append(item)

        elif mult==False and type=='delta E':
            for item in Mol:
                energies1.append(mol.E(item))
                labels.append(mol.label(item))

        elif mult and type=='delta F':
            for item in Mol:
                energies1.append(mol.F(item))
                labels.append(mol.label(item))
                energies2.append(mol.F(item))

        elif mult==False and type=='delta F':
            for item in Mol:
                energies1.append(mol.F(item))
                labels.append(mol.label(item))

        elif mult and type=='delta H':
            for item in Mol:
                energies1.append(mol.H(item))
                labels.append(mol.label(item))
                energies2.append(mol.H(item))

        elif mult==False and type=='delta H':
            for item in Mol:
                energies1.append(mol.H(item))
                labels.append(mol.label(item))

        elif type is not 'delta H' or 'delta F' or 'delta E':
            print(f'Unsupported Energy format')



    def plot_reaction(ax, energies, color, label, linewidth=3, scale=0.32, annotate=True, first_black=False, other_energies=None):

        upper_offset = 0.10  # Offset for upper bars
        lower_offset = -0.40  # Offset for lower bars
        
        for j, energy in enumerate(energies):
            # Skip the black bar for the first coordinate if there's only one reaction
            if other_energies is None:
                bar_color = color  # Use the reaction color for all bars
            else:
                bar_color = "black" if first_black and j == 0 else color

            # Draw Horizontal Bars at Each Energy Level
            ax.plot([(j + 1 - scale), (j + 1 + scale)], [energy, energy], 
                color=bar_color, linewidth=linewidth, zorder=3 if j == 0 and first_black else 2)   

            # Determine if this is the lower or upper bar (for two reactions)
            if other_energies is not None and j < len(other_energies):
                is_lower_bar = energy < other_energies[j]
            else:
                is_lower_bar = True  # Assume current energy is the lower bar if other_energies is shorter

            # Calculate annotation offset based on whether this is the lower or upper bar
            if other_energies is None:
                text_y_offset = lower_offset if is_lower_bar else upper_offset

            # Annotate Energy Values
            ax.text(j + 1, energy + text_y_offset, proper_minus(energy), 
                color='black', fontsize=12, ha='center') if j!=0 else None

            # Draw Dashed Connecting Lines
            if j < len(energies) - 1:
                ax.plot([(j + 1 + scale), (j + 2 - scale)], 
                    [energy, energies[j + 1]], 
                    linestyle=":", color=color, linewidth=linewidth)

        # Invisible plot for the legend label
        ax.plot([], [], color=color, linewidth=linewidth, label=label)

    def plot_reaction_coords(reaction_energies, plot_reaction, energies1,energies2,labels, name_1, name_2):
        relative_energies_1 = [hartree_to_kcal(energies1[0] - e) for e in energies1]
        relative_energies_2 = [hartree_to_kcal(energies2[0] - e) for e in energies2] if energies2 else []

        if type == 'delta E':
            reaction_type = '$\\Delta E$ (kcal $\\cdot$ mol${}^{-1}$)' 
        elif type == 'delta F':
            reaction_type = '$\\Delta F$ (kcal $\\cdot$ mol${}^{-1}$)'
        else:
            reaction_type = None

        # Dynamic Figure Size Based on Number of Reaction Steps
        num_steps = len(labels)
        fig_width = max(6, num_steps * 2)  # Adjust width based on number of steps
        fig, ax = plt.subplots(figsize=(fig_width, 6))

        # Plot the Reaction Coordinate Diagram
        if relative_energies_2:  # If two reactions, set first point of reaction 1 to black
            plot_reaction(ax, relative_energies_1, color_picker('1a'[0]), name_1, first_black=True, other_energies=relative_energies_2)
            plot_reaction(ax, relative_energies_2, color_picker('1a'[1]), name_2, other_energies=relative_energies_1)
        else:
            plot_reaction(ax, relative_energies_1, color_picker('1a'[0]), name_1, first_black=True)

        # Set Y-axis Label
        ax.set_ylabel(reaction_type, fontsize=16)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: proper_minus(x)))

        # Set X-axis Ticks and Labels
        x_ticks = list(range(1, len(relative_energies_1) + 1))
        if len(labels) == len(x_ticks):
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(labels)
        else:
            print(f'Error: reaction_steps ({len(labels)}) and x_ticks ({len(x_ticks)}) do not match')

        # Add X-axis Guide Line at the halfway point
        ax.axhline(0, color="black", linestyle=":", linewidth=1.5, zorder=-4)

        # Adjust Y-axis Limits
        y_min, y_max = min(relative_energies_1 + relative_energies_2), max(relative_energies_1 + relative_energies_2)
        y_range = y_max - y_min
        ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.2 * y_range)

        # Customize Borders
        for spine in ax.spines.values():
            spine.set_visible(False)  # Hide all borders initially

        # Make bottom and left borders thicker and visible
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_linewidth(1.5)  # Thicker bottom border
        ax.spines['left'].set_linewidth(1.5)    # Thicker left border

        # Final Formatting
        ax.tick_params(labelsize=14)
        ax.legend(loc="lower left", frameon=False, fontsize=14)