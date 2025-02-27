import numpy as np
import matplotlib.pyplot as plt
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