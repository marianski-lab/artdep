import numpy as np
import matplotlib.pyplot as plt
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

    def __init__(self, data, labels:list = None, desc:list = None, xtick:int = None, ytick:int = None, xrange:list = None, yrange:list = None, colors:list = None) :
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

        self.fig = fig
        self.ax = ax




