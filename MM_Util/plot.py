import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import colormaps

from matplotlib.ticker import FormatStrFormatter

from .utilities import *
from .mol import *

class Plot():
    """
    A class that creates plots from a numpy array.
    It can make reaction coordinate diagrams, RMSD and RMSF trajectory plots,
    Free energy plots, Scatter Plots, and SNFG Figures.
    """

    def __init__(self, data, labels:list = None, desc:list = None,
                 xtick = None, ytick = None, xrange:list = None, yrange:list = None,
                 colors:list = None, x_extend:float=0, y_extend:float=0) :

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
        self.x_extend = x_extend
        self.y_extend = y_extend

    def cmap(self, color_num: int = 256, offset: float = 0, map: str = 'ice'):
        """
        Generates and processes a colormap with optional offsetting logic.
        :param color_num: (int) Number of discrete colors.
        :param offset: (float) Fractional offset to shift the colormap.
        :param map: (str) Name of the colormap from the colormaps library.
        """
        # Check if the colormap exists in colormaps
        if not hasattr(colormaps, map):
            raise ValueError(f"Colormap '{map}' not found in colormaps library!")

        # Fetch colormap
        colors_obj = getattr(colormaps, map)
        color_num += 1
        # Ensure the colormap has an array of colors
        if not hasattr(colors_obj, 'colors'):
            raise ValueError(f"The selected colormap '{map}' does not have a valid 'colors' attribute!")

        colormap_colors = colors_obj.colors

        # Validating the shape of colormap_colors
        if len(colormap_colors[0]) != 3:
            raise ValueError(f"Expected RGB colors in the colormap, but got shape {np.array(colormap_colors).shape}.")

        # Applying offset manually

        if offset != 0:
            new_colors = []

            for color in colormap_colors:

                new_color = []
                for color_elm in color:
                    color_elm -= offset

                    if color_elm > 1:
                        color_elm = 1

                    if color_elm < 0:
                        color_elm = 0

                    new_color.append(color_elm)
                new_colors.append(new_color)

            colormap_colors = new_colors

        # Discretize the colormap to the required number of colors
        discrete_colors = np.linspace(0, len(colormap_colors) - 1, color_num, dtype=int)
        self.colors = [colormap_colors[i] for i in discrete_colors]

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
        return fig, ax

    def scatter(self):

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

        data_x = data[:, 0]
        data_ys = []

        fig, ax = plt.subplots(1,1, figsize=(5,5))
        ax.set_title(labels[0] if labels is not None else "Scatter Plot")
        ax.set_xlabel(labels[1] if labels is not None else "")
        ax.set_ylabel(labels[2] if labels is not None else "")


        if xtick is not None: ax.set_xticks(xtick)
        if ytick is not None: ax.set_yticks(ytick)
        if xrange is not None:
            xrange[0] -= self.x_extend
            xrange[1] += self.x_extend
            ax.set_xlim(xrange)
        if yrange is not None:
            yrange[0] -= self.y_extend
            yrange[1] += self.y_extend
            ax.set_ylim(yrange)

        ax.tick_params(axis='both', which='both', bottom=True, top=False, labelbottom=True, right=False, left=True,
                       labelleft=True)
        for s in ['top', 'right', 'left', 'bottom']: ax.spines[s].set_visible(False)

        ax.xaxis.set_tick_params(direction='out');
        ax.yaxis.set_tick_params(direction='out')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        for col in range(1, len(data[0,:])):

            data_y = data[:,col]
            data_ys.append(data_y)

            fit = np.polyfit(data_x, data_y, 1)
            val = np.polyval(fit, data_x)

            ax.scatter(data_x, data_y, marker='.', label=desc[col-1], color = colors[col-1])

        xrange = ax.get_xlim() if xrange is None else xrange
        yrange = ax.get_ylim() if yrange is None else yrange

        fig.tight_layout()
        ax.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', borderaxespad=0, frameon=False)
        self.fig = fig
        self.ax = ax
