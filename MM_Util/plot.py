import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import colormaps
import copy

from matplotlib.ticker import NullFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

from .utilities import *
from .mol import *

class Plot():
    """
    A class that creates plots from a numpy array.
    It can make reaction coordinate diagrams, RMSD and RMSF trajectory plots,
    Free energy plots, Scatter Plots, and SNFG Figures.
    """

    def __init__(self, data=None, labels:list = None, desc:list = None,
                 xtick = None, ytick = None, xrange:list = None, yrange:list = None,
                 colors:list = ['b', 'r', 'g', 'c', 'm', 'y', 'k'], x_extend:float=0, y_extend:float=0) :

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

        ax[1].hist(colvar, bins='rice', fc=(0, 0, 1, 0.5), orientation="horizontal") 
        
        # midpt = int(np.round(len(colvar) / 2))
        # ax[1].hist(colvar[0:midpt], bins='rice', fc=(0, 0, 1, 0.3), orientation="horizontal") # First half shown in blue
        # ax[1].hist(colvar[midpt:-1], bins='rice', fc=(0, 0, 1, 0.5), orientation="horizontal") # Second half shown in red
        # ax2.axhline(y=np.average(colvar), color='b', linewidth=2)

        ax[1].set_title(f"average = {np.round(np.average(colvar), 3)}", fontsize = 10)

        ax[1].set_xlabel('structures')

        plt.tight_layout()
        
        self.fig = fig
        self.ax = ax
    
    def puckers_hist(self, mol_pucker, mol_fem):
        """ Plots ring pucker free energy surface. Requires 2 mol objects to run.
        :param mol_pucker: (Mol) Class Mol containing the .xvg file for your ring pucker determination.
        :param mol_fem: (Mol) Class Mol containing the .xvg file for your free energy surface.
        """
        
        def ring_pucker_determination(mol):

            data = mol.data
            
            n = data.shape[1] - 1
            angles = data[:, -n:]
            angles = np.where(angles > 0.0, 180.0 - angles, -angles - 180.0)

            data[:, -n:] = angles
            traj_idx = np.array([str(x) for x in data[:, 0]])

            pucker_table = {
                '1C4': [-35.26, -35.26, -35.26], '4C1': [35.26, 35.26, 35.26],
                '1,4B': [-35.26, 74.20, -35.26], 'B1,4': [35.26, -74.20, 35.26],
                '2,5B': [74.20, -35.26, -35.26], 'B2,5': [-74.20, 35.26, 35.26],
                '3,6B': [-35.26, -35.26, 74.20], 'B3,6': [35.26, 35.26, -74.20],
                '1H2': [-42.16, 9.07, -17.83], '2H1': [42.16, -9.07, 17.83],
                '2H3': [42.16, 17.83, -9.06], '3H2': [-42.16, -17.83, 9.06],
                '3H4': [-17.83, -42.16, 9.07], '4H3': [17.83, 42.16, -9.07],
                '4H5': [-9.07, 42.16, 17.83], '5H4': [9.07, -42.16, -17.83],
                '5H6': [9.07, -17.83, -42.16], '6H5': [-9.07, 17.83, 42.16],
                '6H1': [17.83, -9.07, 42.16], '1H6': [-17.83, 9.07, -42.16],
                '1S3': [0.00, 50.84, -50.84], '3S1': [0.00, -50.84, 50.84],
                '5S1': [50.84, -50.84, 0.00], '1S5': [-50.84, 50.84, 0.00],
                '6S2': [-50.84, 0.00, 50.84], '2S6': [50.84, 0.00, -50.84],
                '1E': [-35.26, 17.37, -35.26], 'E1': [35.26, -17.37, 35.26],
                '2E': [46.86, 0.00, 0.00], 'E2': [-46.86, 0.00, 0.00],
                '3E': [-35.26, -35.26, 17.37], 'E3': [35.26, 35.26, -17.37],
                '4E': [0.00, 46.86, 0.00], 'E4': [0.00, -46.86, 0.00],
                '5E': [17.37, -35.26, -35.26], 'E5': [-17.37, 35.26, 35.26],
                '6E': [0.00, 0.00, 46.86], 'E6': [0.00, 0.00, -46.86]
            }

            pucker_table_list = np.array(list(pucker_table.values()))
            pucker_keys = list(pucker_table.keys())
            len_puck = len(pucker_keys)
            pucker = [] 

            # RMSD calculations
            for ring in angles: 

                dist_matrix  = copy.copy(pucker_table_list)
                dist_matrix -= ring
                l1_norm = np.zeros((len_puck,)) ; l2_norm = np.zeros((len_puck,))

                for i in range(len_puck):
                    l1_norm[i] = 0.333 * np.abs(dist_matrix[i,0] + dist_matrix[i,1] + dist_matrix[i,2])
                    l2_norm[i] = 0.333 * np.sqrt(dist_matrix[i,0]**2 + dist_matrix[i,1]**2 + dist_matrix[i,2]**2)
                #print(l1_norm)
                min_dist_values = np.min(l1_norm)
                min_dist_indices = np.argmin(l1_norm)
                pucker.append(pucker_keys[min_dist_indices])

            # Let us store this information as an attribute:
            # np.array(list(zip(traj_idx, pucker)))

            return pucker

        def load_dihedrals(mol):

            data = mol.data

            phi = data[:, ::2][:, 1:] # Even col
            psi = data[:, 1::2] # Odd col

            return phi.flatten(), psi.flatten()

        def puck_to_id(pucker):

            pucker_table = {
                '1C4': 0, '4C1': 1, '1,4B':2, 'B1,4': 3, '2,5B':4, 'B2,5': 5,  '3,6B':6, 'B3,6': 7,
                '1H2': 8, '2H1':  9, '2H3': 10, '3H2': 11, '3H4': 12, '4H3': 13, '4H5': 14, '5H4': 15,
                '5H6': 16, '6H5': 17,'6H1': 18, '1H6': 19, '1S3': 20, '3S1': 21, '5S1': 22, '1S5': 23,
                '6S2': 24, '2S6': 25, '1E': 26, 'E1': 27,  '2E': 28, 'E2': 29, '3E': 30, 'E3': 31,
                '4E': 32, 'E4': 33, '5E': 34, 'E5': 35, '6E': 36, 'E6': 37 
            }
            return pucker_table[pucker]

        def id_to_puck(_id): 

            pucker_table = ['1C4', '4C1', 
                    '1,4B', 'B1,4', '2,5B', 'B2,5','3,6B', 'B3,6',
                    '1H2', '2H1', '2H3', '3H2', '3H4', '4H3', '4H5', '5H4', '5H6', '6H5','6H1', '1H6',
                    '1S3', '3S1', '5S1', '1S5', '6S2', '2S6',
                    '1E', 'E1', '2E', 'E2', '3E', 'E3', '4E', 'E4', '5E', 'E5', '6E', 'E6']

            return pucker_table[_id]

        pucker = ring_pucker_determination(mol_pucker)
        puck = [puck_to_id(p) for p in pucker]
        puckers_sum = np.zeros((38,))
        
        phi, psi  = load_dihedrals(mol_fem)
        
        limit = 16
        puck1 = 0 ; puck2 = 2 # See id_to_puck for definitions
        Temp = 300.0 ; R = 8.314 # J/K mol

        Hall, x_edge, y_edge = np.histogram2d(phi, psi, bins=72, range=[[-180, 180.0],[-180.0, 180.0]])
        #hmax = max(full_data[:,:])

        Hall = - R * Temp * np.log(Hall)
        hmin = np.min(Hall)
        print(hmin)

        Hpuck, edges  = np.histogramdd((phi, psi, puck), bins=[72,72,38], range=[[-180.0, 180.0],[-180.0, 180.0],[0,38]])
        for i in range(38):
            puckers_sum[i] = np.sum(Hpuck[:,:,i])

        for i in range(38): 
            print("{0:4s}{1:10g}".format(id_to_puck(i), puckers_sum[i]))

        Hpuck = - R* Temp * np.log(Hpuck)


        MHall   = np.ma.masked_greater(0.001*(Hall.T-hmin), limit-1)
        MHpuck1 = np.ma.masked_greater(0.001*(Hpuck[:,:,puck1].T-hmin), limit-1)

        Hpuck[0,0,puck2] = hmin #To get colorbar right
        MHpuck2 = np.ma.masked_greater(0.001*(Hpuck[:,:,puck2].T-hmin), limit-1)

        #Mat = [MHall, MHpuck1]
        Mat = [MHall, MHpuck1, MHpuck2]

        fig, axes = plt.subplots(1,len(Mat), figsize=(4*len(Mat) + (len(Mat)-1)*1,  4), sharex=True, sharey=True)
        colors = self.colors
        colors.reverse()

        for n, ax in enumerate(axes):

          vmin = 0 ; vmax = limit #np.amax(Mat[n])
          ax.set_aspect('equal', adjustable='box')
          ax.grid(True, ls='--', zorder=10.0)

          xmin=-180.0; xmax=180.0 
          xticks = np.linspace(xmin, xmax, 7)
          ax.set_xticks(np.linspace(0, 71, 7))
          ax.set_xticklabels(['{0:d}'.format(int(x)) for x in xticks], fontsize=12) #rotation=45, ha="right", rotation_mode="anchor")
          ax.set_xlabel(r'$\phi$', fontsize=14)

          ymin=-180.0; ymax=180.0 
          yticks = np.linspace(ymax, ymin, 7)
          yticks = yticks[::-1]
          ax.set_yticks(np.linspace(0, 71, 7))
          ax.set_yticklabels(['{0:d}'.format(int(x)) for x in yticks], fontsize=12) #rotation=45, ha="right", rotation_mode="anchor")

          if n==0:  ax.set_ylabel(r'$\psi$', fontsize=14)

          cmap = ListedColormap(colors)
          #plot = ax.imshow(Mat[n],  aspect='auto', interpolation='none', cmap=color_bar[n], vmin=vmin, vmax=vmax)
          plot = ax.contourf(Mat[n],  vmin=0, vmax=limit, cmap=cmap, zorder=1, levels=8)
          plot.set_clim(0,limit)
          cb_ticks = np.linspace(0, limit,  int(limit/2)+1)
          #f n == len(Mat)-1: 
          cb = fig.colorbar(plot, ax=axes[n], ticks=cb_ticks, pad=0.025, aspect=20)
          cb.ax.set_yticklabels([ "{0:3.1f}".format(x) for x in cb_ticks], fontsize=12)

        fig.tight_layout()
        
        self.fig = fig
        self.ax = axes
        
    def rdf(self, mol):
    
        data = mol.data
        data[:,0] = mol.data[:,0] * 10

        blues  = ['#deebf7','#9ecae1','#3182bd']
        reds   = ['#fee0d2','#fc9272','#de2d26']
        greens = ['#e5f5e0','#a1d99b','#31a354']

        fig, ax = plt.subplots(figsize=(6,2))

        xmin = 0; xmax=10

        ax.tick_params(axis='both', which='both', bottom=True, top=False, labelbottom=True, right=False, left=False, labelleft=False)
        ax.spines['top'].set_visible(False) ; ax.spines['right'].set_visible(False) ; ax.spines['left'].set_visible(False)
        ax.xaxis.set_tick_params(direction='out')
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.set_ylim(0,1.5)

        xticks = np.linspace(xmin,xmax,int((xmax-xmin/10)+1))
        ax.set_xticks(xticks)
        ax.set_xticklabels([int(x) for x in xticks], fontsize=10)
        ax.set_xlim(xmin, xmax)

        #TYR
        div = np.amax(data[:,3]) 
        if div == 0: div = 1

        ax.plot(data[:, 0],  np.convolve(data[:,3], np.ones(5)/5, mode='same')/div, color=greens[2])
        ax.fill_between(data[:, 0],  np.convolve(data[:,3], np.ones(5)/5, mode='same')/div, color=greens[0])

        #D or E
        div = np.amax(data[:,2]) 
        if div == 0: div = 1

        ax.plot(data[:, 0],  np.convolve(data[:,2], np.ones(5)/5, mode='same')/div, color=blues[2])
        ax.fill_between(data[:, 0],  np.convolve(data[:,2], np.ones(5)/5, mode='same')/div, color=blues[0])

        #HIS
        div = np.amax(data[:,1]) 
        if div == 0: div = 1

        ax.plot(data[:, 0],  np.convolve(data[:,1], np.ones(5)/5, mode='same')/div, color=reds[2])
        ax.fill_between(data[:, 0],  np.convolve(data[:,1], np.ones(5)/5, mode='same')/div, color=reds[0])
        
        fig.tight_layout()

        self.fig = fig
        self.ax = ax

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
