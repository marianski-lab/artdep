from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import colormaps
import copy

import matplotlib as mpl
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.ticker import FuncFormatter
from numpy.ma.extras import average

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

    def cmap(self, color_num: int = None, offset: float = 0, map: str = 'ice'):
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

        if color_num is not None:
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
        if color_num is not None:
            discrete_colors = np.linspace(0, len(colormap_colors) - 1, color_num, dtype=int)
            self.colors = [colormap_colors[i] for i in discrete_colors]
        else:
            self.colors = colormap_colors.tolist()

    def trajectory(self, mol_list, var_name = 'colvar', col = 1, average:int = 0, title=None, hist=True):
        """ Plots MD trajectory with histogram. Takes in data for CP2K or Gromacs via Mol.
        :param molecule: (Mol / List) Either a Mol object, or a list of moles if you want to overlay data. 
        :param var_name: (list) Name of the collective variable you are plotting on your y-axis.
        :param col: (int) Index of the column containing your colvar data, in the case that you have multiple.
        """
        
        if not isinstance(mol_list, list):
            mol_list = [mol_list]
        
        self.path = mol_list[0].path
        
        # CP2K default timestep unit is in fs, Gromacs is in ps:
        # We convert these to ps and nm respectively:

        if mol_list[0].software == 'cp2k':
            time_unit = 'ps'

        elif mol_list[0].software == 'gromacs':
            time_unit = 'ns'
        
        fig, ax = plt.subplots(1,2, figsize=(11,3), gridspec_kw={'width_ratios': [3.5, 1]})
        
        i = 1
        
        for mol in mol_list:
            
            time = (mol.data[:, 0] / 1000).tolist()  # fs -> ps for CP2K, ps -> ns for GROMACS
            colvar = mol.data[:, col].tolist()

            timestep = np.abs(time[0] - time[1])
            
            color = self.colors

            if average > 1:
                array_len =  len(colvar)
                conv_kernel = np.ones(average)/array_len
                colvar = np.convolve(colvar, conv_kernel, mode='valid').tolist()

                time = time[:-1*average + 1]

            ax[0].plot(time, colvar, linewidth=0.2, color=color[i], alpha=0.8)
            ax[1].hist(colvar, bins='rice', fc=(0, 0, 1, 0.5), orientation="horizontal", color=color[i], alpha=0.5)
            
            if len(mol_list) == 1:
                ax[1].set_title(f"average = {np.round(np.average(colvar), 3)}", fontsize = 10)

            else:
                print(f"mol{i} average = {np.round(np.average(colvar), 3)}")

            i = i+1
        
        ax[0].set_xlabel(f"time ({time_unit}); stepsize = {timestep}{time_unit}")
        ax[0].set_ylabel(var_name)
        
        if title != None:
            ax[0].set_title(f"{title}", fontsize = 10)
            
        if hist == False:
            fig.delaxes(ax[1])

        xmax = ax[0].get_xlim()[1]
        ax[0].set_xlim(0, xmax)
        
        # midpt = int(np.round(len(colvar) / 2))
        # ax[1].hist(colvar[0:midpt], bins='rice', fc=(0, 0, 1, 0.3), orientation="horizontal") # First half shown in blue
        # ax[1].hist(colvar[midpt:-1], bins='rice', fc=(0, 0, 1, 0.5), orientation="horizontal") # Second half shown in red
        # ax2.axhline(y=np.average(colvar), color='b', linewidth=2)

        ax[1].set_xlabel('structures')

        plt.tight_layout()
        
        self.fig = fig
        self.ax = ax
        
    def fes(self, mol, cols=[1,2], temp=300, num_levels = 8, num_ticks = 8):
        """ Plots MD FES. Takes in data for CP2K or Gromacs via Mol.
        :param mol: (Mol) Class Mol.
        :param cols: (int) Index of the 2 columns containing your colvar data, in the case that you have more than 2.
        """
        
        self.path = mol.path
        Temp = temp ; R = 8.314 # J/K mol

        colvar1 = mol.data[:, cols[0]]
        colvar2 = mol.data[:, cols[1]]

        Hall, x_edges, y_edges = np.histogram2d(colvar1, colvar2, bins=72)

        Hall = - R * Temp * np.log(Hall)
        hmin = np.min(Hall)
        
        Hall_rel = 0.001*(Hall.T-hmin)

        vmin, vmax = 0, np.ceil(np.nanmax(Hall_rel[~np.isinf(Hall_rel)]))
        MHall = np.ma.masked_greater(Hall_rel, vmax)

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_aspect('equal', adjustable='box')
        
        # colors = self.colors
        # cmap = ListedColormap(colors)

        num_levels = num_levels
        plot = ax.contourf(x_edges[:-1], y_edges[:-1], MHall, cmap='Blues_r', zorder=1, levels=num_levels)
        # plot = ax.contourf(x_edges[:-1], y_edges[:-1], Hall.T, cmap=cmap, zorder=1, levels=num_levels)

        num_ticks = num_ticks
        x_ticks = np.linspace(x_edges[0], x_edges[-1], num_ticks)
        y_ticks = np.linspace(y_edges[0], y_edges[-1], num_ticks)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_xticklabels([f"{tick:.1f}" for tick in x_ticks], fontsize=12)
        ax.set_yticklabels([f"{tick:.1f}" for tick in y_ticks], fontsize=12)

        cb_ticks = np.linspace(vmin, vmax, 6)
        cb = fig.colorbar(plot, ax=ax, ticks=cb_ticks, pad=0.05, shrink=0.78)
        cb.ax.set_yticklabels([f"{tick:.1f}" for tick in cb_ticks], fontsize=12)
        cb.set_label("\n Free energy [kJ]", fontsize=14)

        # Enable grid that aligns with ticks
        ax.grid(True, ls='--', zorder=10.0)

        # Axis labels and title
        ax.set_xlabel("colvar1", fontsize=14)
        ax.set_ylabel("colvar2", fontsize=14)

        # fig.tight_layout()
        
        self.fig = fig
        self.ax = ax
    
    def puckers_hist(self, mol_pucker, mol_fem, puckers=['1C4', '1,4B'], limit=16, temp=300):
        """ Plots ring pucker free energy surface. Requires 2 mol objects to run.
        :param mol_pucker: (Mol) Class Mol containing the .xvg file for your ring pucker determination.
        :param mol_fem: (Mol) Class Mol containing the .xvg file for your free energy surface.
        """
        
        self.path = mol_fem.path
        
        def ring_pucker_determination(mol):

            data = copy.deepcopy(mol.data)
            
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
    
        Temp = temp ; R = 8.314 # J/K mol

        Hall, x_edge, y_edge = np.histogram2d(phi, psi, bins=72, range=[[-180, 180.0],[-180.0, 180.0]])
        #hmax = max(full_data[:,:])

        Hall = - R * Temp * np.log(Hall)
        hmin = np.min(Hall)
        
        if limit == None:
            Hall_rel = 0.001*(Hall.T-hmin)
            limit = np.ceil(np.nanmax(Hall_rel[~np.isinf(Hall_rel)]))
            
        else:
            limit = limit

        Hpuck, edges  = np.histogramdd((phi, psi, puck), bins=[72,72,38], range=[[-180.0, 180.0],[-180.0, 180.0],[0,38]])
        for i in range(38):
            puckers_sum[i] = np.sum(Hpuck[:,:,i])

        for i in range(38): 
            print("{0:4s}{1:10g}".format(id_to_puck(i), puckers_sum[i]))
            
        Hpuck = - R* Temp * np.log(Hpuck)
        hmin_puck = np.min(Hpuck)
        # Hpuck[0,0,] = hmin #To get colorbar right
        
        
        MHall = np.ma.masked_greater(0.001*(Hall.T-hmin), limit-1)
        Mat = [MHall]
        titles = ['All Puckers']
        
        for p in puckers:
            pid = puck_to_id(p)
            titles.append(p)
            MHpuck = np.ma.masked_greater(0.001*(Hpuck[:,:,pid].T - hmin_puck), limit-1)
            # MHpuck[0,0] = hmin #To get colorbar right

            Mat.append(MHpuck)

        fig, axes = plt.subplots(1,len(Mat), figsize=(4*len(Mat) + (len(Mat)-1)*1,  4), sharex=True, sharey=True)

        color_bar = ['Blues_r']*len(Mat)
        levels = np.linspace(0, limit, 9)  # 8 levels between 0 and limit

        color = self.colors
        color.reverse()
        cmap = ListedColormap(color)

        for n, ax in enumerate(axes):
            ax.set_aspect('equal', adjustable='box')
            ax.set_title(titles[n])
            ax.grid(True, ls='--', zorder=10.0)

            # Set x-axis ticks and labels
            xmin, xmax = -180.0, 180.0
            xticks = np.linspace(xmin, xmax, 7)
            ax.set_xticks(np.linspace(0, 71, 7))
            ax.set_xticklabels(['{0:d}'.format(int(x)) for x in xticks], fontsize=12)
            ax.set_xlabel(r'$\phi$', fontsize=14)

            # Set y-axis ticks and labels
            ymin, ymax = -180.0, 180.0
            yticks = np.linspace(ymax, ymin, 7)[::-1]
            ax.set_yticks(np.linspace(0, 71, 7))
            ax.set_yticklabels(['{0:d}'.format(int(x)) for x in yticks], fontsize=12)

            if n == 0:
                ax.set_ylabel(r'$\psi$', fontsize=14)

            # Create the contourf plot with consistent levels
            plot = ax.contourf(Mat[n], levels=levels, cmap=cmap, zorder=1)

            # Add a color bar with consistent boundaries and ticks
            cb = fig.colorbar(plot, ax=ax, pad=0.025, aspect=20, ticks=levels)
            cb.set_ticklabels(["{0:3.1f}".format(x) for x in levels])


        fig.tight_layout()
        
        self.fig = fig
        self.ax = axes
        
    def puckers_scatter(self, mol, puckers=['1C4', '1,4B']):
        scatter_data = mol.data

        ncol = len(scatter_data[0])

        pucker_data = scatter_data[:, ncol - 1] # Puckers are always last.
        energy = scatter_data[:, ncol - 2].astype(float) # E is always 2nd last.
        psi_phi = scatter_data[:, 1:ncol - 2].astype(float) # Psi_Phi is everything in between. 
        conf = scatter_data[:, 0] # _id is always first.

        psi = psi_phi[:, 1::2]
        phi = psi_phi[:, ::2]

        no_datasets = int(np.shape(psi[1])[0])
        
        puckers.insert(0, 'All Puckers')
        filters = []

        for desired_pucker in puckers:
            if desired_pucker == 'All Puckers':
                pucker_filter = np.full(np.shape(pucker_data), True)
                filters.append(pucker_filter)

            else:
                pucker_filter = np.where(pucker_data == desired_pucker, True, False)
                filters.append(pucker_filter)

        cmaps = ['summer','autumn']

        fig, axes = plt.subplots(1,len(puckers), figsize=(4*len(puckers) + (len(puckers)-1)*1,  4), sharex=True, sharey=True)

        for n, ax in enumerate(axes):
            ax.set_aspect('equal', adjustable='box')
            ax.set_title(puckers[n])
            ax.grid(True, ls='--', zorder=10.0)

            # Set x-axis ticks and labels
            xmin, xmax = -180.0, 180.0
            xticks = np.arange(xmin, xmax + 60, 60)  # Ticks every 60 degrees
            ax.set_xticks(xticks)
            ax.set_xticklabels([f'{int(tick)}' for tick in xticks], fontsize=12)
            ax.set_xlabel(r'$\phi$', fontsize=14)

            # Set y-axis ticks and labels
            ymin, ymax = -180.0, 180.0
            yticks = np.arange(ymin, ymax + 60, 60)  # Ticks every 60 degrees
            ax.set_yticks(yticks)
            ax.set_yticklabels([f'{int(tick)}' for tick in yticks], fontsize=12)

            if n == 0:
                ax.set_ylabel(r'$\psi$', fontsize=14)

            for i in range(no_datasets):
                x = psi[:, i][filters[n]]
                y = phi[:, i][filters[n]]
                z = energy[filters[n]]
                ax.scatter(x, y, c=z, cmap=cmaps[i])

            # Set axis limits
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)

        self.fig = fig
        self.ax = axes

    def rdf(self, mol, xmin = 0, xmax=10):
        """ Plots radial distribution function.
        :param mol: (Mol) Class Mol generated from xvg file
        :param xmin: (Int) Minimum x value for rdf plot
        :param xmax: (Int) Maximum x value for rdf plot
        """
    
        self.path = mol.path
        data = copy.deepcopy(mol.data)
        data[:,0] = mol.data[:,0] * 10

        blues  = ['#deebf7','#9ecae1','#3182bd']
        reds   = ['#fee0d2','#fc9272','#de2d26']
        greens = ['#e5f5e0','#a1d99b','#31a354']

        fig, ax = plt.subplots(figsize=(6,2))

        xmin = xmin; xmax = xmax

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

        color = self.colors

        ax.plot(data[:, 0],  np.convolve(data[:,2], np.ones(5)/5, mode='same')/div, color=color[1])
        ax.fill_between(data[:, 0],  np.convolve(data[:,2], np.ones(5)/5, mode='same')/div, color=color[2])

        #HIS
        div = np.amax(data[:,1]) 
        if div == 0: div = 1

        ax.plot(data[:, 0],  np.convolve(data[:,1], np.ones(5)/5, mode='same')/div, color=reds[2])
        ax.fill_between(data[:, 0],  np.convolve(data[:,1], np.ones(5)/5, mode='same')/div, color=reds[0])
        
        fig.tight_layout()

        self.fig = fig
        self.ax = ax
        
    def contour(self, xpm_mols, limit = 16):
        """ Plots contour maps for a provided list of moles from xpm files.
        :param xpm_mols: (List) List of mol objects generated from xpm files.
        :param limit: (Int) The upper limit on your energy scale
        """
        
        self.path = xpm_mols[0].path
        limit = limit
        
        Mat = []
        
        for xpm_mol in xpm_mols:
            M = xpm_mol.data
            MM = np.ma.masked_greater(M, limit-1)
            Mat.append(MM)
        
#         if xpm_mol1 != None and xpm_mol2 != None and xpm_mol3 != None: 
#             M1 = xpm_mol1.data
#             M2 = xpm_mol2.data
#             M3 = xpm_mol3.data
#             MM1 = np.ma.masked_greater(M1, limit-1)
#             MM2 = np.ma.masked_greater(M2, limit-1) 
#             MM3 = np.ma.masked_greater(M3, limit-1) 
#             Mat = [MM1, MM2, MM3]

#         elif xpm_mol1 != None and xpm_mol2 != None:
#             M1 = xpm_mol1.data
#             M2 = xpm_mol2.data
#              #DiffM = M1-M2
#             MM1 = np.ma.masked_greater(M1, limit-1)
#             MM2 = np.ma.masked_greater(M2, limit-1)  
#             Mat = [MM1, MM2]

#         else:
#             raise ValueError("Give me a correct number of some chunky matrices")
        fig, axes = plt.subplots(1,len(Mat), figsize=(4*len(Mat) + (len(Mat)-1)*1,  4), sharex=True, sharey=True)
        # fig, axes = plt.subplots(1,len(Mat), figsize=(4*len(Mat) + 1.5, 4), sharex=True, sharey=True)
        
        color_bar = ['Blues_r']*len(Mat)
        color = self.colors
        color.reverse()
        cmap = ListedColormap(color)

        levels = np.linspace(0, limit, 9)  # 8 levels between 0 and limit

        for n, ax in enumerate(axes):
            ax.set_aspect('equal', adjustable='box')
            ax.grid(True, ls='--', zorder=10.0)

            # Set x-axis ticks and labels
            xmin, xmax = -180.0, 180.0
            xticks = np.linspace(xmin, xmax, 7)
            ax.set_xticks(np.linspace(0, 71, 7))
            ax.set_xticklabels(['{0:d}'.format(int(x)) for x in xticks], fontsize=12)
            ax.set_xlabel(r'$\phi$', fontsize=14)

            # Set y-axis ticks and labels
            ymin, ymax = -180.0, 180.0
            yticks = np.linspace(ymax, ymin, 7)[::-1]
            ax.set_yticks(np.linspace(0, 71, 7))
            ax.set_yticklabels(['{0:d}'.format(int(x)) for x in yticks], fontsize=12)

            if n == 0:
                ax.set_ylabel(r'$\psi$', fontsize=14)

            # Create the contourf plot with consistent levels
            plot = ax.contourf(Mat[n], levels=levels, cmap=cmap, zorder=1)

            # Add a color bar with consistent boundaries and ticks
            cb = fig.colorbar(plot, ax=ax, pad=0.025, aspect=20, ticks=levels)
            cb.set_ticklabels(["{0:3.1f}".format(x) for x in levels])

        fig.tight_layout()
        
        self.fig = fig
        self.ax = axes

    def scatter(self):

        """
        Generates a scatter plot from data
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

        ax.tick_params(axis='both', which='both', bottom=True, top=False, labelbottom=True, right=False, left=True,
                       labelleft=True)
        for s in ['top', 'right', 'left', 'bottom']: ax.spines[s].set_visible(False)

        ax.xaxis.set_tick_params(direction='out');
        ax.yaxis.set_tick_params(direction='out')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        for col in range(1, len(data[0,:])):

            data_y = data[:,col]
            data_ys.append(data_y)

            fit = np.polyfit(data_x, data_y, 1)
            val = np.polyval(fit, data_x)

            ax.scatter(data_x, data_y, marker='.', label=desc[col-1], color = colors[col-1])

        xrange = list(ax.get_xlim()) if xrange is None else xrange
        yrange = list(ax.get_ylim()) if yrange is None else yrange



        if xrange is not None:
            xrange[0] -= self.x_extend
            xrange[1] += self.x_extend
            ax.set_xlim(xrange)
        if yrange is not None:
            yrange[0] -= self.y_extend
            yrange[1] += self.y_extend
            ax.set_ylim(yrange)

        xtick = list(ax.get_xticks())
        ytick = list(ax.get_yticks())

        minx = round(xtick[1], 1)
        maxx = round(xtick[-2], 1)
        miny = round(ytick[1] ,1)
        maxy = round(ytick[-2] ,1)

        print(minx, maxx, miny, maxy)

        ax.plot([minx-0.05, maxx+0.05], [yrange[0], yrange[0]], color='k')
        ax.plot([xrange[0], xrange[0]], [miny-0.05, maxy+0.05], color='k')

        fig.tight_layout()
        ax.legend(bbox_to_anchor=(-0.5, 0.5), loc='center left', borderaxespad=0, frameon=False)
        self.fig = fig
        self.ax = ax

    def reaction_profile(self, reaction, type, linewidth=3, scale=0.32, annotate=True, color='blue'):
        """
        Plots a reaction coordinate diagram.
        """
        mol_list = reaction.mol_list
        labels = reaction.mol_label
        
        energies = []

        '''
        type = 'delta E' or 'delta F' or 'delta H'
        '''

        for mol in mol_list:
            if type == 'delta E':
                energies.append(mol.E)   
            elif type == 'delta F':
                energies.append(mol.F) 
            elif type == 'delta H':
                energies.append(mol.H) 
            else:
                print("Unsupported Energy Type")
                return  
             
        if not energies:
            raise ValueError("No energies found. Check the input data.")

        # Dynamic Figure Size Based on Number of Reaction Steps
        num_steps = len(mol_list)
        fig_width = max(6, num_steps * 2)  # Adjust width based on number of steps
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        

        relative_energies = [hartree_to_kcal(e - energies[0]) for e in energies]

        annotation_offset = 0.13

        for j, energy in enumerate(relative_energies):
            # Draw Horizontal Bars at Each Energy Level
            ax.plot([(j + 1 - scale), (j + 1 + scale)], [energy, energy],
                    color=color, linewidth=linewidth)

            # Annotate Energy Values
            if annotate and j != 0:
                ax.text(j + 1, energy + annotation_offset, f"{energy:.2f}", fontsize=12, ha='center', color='black')

            # Draw Dashed Connecting Lines
            if j < len(relative_energies) - 1:
                ax.plot([(j + 1 + scale), (j + 2 - scale)],
                        [energy, relative_energies[j + 1]],
                        linestyle=":", color=color, linewidth=linewidth)

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

        
        ax.set_xticks(range(1, len(energies) + 1))
        ax.set_xticklabels(labels) 
            
        
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_linewidth(1.5)  
        ax.spines['left'].set_linewidth(1.5)    

        # Final Formatting
        ax.tick_params(labelsize=14)
        ax.legend(loc="lower left", frameon=False, fontsize=14)

        self.fig = fig;
        self.ax = ax
        
    def savefig(self, filename='fig', format:str='png'):
        self.fig.savefig(f"{self.path}/{filename}.{format}", dpi=300, bbox_inches='tight')
