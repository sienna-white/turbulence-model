import numpy as np
import matplotlib.pyplot as plt
from constants import *
import matplotlib as mpl
#
# ## PLOTTING UTILITIES ### 


def plot_shear(self, skip, passed_string='', show=True, ):
    fig, ax = plt.figure(figsize=(6,6)), plt.gca()
    plot_index, legend_ind = self.get_plot_indices (skip)
    velocity = self.dataset['U'].values 
    shear = np.gradient(velocity, axis=0)/self.dz

    # assert(False)
    ls = '-'
    z = self.z
    for i, ind in reversed(list(enumerate(plot_index))):
            time = self.dataset.time.values[ind]
            if time == 0:
                    ax.plot(shear[:,0], 
                    z, '--',
                    color = 'k',
                    linewidth = 2, 
                    label='Initial condition')
            else:
                if i%legend_ind==0: 
                    ax.plot(shear[:,ind], 
                            z, ls,
                            color = mpl.cm.viridis(i/len(plot_index)),
                            linewidth = 2.5, 
                            alpha = 0.6,
                            label='t = %2.1f hr' % seconds2hours(time))
                else:
                    ax.plot(shear[:,ind], 
                            z, ls,
                            color = mpl.cm.viridis(i/len(plot_index)),
                            linewidth = 2.5, 
                            alpha = 0.6)
        
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax.grid(alpha = 0.5)
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel("1/s")
    ax.set_xlim(-0.1, 0.1)
    ax.set_title("Shear " + passed_string)

    plt.tight_layout()
    if show:
        plt.show()
    return fig, ax  

def seconds2hours(seconds):
    return seconds/3600

def get_plot_indices(self, skip):
    plots = np.arange(0, len(self.dataset.time), skip)
    if len(plots)>10: 
        legend_ind = np.floor(len(plots)/10)
    else:
        legend_ind = 1
    return plots, legend_ind
    
def plot_profiles(self, variable, skip=1, passed_string='', show=True, ax=None):
    
    if ax is None: 
        fig, ax = plt.figure(figsize=(8,4)), plt.gca()
    else: 
        fig = plt.gcf()
    plot_index, legend_ind = self.get_plot_indices(skip)
    ls = '-'

    self.add_profile_to_axis(variable, plot_index, legend_ind, ax)
    ax.set_xlabel(variable2units[variable])
    ax.set_title(variable2name[variable] + ' ' +  passed_string)
    plt.tight_layout()
    if show:
        plt.show()
    return fig, ax  

def add_time_series_to_axis(self, variable, label, ax):
    ''' Add a time series to a given axis. ''' 
    ds = self.dataset[variable].dropna(dim="time") # Drop nan values 
    if label is None:
        label = variable2name[variable] 

    ax.plot(seconds2hours(ds.time), ds.values, '-o', markersize = 5, label = label)
    ax.grid(alpha = 0.5)
    ax.set_xlabel('Time (hr)')
    return ax 

def plot_time_series(self, variable, passed_string, show=True):
    fig, ax = plt.figure(figsize=(8,4)), plt.gca()
    ax = self.add_time_series_to_axis(variable, None, ax)
    ax.set_title(variable2name[variable])
    ax.legend()
    plt.tight_layout()
    if show:
        plt.show()
    return fig, ax

# def plot_profiles(self, variable, skip=1, passed_string='', show=True):
#     ''' plot profiles ''' 
#     fig, ax = plt.figure(figsize=(8,4)), plt.gca()
#     plot_index, legend_ind = self.get_plot_indices(skip)
#     ls = '-'
#     self.add_profile_to_axis(variable, plot_index, legend_ind, ax)
#     ax.set_xlabel(variable2units[variable])
#     ax.set_title(variable2name[variable] + ' ' +  passed_string)
#     plt.tight_layout()
#     if show:
#         plt.show()
#     return fig, ax  

def add_profile_to_axis(self, variable, plot_index, legend_ind, ax):
    ''' Add a vertical profile given by the key variable to a given axis'''
    ls  = '-'
    ds = self.dataset[variable] * variable2convert[variable]
    for i, ind in reversed(list(enumerate(plot_index))):
        time = ds.time.values[ind]
        if time == 0:
                ax.plot(ds.isel(time=0), 
                ds.z, '--',
                color = 'k',
                linewidth = 2, 
                label='Initial condition')
        else:
            if i%legend_ind==0: 
                ax.plot(ds.isel(time=ind), 
                        self.z, ls,
                        color = mpl.cm.viridis(i/len(plot_index)),
                        linewidth = 2.5, 
                        alpha = 0.6,
                        label='t = %d hr' % seconds2hours(time))
            else:
                ax.plot(ds.isel(time=ind), 
                        self.z, ls,
                        color = mpl.cm.viridis(i/len(plot_index)),
                        linewidth = 2.5, 
                        alpha = 0.6)
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax.grid(alpha = 0.5)
    ax.set_ylabel('Depth (m)')
    
    return ax 


def plot_phasing(self, ListOfSpecies, pressure, light, wind, passed_string='', skip=1, show=True):
    ''' 
        1st row:    biomass 1, biomass 2
        2nd row :   concentration 1, concentration 2 
        3rd row:    diurnal light, tidal forcing 
    '''

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize = (18, 10))
    plot_index, legend_ind = self.get_plot_indices(skip)
    axs = axs.flatten()

    # Add biomass 
    axs[0].set_title("Biomass over time")
    label1 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[0].name, ListOfSpecies[0].pmax, ListOfSpecies[0].ws)
    label2 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[1].name, ListOfSpecies[1].pmax, ListOfSpecies[1].ws)
    axs[0] = self.add_time_series_to_axis('biomass1', label1, axs[0])
    axs[0] = self.add_time_series_to_axis('biomass2', label2, axs[0])
    axs[0].set_ylim(0,1)
    axs[0].legend()

    # Add forcings
    time = self.get_time_steps()/3600 
    length = len(time)

    # Add concentrations 
    axs[1] = self.add_profile_to_axis( 'algae1', plot_index, legend_ind, axs[1])
    axs[1].set_title("Concentration of %s" % ListOfSpecies[0].name)

    axs[3] = self.add_profile_to_axis( 'algae2', plot_index, legend_ind, axs[3])
    axs[3].set_title("Concentration of %s" % ListOfSpecies[1].name)


    axs[2].plot(time, light/max(light), label='Light', linewidth = 5, color='#f3b151')
    axs[2].plot(time,pressure/max(pressure), label='Pressure', linewidth = 3, color='#0e2d74')
    axs[2].plot(time, wind/max(wind), '--', label='Wind', linewidth = 2, color="#423af3")

    axs[2].legend()
    axs[2].grid(alpha = 0.5)
    axs[2].set_title("Normalized temporal forcings")
    axs[2].set_ylabel("[-]")
    axs[2].set_xlabel("Time (hr)")

    return fig, axs