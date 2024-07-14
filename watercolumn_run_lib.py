import pandas as pd 
import xarray as xr 
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd 
import os

['U', 'C', 'Q2', 'Q2L', 'rho', 'L', 'nu_t', 'Kz', 'Kq', 'N_BV']
variable2name = {} 
variable2name['U'] = 'velocity'
variable2name['C'] = 'temperature'
variable2name['Q2'] = 'turbulent kinetic energy'
variable2name['Q2L'] = 'turbulent kinetic energy times a length scale'
variable2name['rho'] = 'density'
variable2name['L'] = 'length scale'
variable2name['nu_t'] = r'turbulent viscosity ($\nu_t$)'
variable2name['Kz'] = r'turbulent diffusivity ($\kappa_z$)'
variable2name['Kq'] = r'turbulent diffusivity ($\kappa_q$)'
variable2name['N_BV'] = r'Brunt-Väisälä frequency ($N_{BV}$)'
variable2name['algae'] = r'Algae concentration'
variable2name['biomass'] = r'Total algae biomass'
variable2name['net_growth'] = r'Net Algae Growth (loss + growth)'

variable2name['algae1'] = r'Algae concentration Diatoms'
variable2name['biomass1'] = r'Total algae biomass (species 1)'
variable2name['net_growth1'] = r'Net Algae Growth (loss + growth)'
variable2name['algae2'] = r'Algae concentration HABs'
variable2name['biomass2'] = r'Total algae biomass (species 2)'
variable2name['net_growth2'] = r'Net Algae Growth (loss + growth)'

variable2units = {}
variable2units['U'] = 'm/s'
variable2units['C'] = 'deg C'
variable2units['Q2'] = r'tke'
variable2units['Q2L'] = r'tke*L'
variable2units['rho'] = r'kg/m$^3$'
variable2units['L'] = r'm'
variable2units['nu_t'] = r'm$^2$/s'
variable2units['Kz'] = r'm$^2$/s'
variable2units['Kq'] = r'm$^2$/s'
variable2units['N_BV'] = r'1/s'
variable2units['algae'] = r'10$^6$ cells/m$^2$'
variable2units['algae1'] = r'10$^6$ cells/m$^2$'
variable2units['algae2'] = r'10$^6$ cells/m$^2$'
variable2units['biomass'] = r'10$^6$ cells'
variable2units['net_growth'] = r'hour$^{-1}$'



class WCRun:
    def __init__(self, N, z, save_output=True):
        self.N = N
        self.z = z
        self.first_data = True
        self.save_output = save_output
    def save_run_info(self, **kwargs):
        # Add to attributes
        self.dataset = self.dataset.assign_attrs(**kwargs)

    def read_from_file(self, filename):
        self.dataset = xr.open_dataset(filename)
        # df = self.dataset.to_dataframe()
        # df = df.dropna(how="all")
        # print(df)
        # assert(False)
        self.z = self.dataset.z.values

    def save_1d_data(self, time, **kwargs):
        ''' 1D data inclues variables like biomass, forcings & photic depth
            that varies only with time (and not with depth)'''
        if self.save_output:
            df = pd.DataFrame(index=[time]) 
            df.index.name = 'time'
            for key, value in kwargs.items():
                df[key] = value
            dfnc = df.to_xarray() 
            self.dataset = xr.concat([dfnc, self.dataset], dim='time')
        else: 
            return 
    
    def save_2d_data(self, time, **kwargs):
        if self.save_output:
            time_index = [time] * self.N
            index = pd.MultiIndex.from_arrays([self.z, time_index], names=["z", "time"])
            df = pd.DataFrame(index=index) 
            for key, value in kwargs.items():
                df[key] = value
            dfnc = df.to_xarray()
            # dfnc = (df.unstack('z').to_xarray().to_array('time'))

            if self.first_data: 
                self.dataset= dfnc 
                self.first_data = False
            else:
                self.dataset = xr.concat([dfnc, self.dataset], dim='time')
        else:
            return

    def save_dataset(self, filename):
        if self.save_output:
            self.dataset.to_netcdf(filename)

    ### PLOTTING UTILITIES ### 
    def seconds2hours(self, seconds):
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

        ax.plot(self.seconds2hours(ds.time), ds.values, '-o', markersize = 5, label = label)
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
        ds = self.dataset[variable]
        for i, ind in reversed(list(enumerate(plot_index))):
            time = ds.time.values[ind]
            if time == 0:
                 ax.plot(ds[variable].values[:,ind], 
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
                            label='t = %2.1f hr' % self.seconds2hours(time))
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
    


    def plot_phasing(self, ListOfSpecies, passed_string='', skip=1, show=True):
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
        axs[0].legend()

        # Add forcings
        time = self.dataset.time
        hours = self.seconds2hours(time)
        I_in = self.dataset.attrs['I_in']
        Px0 = self.dataset.attrs['Px0']
        T_Px = self.dataset.attrs['T_Px']
        diurnal = [diurnal_light(t, I_in, True) for t in time]
        ax2 = axs[2].twinx()
        axs[2].plot(hours[0:-2], diurnal[0:-2], label='Diurnal light', linewidth = 3, color='yellow')
        axs[2].plot([],[], 'o', label='Tidal forcing', alpha = 0.4, linewidth=3, color='skyblue')
        if self.dataset.attrs['T_Px'] != 0:
            pressure = [Px0*math.cos(2*math.pi*t/(3600*T_Px)) for t in time]
            ax2.plot(hours[0:-2], pressure[0:-2], '-o',  alpha = 0.4, label='Tidal forcing', color='skyblue')
        axs[2].legend()
        axs[2].grid(alpha = 0.5)
        axs[2].set_title("Temporal forcings")

        # Add concentrations 
        axs[1] = self.add_profile_to_axis( 'algae1', plot_index, legend_ind, axs[1])
        axs[1].set_title("Concentration of %s" % ListOfSpecies[0].name)

        axs[3] = self.add_profile_to_axis( 'algae2', plot_index, legend_ind, axs[3])
        axs[3].set_title("Concentration of %s" % ListOfSpecies[1].name)
        axs[3].set_xlim(0, 25)
        axs[1].set_xlim(0, 25)

        return fig, axs

def diurnal_light(t, I_max, diurnal):
    '''Function to generate estimate of light according to diurnal cycle.
     Inputs: t (in seconds); I_max (maximum light intensity/ light at noon)'''
    if diurnal:
        hour = t/3600
        period = (2*np.pi)/24
        phase_shift = 12 
        light = I_max * np.cos(period * (hour - phase_shift))
        light = light.clip(min=0) # During night, light is zero
        return light
    else:
        return I_max