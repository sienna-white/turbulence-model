import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd 
import os
import xarray as xr 

from constants import *

import advance_variables as advance 
import plotting_functions as plotting 
import calculate_physical_variables as cpv 
import save_output as save
import initial_condition as init 
import forcings as forcings
# ['U', 'C', 'Q2', 'Q2L', 'rho', 'L', 'nu_t', 'Kz', 'Kq', 'N_BV2']



def initialize_abcd(N):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    return a, b, c, d



def add_text(fig, ax, leg, species):
    text = "Algal Species : %s \n pmax = %2.2e \n ws = %2.2e" % (species.name, species.pmax, species.ws)
    # ax.text(0.6, 0.1, text, transform=ax.transAxes)
    ax.text(1.05, 0.05, text, transform=ax.transAxes)


def save_output(csv_file, var1, var2, output1, output2, header):

    # Check if the CSV file exists
    file_exists = os.path.isfile(csv_file)

    data = {header[0]: [var1], header[1]: [var2], header[2]: [output1], header[3]: [output2]}
    df = pd.DataFrame(data)

    # Append the DataFrame to the CSV file
    if file_exists:
        mode = 'a'
        header=False
    else:
        mode = 'w'
        header=True

    df.to_csv(csv_file, mode=mode, index=False, header=header)
    



class WCModel():


    ''' 
    ***************************************************************************
    (1) Initialize the model + set various parameters 
    ***************************************************************************
    '''   
    def __init__(self, N, H, dt, N_time_steps, base_temp, save_output=True):
        self.N = N
        self.H = H
        self.dz = H/N
        self.dt = dt
        self.M = N_time_steps
        self.base_temp = base_temp
        
        # Create shorthand beta for use in discretization 
        self.beta = (dt/self.dz**2)
        self.top =  N-1
        self.z = np.array([(-self.H + self.dz*(i + 0.5)) for i in range(self.N)])
        self.first_data = True
        self.save_output = save_output
        self.t = self.get_time_steps()
        print("Model run time set to %2.2f hours" % (self.M*dt/3600))

    def set_pressure_parameters(self, Px0, T_Px):
        print("Setting pressure parameters...")
        self.Px0 = Px0
        self.T_Px = T_Px

    def get_time_steps(self):
        # Create a vector of time steps 
        t = np.zeros((self.M))
        t[1:self.M] = self.dt * (np.arange(1,self.M) - 1)
        return t


    ''' 
    ***************************************************************************
    (2) Define utilities for saving data
    ***************************************************************************
    '''   
    # Save run information
    save_run_info = save.save_run_info

    # Save 1/2D data during the run
    save_1d_data = save.save_1d_data
    save_2d_data = save.save_2d_data

    # Export the dataset to a file 
    save_dataset = save.save_dataset
   
    ''' 
    ***************************************************************************
    (3) Setting initial condition
    ***************************************************************************
    '''    
    check_initial_condition = init.check_initial_condition
    temp_profile = init.temp_profile
    initialize_arrays = init.initialize_arrays
    initialize_N_BV = init.initialize_N_BV
    initialize_turbulent_functions = init.initialize_turbulent_functions
    ''' 
    ***************************************************************************
    (3) Calculate physical parameters
    ***************************************************************************
    '''    
    calculate_lengthscale = cpv.calculate_lengthscale
    calculate_rho = cpv.calculate_rho
    calculate_Gh = cpv.calculate_Gh
    calculate_ustar = cpv.calculate_ustar
    calculate_turbulent_functions = cpv.calculate_turbulent_functions
    calculate_sm = cpv.calculate_sm
    calculate_sh = cpv.calculate_sh
    calculate_brunt_vaisala = cpv.calculate_brunt_vaisala
    calculate_photic_depth = cpv.calculate_photic_depth
    add_noise_floor = cpv.add_noise_floor
    
    ''' 
    ***************************************************************************
    (3) Retrieve forcings @ timestep
    ***************************************************************************
    '''   
    analytical_temperature_profile = forcings.analytical_temperature_profile
    get_pressure_at_timestep = forcings.get_pressure_at_timestep
    read_forcings_from_file = forcings.read_forcings_from_file
    air_temperature = forcings.air_temperature
    temperature = forcings.temperature
    diurnal_light = forcings.diurnal_light
    wind = forcings.wind
    wind_speed = forcings.wind_speed


    ''' 
    ***************************************************************************
    (3) Calculate hydrodynamic + scalar equations
    ***************************************************************************
    '''       
    advance_velocity = advance.advance_velocity
    advance_algae = advance.advance_algae
    advance_scalar = advance.advance_scalar
    advance_Q2 = advance.advance_Q2
    advance_Q2L = advance.advance_Q2L


    ''' 
    ***************************************************************************
    (3) Plotting
    ***************************************************************************
    '''       
    plot_profiles = plotting.plot_profiles
    plot_shear = plotting.plot_shear
    get_plot_indices = plotting.get_plot_indices
    plot_time_series =  plotting.plot_time_series
    plot_phasing=  plotting.plot_phasing
    add_profile_to_axis = plotting.add_profile_to_axis
    add_time_series_to_axis = plotting.add_time_series_to_axis



#***************************************************************************
#   Define supporting functions













# C = base_temp + delC*(z - zdelC + 0.5*dzdelC)/dzdelC
# #   For values of z BELOW the thermocline, set C = base_temp
# C[(z <= (zdelC - half_height_thermocline))] = base_temp
# #   For values of z ABOVE the thermocline, set C = base_temp + delC
# C[(z > (zdelC + half_height_thermocline))] = base_temp + delC




# class SavedProfiles:
#     def __init__(self,n_profiles, variables_to_save, N, isave):
#         self.n_profiles = n_profiles
#         self.variables_to_save = variables_to_save
#         self.N = N
#         self.isave = isave
#         self.saved_profiles = self.initialize_saved_profiles()
#         self.initialize_saved_profiles()


#     def initialize_saved_profiles(self):
#         saved_profiles = {} 
#         for variable in self.variables_to_save:
#             saved_profiles[variable] = np.zeros((self.N, self.n_profiles+1))
#         saved_profiles['time'] = np.zeros((self.n_profiles+1))
#         self.saved_profiles = saved_profiles

#     def store_z(self, z):
#         self.z = z

#     def save_profile_at_timestep(self, profile_num, time, **kwargs):
#         #U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae, biomass):
#         profile_num = profile_num//self.isave
#         self.saved_profiles['time'][profile_num] = time
#         for key, value in kwargs.items():
#             self.saved_profiles[key][:,profile_num] = value
#         self.profile_num = profile_num

#     def output_final_to_csv(self, csv_name):
#         df = pd.DataFrame() 
#         for key in self.saved_profiles.keys():
#             print(key)
#             if key=="time":
#                 continue
#             df[key] = self.saved_profiles[key][:,self.profile_num]
#         print(csv_name)
#         df.to_csv(csv_name)
            

#     def seconds2hours(self, seconds):
#         return seconds/3600

#     def plot_phasing(self, ListOfSpecies, passed_string='', skip=1, show=True):

#     #     # first row: biomass 1, biomass 2
#     #     # second row : concentration 1, concentration 2 
#     #     # third row: diurnal light, tidal forcing 
#         fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize = (18, 10))
#         plot_index, legend_ind = self.get_plot_indices(skip)

#         axs = axs.flatten()

#         # Add biomass 
#         axs[0].set_title("Biomass over time")
#         label1 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[0].name, ListOfSpecies[0].pmax, ListOfSpecies[0].ws)
#         label2 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[1].name, ListOfSpecies[1].pmax, ListOfSpecies[1].ws)
#         axs[0] = self.add_time_series_to_axis('biomass1', label1, axs[0])
#         axs[0] = self.add_time_series_to_axis('biomass2', label2, axs[0])
#         axs[0].legend()

#         # Add forcings
#         time = self.saved_profiles['time'][0:] 
#         hours = self.seconds2hours(time)

#         diurnal_light = self.forcings["I_in"] #[diurnal_light(t, self.I_in, True) for t in time]
#         temperature = self.forings["T_in"] 
#         ax2 = axs[2].twinx()
#         axs[2].plot(hours[0:-2], diurnal_light[0:-2], label='Diurnal light', linewidth = 3, color='yellow')
#         axs[2].plot(hours[0:-2], temperature[0:-2], label='temp', linewidth = 3, color='crimson')

#         axs[2].plot([],[], 'o', label='Tidal forcing', alpha = 0.4, linewidth=3, color='skyblue')
#         if self.T_Px != 0:
#             pressure = [self.Px0*math.cos(2*math.pi*t/(3600*self.T_Px)) for t in time]
#             ax2.plot(hours[0:-2], pressure[0:-2], '-o',  alpha = 0.4, label='Tidal forcing', color='skyblue')
#             print(hours, pressure)
#         axs[2].legend()
#         axs[2].grid(alpha = 0.5)
#         axs[2].set_title("Temporal forcings")

#         # Add concentrations 
#         axs[1] = self.add_profile_to_axis( 'algae1', plot_index, legend_ind, axs[1])
#         axs[1].set_title("Concentration of %s" % ListOfSpecies[0].name)

#         axs[3] = self.add_profile_to_axis( 'algae2', plot_index, legend_ind, axs[3])
#         axs[3].set_title("Concentration of %s" % ListOfSpecies[1].name)
#         axs[3].set_xlim(0, 25)
#         axs[1].set_xlim(0, 25)

#         return fig, axs


#     def plot_concentration(self, ListOfSpecies, ListOfKeys, passed_string='', skip=1, show=True):
    
#         fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize = (15, 5))
#         plot_ind, legend_ind = self.get_plot_indices(skip)

#         def get_color(k, i):
#             if k==0:
#                 return mpl.cm.viridis(i/len(plot_ind))
#             if k==1:
#                 return mpl.cm.plasma(i/len(plot_ind))
            
#         ls = '-'

#         for k, key in enumerate(ListOfKeys):
#             self.add_profile_to_axis(key, plot_ind, legend_ind, axs[k])

#             leg = axs[k].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
#             axs[k].grid(alpha = 0.5)
#             axs[k].set_ylabel('Depth (m)')
#             axs[k].set_xlabel('%s (%s)' % (variable2name['algae'], variable2units['algae']))
#             axs[k].set_title("%s" %  ListOfSpecies[k].name + passed_string)
#             add_text(fig, axs[k], leg, ListOfSpecies[k])

#         plt.tight_layout()

#         if show:
#             plt.show()
#         return fig, axs


#     def add_time_series_to_axis(self, variable, label0, ax):
#         series = [self.saved_profiles[variable][0,i] for i in range(self.n_profiles)]
#         time = self.seconds2hours(self.saved_profiles['time'][0:len(series)])

#         label = label0 
#         ax.plot(time, series, '-o', markersize = 5, label = label)
#         ax.grid(alpha = 0.5)
#         ax.set_xlabel('Time (hr)')
#         return ax 


#     def plot_biomass(self, ListOfSpecies, ListOfKeys, passed_string='', show=True):

#         def seconds2hours(seconds):
#             return seconds/3600
        
#         fig, ax = plt.figure(figsize=(8,4)), plt.gca()

#         for i, key in enumerate(ListOfKeys):
#             label0 = "Biomass of %s" % ListOfSpecies[i].name 
#             ax = self.add_time_series_to_axis(key, label0, ax)

#         ax.set_ylabel('Algal biomass (%s)' % variable2units["biomass"])
#         ax.legend()
#         ax.set_title(variable2name["biomass"] + passed_string)
#         plt.tight_layout()
#         if show:
#             plt.show()
#         return fig, ax  
    

#     def does_biomass_increase(self):
#         return sum(self.saved_profiles['algae'][-1]) > sum(self.saved_profiles['algae'][0])
    

#     def add_profile_to_axis(self, variable, plot_index, legend_ind, ax):
#         ''' Add a vertical profile given by the key variable to a given axis'''
#         ls  = '-'
#         for i, ind in enumerate(plot_index):
#             if self.saved_profiles['time'][ind] == 0:
#                  ax.plot(self.saved_profiles[variable][:,ind], 
#                     self.z, '--',
#                     color = 'k',
#                     linewidth = 2, 
#                     label='Initial condition')
#             else:
#                 if i%legend_ind==0: 
#                     ax.plot(self.saved_profiles[variable][:,ind], 
#                             self.z, ls,
#                             color = mpl.cm.viridis(i/len(plot_index)),
#                             linewidth = 2.5, 
#                             alpha = 0.6,
#                             label='t = %2.1f hr' % self.seconds2hours(self.saved_profiles['time'][ind]))
#                 else:
#                     ax.plot(self.saved_profiles[variable][:,ind], 
#                             self.z, ls,
#                             color = mpl.cm.viridis(i/len(plot_index)),
#                             linewidth = 2.5, 
#                             alpha = 0.6)
        
#         ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
#         ax.grid(alpha = 0.5)
#         ax.set_ylabel('Depth (m)')
        
#         return ax 
    
#     def get_plot_indices(self, skip):
#         plots = np.arange(0, self.n_profiles, skip)
#         if len(plots)>10: 
#             legend_ind = np.floor(len(plots)/10)
#         else:
#             legend_ind = 1
#         return plots, legend_ind
    
#     def plot_profiles(self, variable, skip=1, passed_string='', show=True):
        
#         fig, ax = plt.figure(figsize=(8,4)), plt.gca()
#         plot_index, legend_ind = self.get_plot_indices(skip)
    
#         ls = '-'

#         self.add_profile_to_axis(variable, plot_index, legend_ind, ax)
#         ax.set_xlabel(variable2units[variable])
#         ax.set_title(variable2name[variable] + ' ' +  passed_string)
#         plt.tight_layout()
#         if show:
#             plt.show()
#         return fig, ax  
