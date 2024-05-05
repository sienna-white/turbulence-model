import math
import numpy as np
import matplotlib.pyplot as plt
from turbulence_model import Turbulence_Model as TM 
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
variable2name['N_BV'] = r'Brunt-Vaisala frequency ($N_{BV}$)'
variable2name['algae'] = r'Algae concentration'
variable2name['biomass'] = r'Total algae biomass'
variable2name['net_growth'] = r'Net Algae Growth (loss + growth)'

variable2name['algae1'] = r'Algae concentration'
variable2name['biomass1'] = r'Total algae biomass'
variable2name['net_growth1'] = r'Net Algae Growth (loss + growth)'
variable2name['algae2'] = r'Algae concentration'
variable2name['biomass2'] = r'Total algae biomass'
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
variable2units['biomass'] = r'10$^6$ cells'
variable2units['net_growth'] = r'hour$^{-1}$'

def check_initial_condition(Px0):
    csv_name='./initial_condition/initial_condition-pressure=%2.2e.csv' % Px0
    if os.path.isfile(csv_name):
        print("Using initial condition from %s" % csv_name)
        ic = pd.read_csv(csv_name)
        Q2 = ic['Q2'].values
        Q2L = ic['Q2L'].values
        rho = ic['rho'].values
        L = ic['L'].values
        nu_t = ic['nu_t'].values
        Kz = ic['Kz'].values
        Kq = ic['Kq'].values
        N_BV = ic['N_BV'].values
        U = ic['U'].values 
        return Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, U
    else: 
        raise("No initial condition for this pressure gradient --> comment out this line!")

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

def initialize_abcd(N):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    return a, b, c, d

def TDMA(aX, bX, cX, dX, N):
    # Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
    # a = Lower Diag, b = Main Diag, c = Upper Diag, d = solution vector
    x = np.zeros(N)
    for i in range(1, N):
        bX[i] = bX[i] - aX[i]/bX[i-1]*cX[i-1]
        dX[i] = dX[i] - aX[i]/bX[i-1]*dX[i-1]
    x[-1] = dX[-1]/bX[-1]
    for i in range(N-2, -1, -1):
        x[i] = (1/bX[i])*(dX[i] - cX[i]*x[i+1])
    return x

class SavedProfiles:
    def __init__(self,n_profiles, variables_to_save, N, isave):
        self.n_profiles = n_profiles
        self.variables_to_save = variables_to_save
        self.N = N
        self.isave = isave
        self.saved_profiles = self.initialize_saved_profiles()
        self.initialize_saved_profiles()

    def initialize_saved_profiles(self):
        saved_profiles = {} 
        for variable in self.variables_to_save:
            saved_profiles[variable] = np.zeros((self.N, self.n_profiles+1))
        saved_profiles['time'] = np.zeros((self.n_profiles+1))
        self.saved_profiles = saved_profiles

    def store_z(self, z):
        self.z = z

    def save_profile_at_timestep(self, profile_num, time, **kwargs):
        #U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae, biomass):
        profile_num = profile_num//self.isave
        self.saved_profiles['time'][profile_num] = time
        for key, value in kwargs.items():
            self.saved_profiles[key][:,profile_num] = value
        self.profile_num = profile_num

    def output_final_to_csv(self, csv_name):
        df = pd.DataFrame() 
        for key in self.saved_profiles.keys():
            print(key)
            if key=="time":
                continue
            df[key] = self.saved_profiles[key][:,self.profile_num]
        print(csv_name)
        df.to_csv(csv_name)
            

    def seconds2hours(self, seconds):
        return seconds/3600

    def plot_phasing(self, ListOfSpecies, passed_string='', skip=1, show=True):

    #     # first row: biomass 1, biomass 2
    #     # second row : concentration 1, concentration 2 
    #     # third row: diurnal light, tidal forcing 
        fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize = (18, 10))
        plot_index, legend_ind = self.get_plot_indices(skip)

        axs = axs.flatten()

        # Add biomass 
        axs[0].set_title("Biomass over time")
        label1 = "Biomass of %s" % ListOfSpecies[0].name 
        label2 = "Biomass of %s" % ListOfSpecies[1].name
        axs[0] = self.add_time_series_to_axis('biomass1', label1, axs[0])
        axs[0] = self.add_time_series_to_axis('biomass2', label2, axs[0])
        axs[0].legend()

        # Add forcings
        time = self.saved_profiles['time'][0:] 
        hours = self.seconds2hours(time)
        diurnal = [diurnal_light(t, self.I_in, True) for t in time]
        ax2 = axs[1].twinx()
        axs[1].plot(hours[0:-2], diurnal[0:-2], label='Diurnal light', linewidth = 3, color='yellow')
        print(hours, diurnal)
        axs[1].plot([],[], 'o', label='Tidal forcing', alpha = 0.4, linewidth=3, color='skyblue')
        if self.T_Px != 0:
            pressure = [self.Px0*math.cos(2*math.pi*t/(3600*self.T_Px)) for t in time]
            ax2.plot(hours[0:-2], pressure[0:-2], '-o',  alpha = 0.4, label='Tidal forcing', color='skyblue')
            print(hours, pressure)
        axs[1].legend()
        axs[1].grid(alpha = 0.5)
        axs[1].set_title("Temporal forcings")

        # Add concentrations 
        axs[2] = self.add_profile_to_axis( 'algae1', plot_index, legend_ind, axs[2])
        axs[2].set_title("Concentration of %s" % ListOfSpecies[0].name)

        axs[3] = self.add_profile_to_axis( 'algae2', plot_index, legend_ind, axs[3])
        axs[3].set_title("Concentration of %s" % ListOfSpecies[1].name)
        axs[3].set_xlim(0, 25)
        axs[2].set_xlim(0, 25)

        return fig, axs


    def plot_concentration(self, ListOfSpecies, ListOfKeys, passed_string='', skip=1, show=True):
    
        fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize = (15, 5))
        plot_ind, legend_ind = self.get_plot_indices(skip)

        def get_color(k, i):
            if k==0:
                return mpl.cm.viridis(i/len(plot_ind))
            if k==1:
                return mpl.cm.plasma(i/len(plot_ind))
            
        ls = '-'

        for k, key in enumerate(ListOfKeys):
            self.add_profile_to_axis(key, plot_ind, legend_ind, axs[k])

            leg = axs[k].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            axs[k].grid(alpha = 0.5)
            axs[k].set_ylabel('Depth (m)')
            axs[k].set_xlabel('%s (%s)' % (variable2name['algae'], variable2units['algae']))
            axs[k].set_title("%s" %  ListOfSpecies[k].name + passed_string)
            add_text(fig, axs[k], leg, ListOfSpecies[k])

        plt.tight_layout()

        if show:
            plt.show()
        return fig, axs


    def add_time_series_to_axis(self, variable, label0, ax):
        series = [self.saved_profiles[variable][0,i] for i in range(self.n_profiles)]
        time = self.seconds2hours(self.saved_profiles['time'][0:len(series)])

        label = label0 
        ax.plot(time, series, '-o', markersize = 5, label = label)
        ax.grid(alpha = 0.5)
        ax.set_xlabel('Time (hr)')
        return ax 


    def plot_biomass(self, ListOfSpecies, ListOfKeys, passed_string='', show=True):

        def seconds2hours(seconds):
            return seconds/3600
        
        fig, ax = plt.figure(figsize=(8,4)), plt.gca()

        for i, key in enumerate(ListOfKeys):
            label0 = "Biomass of %s" % ListOfSpecies[i].name 
            ax = self.add_time_series_to_axis(key, label0, ax)

        ax.set_ylabel('Algal biomass (%s)' % variable2units["biomass"])
        ax.legend()
        ax.set_title(variable2name["biomass"] + passed_string)
        plt.tight_layout()
        if show:
            plt.show()
        return fig, ax  
    

    def does_biomass_increase(self):
        return sum(self.saved_profiles['algae'][-1]) > sum(self.saved_profiles['algae'][0])
    

    def add_profile_to_axis(self, variable, plot_index, legend_ind, ax):
        ''' Add a vertical profile given by the key variable to a given axis'''
        ls  = '-'
        for i, ind in enumerate(plot_index):
            if self.saved_profiles['time'][ind] == 0:
                 ax.plot(self.saved_profiles[variable][:,ind], 
                    self.z, '--',
                    color = 'k',
                    linewidth = 2, 
                    label='Initial condition')
            else:
                if i%legend_ind==0: 
                    ax.plot(self.saved_profiles[variable][:,ind], 
                            self.z, ls,
                            color = mpl.cm.viridis(i/len(plot_index)),
                            linewidth = 2.5, 
                            alpha = 0.6,
                            label='t = %2.1f hr' % self.seconds2hours(self.saved_profiles['time'][ind]))
                else:
                    ax.plot(self.saved_profiles[variable][:,ind], 
                            self.z, ls,
                            color = mpl.cm.viridis(i/len(plot_index)),
                            linewidth = 2.5, 
                            alpha = 0.6)
        
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        ax.grid(alpha = 0.5)
        ax.set_ylabel('Depth (m)')
        return ax 
    
    def get_plot_indices(self, skip):
        plots = np.arange(0, self.n_profiles, skip)
        if len(plots)>10: 
            legend_ind = np.floor(len(plots)/10)
        else:
            legend_ind = 1
        return plots, legend_ind
    
    def plot_profiles(self, variable, skip=1, passed_string='', show=True):
        
        fig, ax = plt.figure(figsize=(8,4)), plt.gca()
        plot_index, legend_ind = self.get_plot_indices(skip)
    
        ls = '-'

        self.add_profile_to_axis(variable, plot_index, legend_ind, ax)
        ax.set_title(variable2name[variable] + passed_string)
        plt.tight_layout()
        if show:
            plt.show()
        return fig, ax  


def add_text(fig, ax, leg, species):
    text = "Algal Species : %s \n pmax = %2.2e \n ws = %2.2e" % (species.name, species.pmax, species.ws)
    # ax.text(0.6, 0.1, text, transform=ax.transAxes)
    ax.text(1.05, 0.05, text, transform=ax.transAxes)



def save_output(csv_file, var1, var2, output, header):

    # Check if the CSV file exists
    file_exists = os.path.isfile(csv_file)

    data = {header[0]: [var1], header[1]: [var2], header[2]: [output]}
    df = pd.DataFrame(data)

    # Append the DataFrame to the CSV file
    if file_exists:
        mode = 'a'
        header=False
    else:
        mode = 'w'
        header=True

    df.to_csv(csv_file, mode=mode, index=False, header=header)
    

def advance_algae(ws, wsdtdz, gamma, beta, Kzp, Ap, N, top, dt):
    aA, bA, cA, dA = initialize_abcd(N)

    # If settling speed is UPWARD (swimming!)
    if ws>0:
        aA[1:top]  = -wsdtdz - beta/2 * (Kzp[0:top-1]+ Kzp[1:top]) 
        bA[1:top]  = 1 + wsdtdz - gamma[1:top]*dt  + beta/2*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1]) 
        cA[1:top]  = -beta/2 * (Kzp[1:top] + Kzp[2:top+1])
        dA = Ap

        # Bottom-Boundary: no flux for scalars
        bA[0] =  1 - (gamma[0]*dt) + beta/2*(Kzp[1] + Kzp[0]) + wsdtdz
        cA[0] = -beta/2 * (Kzp[1] + Kzp[0])
        dA[0] =  Ap[0]

        # Top-Boundary: no flux for scalars
        aA[top] = -wsdtdz - beta/2 * (Kzp[top] + Kzp[top-1])
        bA[top] = 1  - gamma[top]*dt + beta/2 * (Kzp[top] + Kzp[top-1]) # removed + wsdtdz 
        dA[top] = Ap[top]

    # If settling speed is DOWNWARD (sinking!)
    if ws<=0:
        aA[1:top]  = -beta/2 * (Kzp[0:top-1] + Kzp[1:top]) 
        bA[1:top]  = 1 + wsdtdz - gamma[1:top]*dt  + beta/2*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1]) 
        cA[1:top]  = -wsdtdz - beta/2 * (Kzp[1:top] + Kzp[2:top+1])
        dA = Ap

        # Bottom-Boundary: no flux for scalars
        bA[0] =  1 + wsdtdz - (gamma[0]*dt) + beta/2*(Kzp[1] + Kzp[0]) 
        cA[0] =  -wsdtdz -beta/2 * (Kzp[1] + Kzp[0])
        dA[0] =  Ap[0]

        # Top-Boundary: no flux for scalars
        aA[top] =  -beta/2 * (Kzp[top] + Kzp[top-1])
        bA[top] =  1 - gamma[top]*dt + beta/2 * (Kzp[top] + Kzp[top-1]) + wsdtdz # okay adding this here 
        dA[top] = Ap[top]

    return aA, bA, cA, dA




def advance_velocity(Up, N, top, beta, nu_tp, Px, C_D, kappa, dt, W=None):

    aU, bU, cU, dU = initialize_abcd(N)

    aU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[0:top-1])
    bU[1:top] = 1 + beta/2*(nu_tp[2:top+1] + 2*nu_tp[1:top] + nu_tp[0:top-1])
    cU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[2:top+1])
    dU[1:top] = Up[1:top] - dt*Px[1:top]

    # Bottom boundary: log-law
    bU[0] = 1 + beta/2*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
    cU[0] = -beta/2*(nu_tp[1] + nu_tp[0])
    dU[0] = Up[0] - dt*Px[0]

    # aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
    # bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
    # dU[top] = Up[top] - dt*Px[top] + beta*(nu_tp[top]/2)*W

    # Top boundary: no stress
    aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
    bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
    dU[top] = Up[top] - dt*Px[top]

    return aU, bU, cU, dU