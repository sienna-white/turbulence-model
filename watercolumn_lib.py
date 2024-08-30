import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd 
import os
from numba import jit 

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
    csv_name='initial_condition/initial_condition-pressure=%2.2e.csv' % Px0
    if os.path.isfile(csv_name): 
        pass
    else:
        csv_name='initial_condition/initial_condition-pressure=2.00e-07.csv' 
        print("No initial condition for this pressure gradient --> using default!")

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

def air_temperature(t, max_temp, diurnal):
    '''Function to generate estimate of heat flux according to diurnal cycle.
     Inputs: t (in seconds); flux_max (maximum heat intensity/ light at noon)'''
    if diurnal:
        hour = t/3600
        period = (2*np.pi)/24
        phase_shift = 12 
        temp = max_temp * np.cos(period * (hour - phase_shift))
        temp = temp.clip(min=18) # During night, light is zero
        return temp
    else:
        return temp
    

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



@jit 
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
        label1 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[0].name, ListOfSpecies[0].pmax, ListOfSpecies[0].ws)
        label2 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[1].name, ListOfSpecies[1].pmax, ListOfSpecies[1].ws)
        axs[0] = self.add_time_series_to_axis('biomass1', label1, axs[0])
        axs[0] = self.add_time_series_to_axis('biomass2', label2, axs[0])
        axs[0].legend()

        # Add forcings
        time = self.saved_profiles['time'][0:] 
        hours = self.seconds2hours(time)
        diurnal_light = [diurnal_light(t, self.I_in, True) for t in time]
        ax2 = axs[2].twinx()
        axs[2].plot(hours[0:-2], diurnal_light[0:-2], label='Diurnal light', linewidth = 3, color='yellow')
        axs[2].plot([],[], 'o', label='Tidal forcing', alpha = 0.4, linewidth=3, color='skyblue')
        if self.T_Px != 0:
            pressure = [self.Px0*math.cos(2*math.pi*t/(3600*self.T_Px)) for t in time]
            ax2.plot(hours[0:-2], pressure[0:-2], '-o',  alpha = 0.4, label='Tidal forcing', color='skyblue')
            print(hours, pressure)
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
        ax.set_xlabel(variable2units[variable])
        ax.set_title(variable2name[variable] + ' ' +  passed_string)
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
    


##########################################################################################
# Physical parameters 
z0 = 0.01         # Bottom roughness [m]
zb = 10*z0        # Bottom height [m]
g  = 9.81         # Gravity [m/s^2]
C_D = 0.0025      # Friction coefficient 
SMALL = 1e-6      # Noise floor for turbulence quantities [m/s?]
kappa = 0.4       # Von Karman constant
nu = 1e-6         # Kinematic viscosity [m^2/s]
rho0 = 1000       # Water density [kg/m^3]
alpha  =  2.1e-4      # Thermal expansivity, set to zero for passive scalar case

##########################################################################################
# Mellor-Yamada closure parameters. All dimensionless --> no need to change  
A1=0.92 # [-]
A2=0.74 # [-]
B1=16.6 # [-]
B2=10.1 # [-]
C1=0.08 # [-]
E1=1.8  # [-]
E2=1.33 # [-]
E3=0.25 # [-]
Sq=0.2  # [-]


class WCModel():
    def __init__(self, N, H, dt, N_time_steps, base_temp):
        self.N = N
        self.H = H
        self.dz = H/N
        self.dt = dt
        self.M = N_time_steps
        self.base_temp = base_temp
        print("Model run time set to %2.2f hours" % (self.M*dt/3600))
        
        # Create shorthand beta for use in discretization 
        self.beta = (dt/self.dz**2)
        self.top =  N-1
        self.z = self.get_z()
    
    def set_pressure_parameters(self, Px0, T_Px):
        self.Px0 = Px0
        self.T_Px = T_Px

    def get_time_steps(self):
        # Create a vector of time steps 
        t = np.zeros((self.M))
        t[1:self.M] = self.dt * (np.arange(1,self.M) - 1)
        return t
    
    def get_z(self):
        return np.array([(-self.H + self.dz*(i + 0.5)) for i in range(self.N)]) 

    def temp_profile(self, dtemp, stretch, STRATIFIED_INIT_TEMP):
        z = self.get_z() 
        if STRATIFIED_INIT_TEMP:
            print("Initializing stratified temperature profile...")
            centered_z = 2*z + self.H # center z vector around zero 
            return np.tanh(centered_z * stretch)*dtemp + self.base_temp
        else: 
            print("Initializing constant temperature profile...")
            return np.zeros(self.N) + self.base_temp
    
    def calculate_rho(self, C):
        rho = rho0*(1 - alpha*(C - self.base_temp))  # Single scalar, linear equation of state
        return rho 

    def initialize_arrays(self):
        Q2  = SMALL*np.ones(self.N)   # "seed" the turbulent field with small values, then let it evolve
        Q2L = SMALL*np.ones(self.N)
        Q = np.sqrt(Q2)
        z = self.get_z()
        empty_arrays = [np.zeros(self.N) for i in range(5)]
        Sm, Sh, nu_t, Kq, Kz = empty_arrays
        return Q2, Q2L, L, Q, Sm, Sh, nu_t, Kq, Kz
    
    def initialize_N_BV(self, rho):
        top = self.top 
        dpdz = (rho[1:top+1]-rho[0:top])/ self.dz 
        N_BV = np.zeros(self.N,)
        N_BV[0:top]  = np.sqrt(abs((-g/rho0)*dpdz))
        N_BV[top] = np.sqrt(abs((-g/rho0)*(rho[top] - rho[top-1])/(self.dz)))
        return N_BV

    def calculate_Gh(self, N_BV, L, Q):
        # Gh = -((N_BV*L)/(Q + SMALL))**2
        Gh = -(N_BV*L**2)/(Q + SMALL)**2
        Gh = np.clip(Gh, -0.28, 0.0233)
        return Gh 
    
    def initialize_turbulent_functions(self, N_BV):
        Q2  = SMALL*np.ones(self.N)   # "seed" the turbulent field with small values, then let it evolve
        Q2L = SMALL*np.ones(self.N)
        Q = np.sqrt(Q2)
        z = self.get_z()
        L = -kappa*self.H*(z/self.H)*(1-(z/self.H)) # Q2L(n,1)/Q2(n,1) = 1 at initialization

        # Initialize Gh (stratification correction)
        Gh = -((N_BV*L)/(Q + SMALL))**2
        Gh = np.clip(Gh, -0.28, 0.0233)
        nu_t, Kq, Kz  = self.calculate_turbulent_functions(Gh, Q, L)
        return Q2, Q2L,Q, L, Gh, nu_t, Kq, Kz

    def get_pressure_at_timestep(self, time):
        Px = np.zeros(self.N)
        # Update pressure forcing term for the current timestep
        if self.T_Px == 0.0:
            Px = Px + self.Px0 # Steady and constant forcing for now
        else: 
            Px = Px + self.Px0*math.cos((2*math.pi*time - 3600*3)/(3600*self.T_Px)) 
        return Px 

    def calculate_ustar(self, Up0):
        return abs(Up0)*math.sqrt(C_D)
    
    def calculate_turbulent_functions(self, Gh, Q, L):
        Sm = self.calculate_sm(Gh) 
        Sh = self.calculate_sh(Gh) 
        nu_t = (Sm * Q * L) + nu # Turbulent diffusivity for Q2
        Kq = (Sq * Q * L) + nu   # Turbulent viscosity
        Kz = (Sh * Q * L) + nu
        Kz = Kz.clip(SMALL,) 
        return nu_t, Kq, Kz
    
    
    def calculate_sm(self, gh):
        num = B1**(-1/3) - A1*A2*gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))
        dem = (1-3*A2*gh*(B2+6*A1))*(1-9*A1*A2*gh)
        Sm  = num/dem
        return Sm
    
    def calculate_sh(self, gh):
        return A2*(1-6*A1/B1)/(1-3*A2*gh*(B2+6*A1))

    def calculate_brunt_vaisala(self, rho):
        # when the density gradient is positive at night, things get 
        # funky b/c dp/dz should be reversing sign -- taking abs is causing 
        # our issue 
        N_BV = np.zeros(self.N)
        dpdz = np.zeros(self.N)
        for i in range(0,self.top):
            dpdz[i] = (rho[i+1] - rho[i])/self.dz     # Density gradient 
        N_BV = (-g/rho0)* dpdz
        N_BV[self.top] = (-g/rho0)*(rho[self.top] - rho[self.top-1])/(self.dz)
        # N_BV = np.sqrt(abs((-g/rho0)* dpdz))
        # N_BV[self.top] = np.sqrt(abs((-g/rho0)*(rho[self.top] - rho[self.top-1])/(self.dz)))
        return N_BV

    
    def calculate_photic_depth(self, Light, light_at_z):
        if Light<1:
            photic_depth = 0
        else:
            photic_depth = self.z[light_at_z>(0.1 * Light)][0]
        return photic_depth

    # **************************************************************************

    def advance_velocity(self, Up, nu_tp, Px, W=None):

        aU, bU, cU, dU = initialize_abcd(self.N) 
        top = self.top
        beta = self.beta 
        dt  = self.dt 

        aU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[0:top-1])
        bU[1:top] = 1 + beta/2*(nu_tp[2:top+1] + 2*nu_tp[1:top] + nu_tp[0:top-1])
        cU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[2:top+1])
        dU[1:top] = Up[1:top] - dt*Px[1:top]

        # Bottom boundary: log-law
        bU[0] = 1 + beta/2*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
        cU[0] = -beta/2*(nu_tp[1] + nu_tp[0])
        dU[0] = Up[0] - dt*Px[0]

        aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
        bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
        dU[top] = Up[top] - dt*Px[top] + W # beta*(nu_tp[top]/2)*W

        # # Top boundary: no stress
        # aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
        # bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
        # dU[top] = Up[top] - dt*Px[top]
        U = TDMA(aU, bU, cU, dU, self.N)
        return U


    def advance_algae(self, ws, gamma, Kzp, Ap):
        aA, bA, cA, dA = initialize_abcd(self.N)
        wsdtdz = abs(ws* self.dt)/self.dz
        top = self.top
        beta = self.beta 
        dt  = self.dt 

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

        A = TDMA(aA, bA, cA, dA, self.N)   
        return A


    def advance_scalar(self, Kzp, Cp, heat_flux):
        # Initialize tridiagonal arrays for C/temperature. dC is the RHS vector
        aC, bC, cC, dC = initialize_abcd(self.N)
        top = self.top
        beta = self.beta 
        dt  = self.dt 

        aC[1:top] = -0.5*beta*(Kzp[1:top] + Kzp[0:top-1])
        bC[1:top] = 1 + 0.5*beta*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1])
        cC[1:top] = -0.5*beta*(Kzp[1:top] + Kzp[2:top+1])
        dC[1:top] = Cp[1:top]

        # Bottom-Boundary: no flux for scalars
        bC[0] = 1 + beta/2*(Kzp[1] + Kzp[0])
        cC[0] = -beta/2*(Kzp[1] + Kzp[0])
        dC[0] =  Cp[0] 

        # Top-Boundary: no flux for scalars
        aC[top] = -0.5*beta*(Kzp[top] + Kzp[top-1])
        bC[top] = 1+0.5*beta*(Kzp[top] + Kzp[top-1])
        dC[top] = Cp[-1] + heat_flux # TEMP test
        
        # print('heat flux= ', heat_flux * self.dz )
        return aC, bC, cC, dC
    
    def advance_Q2(self, Q2p, Lp, Kqp, nu_tp, Up, Kzp, N_BVp, ustar):
        # Initialize tridiagonal arrays for turbulent kinetic energy
        aQ2, bQ2, cQ2, dQ2 = initialize_abcd(self.N)
        top = self.top
        beta = self.beta 
        dt  = self.dt 

        # Dissipation is a size (N-2) vector used for the non-boundary terms in the Q2 equation 
        diss = (2 * dt *(Q2p[1:top]**0.5))/(B1*Lp[1:top])
        aQ2[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[0:top-1])
        bQ2[1:top] = 1 + 0.5*beta*(Kqp[2:top+1] + 2*Kqp[1:top] + Kqp[0:top-1]) + diss 
        cQ2[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[2:top+1])                                  # buoyancy production term  (should be negative such that this term is adding TKE when density is unstable)
        dQ2[1:top] = Q2p[1:top] + 0.25*beta*nu_tp[1:top]*(Up[2:top+1]-Up[0:top-1])**2 - dt*Kzp[1:top]*(N_BVp[1:top])

        # Bottom-Boundary Condition 
        Q2bot = B1**(2/3) * ustar**2
        bdryterm = 0.5*beta*Kqp[0]*Q2bot
        dissipation = 2 * dt *((Q2p[0]**0.5)/(B1*Lp[0]))
        bQ2[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + dissipation
        cQ2[0] = -0.5*beta*(Kqp[1] + Kqp[0])
        dQ2[0] = Q2p[0] + dt*((ustar**4)/nu_tp[0]) - dt*Kzp[0]*(N_BVp[0]) + bdryterm

        # Top boundary condition
        dissipation =  2 * dt *((Q2p[top]**0.5)/(B1*Lp[top]))
        aQ2[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
        bQ2[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kqp[top-1]) + dissipation # sw note --> fixed typo w kqp
        dQ2[top] = Q2p[top] + 0.25*beta*nu_tp[top]*((Up[top] - Up[top-1])**2) -4*dt*Kzp[top]*(N_BVp[top])
        
        Q2  = TDMA(aQ2, bQ2, cQ2, dQ2, self.N) 

        return Q2

    def advance_Q2L(self, Q2p, Q2Lp, Lp, Kqp, nu_tp, Up, Kzp, N_BVp, ustar):

        # Initialize tridiagonal arrays for Q^2 * L (turbulent kinetic energy times a lengthscale)
        aQ2L, bQ2L, cQ2L, dQ2L = initialize_abcd(self.N)
        top = self.top
        beta = self.beta 
        dt  = self.dt 
        z = self.get_z()
        H = self.H

        diss = 2*dt*((Q2p[1:top]**0.5) / (B1*Lp[1:top]))*(1+E2*(Lp[1:top]/(kappa*abs(-H-z[1:top])))**2 \
                                                      + E3*(Lp[1:top]/(kappa*abs(z[1:top])))**2)

        aQ2L[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[0:top-1])
        bQ2L[1:top] = 1 + 0.5*beta*(Kqp[2:top+1] + 2*Kqp[1:top] + Kqp[0:top-1]) + diss
        cQ2L[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[2:top+1]) 
        dQ2L[1:top] = Q2Lp[1:top] + 0.25*beta*nu_tp[1:top]*E1*Lp[1:top] * (Up[2:top+1]-Up[0:top-1])**2 \
                                    - 2*dt*Lp[1:top]*E1*Kzp[1:top]*(N_BVp[1:top])

        # Bottom boundary Condition
        q2lbot = B1**(2/3) * (ustar**2) * kappa * zb
        bdryterm = 0.5*beta*Kqp[0]*q2lbot
        diss =  2 * dt *(Q2p[0]**0.5)/(B1*Lp[0])*(1+E2*(Lp[0]/(kappa*abs(-H-z[0])))**2 + E3*(Lp[0]/(kappa*abs(z[0])))**2)
        bQ2L[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + diss
        cQ2L[0] = -0.5*beta*(Kqp[1] + Kqp[0])
        dQ2L[0] = Q2Lp[0] + dt*((ustar**4)/nu_tp[0])*E1*Lp[0] - dt*Lp[0]*E1*Kzp[0]*(N_BVp[0]) + bdryterm

        # Top boundary condition
        dissipation =  2 * dt *(Q2p[top]**0.5)/(B1*Lp[top])*(1+E2*(Lp[top]/(kappa*abs(-H-z[top])))**2 \
                                                    + E3*(Lp[top]/(kappa*abs(z[top])))**2)
        aQ2L[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
        bQ2L[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kqp[top-1]) + dissipation # Are we using kq or kqp here?
        dQ2L[top] = Q2Lp[top] + 0.25*beta*nu_tp[top]*E1*Lp[top]*(Up[top]-Up[top-1])**2 - 2*dt*Lp[-1]*E1*Kzp[top]*(N_BVp[top])
        
        Q2L  = TDMA(aQ2L, bQ2L, cQ2L, dQ2L, self.N)
        return  Q2L 
    
    def add_noise_floor(self, vector):
        vector[vector < 0] = SMALL
        return vector 
    
    def calculate_lengthscale(self, Q2, Q2L, N_BV):
        L = Q2L/(Q2 + SMALL)
        # Check length scale 
        ind = ((L**2)*(N_BV**2)) > (0.281*Q2) # Vectorized if-statement 
        # Double check this if-statement doesn't get executed when 
        # NBV^2 is negative ! 
        if sum(ind) > 0: 
            Q2L[ind] = Q2[ind]*np.sqrt(0.281*Q2[ind]/(N_BV[ind] + SMALL))
            L[ind] = Q2L[ind] / Q2[ind]
        L[abs(L) <= zb] = zb
        return L 
    
#***************************************************************************
#   Define supporting functions













# C = base_temp + delC*(z - zdelC + 0.5*dzdelC)/dzdelC
# #   For values of z BELOW the thermocline, set C = base_temp
# C[(z <= (zdelC - half_height_thermocline))] = base_temp
# #   For values of z ABOVE the thermocline, set C = base_temp + delC
# C[(z > (zdelC + half_height_thermocline))] = base_temp + delC