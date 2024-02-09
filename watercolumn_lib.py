import math
import numpy as np
import matplotlib.pyplot as plt
from turbulence_model import Turbulence_Model as TM 
import matplotlib as mpl


['U', 'C', 'Q2', 'Q2L', 'rho', 'L', 'nu_t', 'Kz', 'Kq', 'N_BV']
variable2name = {} 
variable2name['U'] = 'velocity'
variable2name['C'] = 'temperature'
variable2name['Q2'] = 'turbulent kinetic energy'
variable2name['Q2L'] = 'turbulent kinetic energy times a length scale'
variable2name['rho'] = 'density'
variable2name['L'] = 'length scale'
variable2name['nu_t'] = r'turbulent viscosity ($\nu_t$)'
variable2name['Kz'] = r'turbulent diffusivity ($K_z$)'
variable2name['Kq'] = r'turbulent diffusivity ($K_q$)'
variable2name['N_BV'] = r'Brunt-Vaisala frequency ($N_{BV}$)'
variable2name['algae'] = r'Algae Concentration'

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
variable2units['algae'] = r'?'


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
            saved_profiles[variable] = np.zeros((self.N, self.n_profiles))
        saved_profiles['time'] = np.zeros((self.n_profiles))
        self.saved_profiles = saved_profiles

    def store_z(self, z):
        self.z = z
    def save_profile_at_timestep(self, profile_num, time, U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae):
        profile_num = profile_num//self.isave
        self.saved_profiles['U'][:,profile_num] = U
        self.saved_profiles['C'][:,profile_num] = C
        self.saved_profiles['Q2'][:,profile_num] = Q2
        self.saved_profiles['Q2L'][:,profile_num] = Q2L
        self.saved_profiles['rho'][:,profile_num] = rho
        self.saved_profiles['L'][:,profile_num] = L
        self.saved_profiles['nu_t'][:,profile_num] = nu_t
        self.saved_profiles['Kz'][:,profile_num] = Kz
        self.saved_profiles['Kq'][:,profile_num] = Kq
        self.saved_profiles['N_BV'][:,profile_num] = N_BV
        self.saved_profiles['algae'][:,profile_num] = algae
        self.saved_profiles['time'][profile_num] = time

    def plot_profiles(self, variable, skip=1):
        fig, ax = plt.figure(), plt.gca()
        plots = np.arange(0, self.n_profiles, skip)
        for i, ind in enumerate(plots):
            if self.saved_profiles['time'][ind] == 0:
                 ax.plot(self.saved_profiles[variable][:,ind], 
                    self.z, '--',
                    color = 'k',
                    linewidth = 2, 
                    label='Initial condition')
            else:
                ax.plot(self.saved_profiles[variable][:,ind], 
                        self.z, '-',
                        color = mpl.cm.viridis(i/len(plots)),
                        linewidth = 2.5, 
                        alpha = 0.6,
                        label='t = %d sec' % self.saved_profiles['time'][ind],)
        ax.legend()
        ax.grid(alpha = 0.5)
        ax.set_ylabel('Depth (m)')
        ax.set_xlabel('%s (%s)' % (variable2name[variable], variable2units[variable]))
        # ax.hlines(0, color = 'k', linestyle = '--')
        ax.set_title(variable2name[variable])
        plt.show()