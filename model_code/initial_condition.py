import numpy as np 
from constants import * 
import os 
import pandas as pd

def temp_profile(self, dtemp, stretch, STRATIFIED_INIT_TEMP):
    z = self.z
    if STRATIFIED_INIT_TEMP:
        print("Initializing stratified temperature profile...")
        centered_z = 2*z + self.H # center z vector around zero 
        return np.tanh(centered_z * stretch)*dtemp + self.base_temp
    else: 
        print("Initializing constant temperature profile...")
        return np.zeros(self.N) + self.base_temp

def initialize_arrays(self):
    Q2  = SMALL*np.ones(self.N)   # "seed" the turbulent field with small values, then let it evolve
    Q2L = SMALL*np.ones(self.N)
    Q = np.sqrt(Q2)
    z = self.z
    empty_arrays = [np.zeros(self.N) for i in range(5)]
    Sm, Sh, nu_t, Kq, Kz = empty_arrays
    return Q2, Q2L, L, Q, Sm, Sh, nu_t, Kq, Kz

def initialize_N_BV(self, rho):
    top = self.top 
    dpdz = (rho[1:top+1]-rho[0:top])/ self.dz 
    N_BV2 = np.zeros(self.N,)
    N_BV2[0:top]  = np.sqrt(abs((-g/rho0)*dpdz))
    N_BV2[top] = np.sqrt(abs((-g/rho0)*(rho[top] - rho[top-1])/(self.dz)))
    return N_BV2

def initialize_turbulent_functions(self, N_BV2):
    Q2  = SMALL*np.ones(self.N)   # "seed" the turbulent field with small values, then let it evolve
    Q2L = SMALL*np.ones(self.N)
    Q = np.sqrt(Q2)
    z = self.z
    L = -kappa*self.H*(z/self.H)*(1-(z/self.H)) # Q2L(n,1)/Q2(n,1) = 1 at initialization

    # Initialize Gh (stratification correction)
    Gh = -((N_BV2*L)/(Q + SMALL))**2
    Gh = np.clip(Gh, -0.28, 0.0233)
    nu_t, Kq, Kz  = self.calculate_turbulent_functions(Gh, Q, L)
    return Q2, Q2L,Q, L, Gh, nu_t, Kq, Kz


def check_initial_condition(self, Px0):
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
    N_BV2 = ic['N_BV'].values
    U = ic['U'].values 
    return Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, U
