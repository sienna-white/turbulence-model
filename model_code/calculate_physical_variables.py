from constants import * 
import numpy as np
import math 

def calculate_rho(self, C):
    rho = rho0*(1 - alpha*(C - self.base_temp))  # Single scalar, linear equation of state
    return rho 

def calculate_Gh(self, N_BV2, L, Q):
    # Gh = -((N_BV2*L)/(Q + SMALL))**2
    Gh = -(N_BV2*L**2)/(Q + SMALL)**2
    Gh = np.clip(Gh, -0.28, 0.0233)
    return Gh 

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
    N_BV2 = np.zeros(self.N)
    dpdz = np.zeros(self.N)
    for i in range(0,self.top):
        dpdz[i] = (rho[i+1] - rho[i])/self.dz     # Density gradient 
    N_BV2 = (-g/rho0)* dpdz
    N_BV2[self.top] = (-g/rho0)*(rho[self.top] - rho[self.top-1])/(self.dz)
    # N_BV = np.sqrt(abs((-g/rho0)* dpdz))
    # N_BV[self.top] = np.sqrt(abs((-g/rho0)*(rho[self.top] - rho[self.top-1])/(self.dz)))
    return N_BV2


def calculate_photic_depth(self, Light, light_at_z):
    if Light<1:
        photic_depth = 0
    else:
        photic_depth = self.z[light_at_z>(0.1 * Light)][0]
    return photic_depth

def add_noise_floor(self, vector):
    vector[vector < 0] = SMALL
    return vector 

def calculate_lengthscale(self, Q2, Q2L, N_BV2):
    L = Q2L/(Q2 + SMALL)
    # Check length scale 
    ind = ((L**2)*(N_BV2**2)) > (0.281*Q2) # Vectorized if-statement 
    # Double check this if-statement doesn't get executed when 
    # NBV^2 is negative ! 
    if sum(ind) > 0: 
        Q2L[ind] = Q2[ind]*np.sqrt(0.281*Q2[ind]/(N_BV2[ind] + SMALL))
        L[ind] = Q2L[ind] / Q2[ind]
    L[abs(L) <= zb] = zb
    return L 