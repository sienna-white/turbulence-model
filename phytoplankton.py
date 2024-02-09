

import importlib  
import numpy as np 


N = 80    # number of grid points
H = 20    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
z = [(-H + dz*(i + 0.5)) for i in range(N)]
z = np.array(z)

class Algae_Species:
    def __init__(self, k=0.1, pmax = 1.0):
        self.k = 0.1
        self.pmax = 1.0 
        self.c = np.zeros(N) + 2
    def set_initial_concentration(self, N):
        self.c = np.zeros(N) + 2


# diatoms 
a1 = Algae_Species(k = 0.7e-6, pmax = 0.05)

# microcystis 
a2 = Algae_Species(k = 0.034e-6, pmax = 0.008)

turb = 0.6
I_in = 350 
# Initialize z vector --> bottom at z[0]; top at z[N-1] or z[top] .. or zzTop .. just kidding


nspecies = 2 


def self_shading(ListofAlgae, I_in, turbidity):
    sum1 = [(species.k * species.c) for species in ListofAlgae][0]
    corrected_z = H + z 
    I = np.zeros(N)
    vector = sum1 - turbidity*corrected_z
    print(vector)
    for i in range(N):
        I[i] = np.sum(vector[0:i])*corrected_z[i]

    I = I_in * np.exp(-I)   
    I[I>1000] = 1000

    print(I)
 
    # I = I_in * np.exp(np.sum(-sum1 - turbidity*corrected_z)
    # for zi in z:
    #     corrected_z = H - zi
    #     sum2 = sum1 - turbidity*zi 
    # i = I_in * np.exp(-sum2)
    # print(i)

self_shading([a1, a2], I_in=I_in, turbidity=turb)
    # I = I_in * np.exp(-k * C)
    # return I


# algae = np.zeros(N, nspecies)

# pi = p_max(i) * I/(Hi + I)

# I(z,t) = I_in * np.exp( k(i) * C())