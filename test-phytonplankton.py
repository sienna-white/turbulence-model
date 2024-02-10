import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import watercolumn_lib as lib

from phytoplankton import Algae_Species

# Spatial Parameters 
N = 80    # number of grid points
H = 20    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
dt = 10 #60   # (seconds) size of time step 
M = 400# 400  # number of time steps 

# Algae parameters 
diatoms = Algae_Species(k = 0.7,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                   pmax = 0.05,     # maximum specific growth rate [1/hour]
                   ws = 0,          # vertical velocity [cm/hour]
                   Hi = 40,         # half-saturation of light-limited growth [mu mol photons * m^2/s]
                   Li = 0.006,      # specific loss rate [1/hour]
                   name = "Diatoms")

diatoms.set_initial_concentration(N, init=10)
diatoms.set_vertical_grid(H, N, dz)

background_turbidity =  0.16
I_in = 350 

fig2 = plt.figure()
ax2 = plt.gca()
ax2.set_xlabel(r'Growth Rate')
ax2.set_ylabel('Depth (m)')
# ax.hlines(photic_depth, 0, I_in, color = 'black', linestyle = '--', label = 'Photic Depth')
ax2.set_ylim(-H,0)
ax2.grid(alpha = 0.5)


fig = plt.figure()
ax1 = plt.gca()
ax1.set_xlabel(r'Concentration of Algae')
ax1.set_ylabel('Depth (m)')
# ax.hlines(photic_depth, 0, I_in, color = 'black', linestyle = '--', label = 'Photic Depth')
ax1.set_ylim(-H,0)
ax1.grid(alpha = 0.5)
   

for t in range(10):
    past_concentration = diatoms.c
    gamma = diatoms.get_loss_and_growth(I_in = 350, current_concentration = past_concentration)
    diatoms.c = past_concentration + gamma * past_concentration * dt 
    ax1.plot(diatoms.c, diatoms.z, linewidth = 2, label = 't=%d' % t)
    ax2.plot(gamma, diatoms.z, linewidth = 2, label = 't=%d' % t)

ax1.legend()
ax2.legend()

plt.show()
