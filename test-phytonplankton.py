import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import watercolumn_lib as lib

from phytoplankton import Algae_Species

# Spatial Parameters 
N = 100    # number of grid points
H = 20    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
dt = 60 #60   # (seconds) size of time step 
M = 400# 400  # number of time steps 

background_turbidity =  0.16
I_in = 350 

# Algae parameters 
diatoms = Algae_Species(k = 0.7,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                   pmax = 0.05,     # maximum specific growth rate [1/hour]
                   ws = 0,          # vertical velocity [cm/hour]
                   Hi = 40,         # half-saturation of light-limited growth [mu mol photons * m^2/s]
                   Li = 0.006,      # specific loss rate [1/hour]
                   name = "Diatoms")

diatoms.set_initial_concentration(N, init=5)
diatoms.set_vertical_grid(H, N, dz)


def diurnal_light(t, I_max):
    '''Function to generate estimate of light according to diurnal cycle.
     Inputs: t (in seconds); I_max (maximum light intensity/ light at noon)'''
    hour = t/3600
    period = (2*np.pi)/24
    phase_shift = 12 
    light = I_max * np.cos(period * (hour - phase_shift))
    light = light.clip(min=0) # During night, light is zero
    return light

# Set up plot
fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (14,5)) 

# Light 
ax[0].set_xlabel(r'Light Intensity ($\mu$mol photons/m$^2$/s)')

# Growth Rate
ax[1].set_xlabel('Growth Rate')

# Algae Concentration
ax[2].set_xlabel(r'Concentration of Algae')
   

total_time_steps = 4500
interval = 500 
for t in range(total_time_steps):
    past_concentration = diatoms.c
    I_in = diurnal_light(t*dt, I_max=350)
    I, photic_depth = diatoms.get_light_intensity(I_in)
    gamma = diatoms.get_loss_and_growth(I_in =I, current_concentration = past_concentration)
    diatoms.c = past_concentration + (gamma * past_concentration * dt) 

    if t % interval == 0:
        T = t/60 
        c = mpl.cm.cividis((t/interval)/(total_time_steps/interval))
        if I_in == 0:
            text = "t=%d hr (night)" % T 
        else:
            text = "t=%d hr (day)" % T 
        ax[2].plot(diatoms.c, diatoms.z, linewidth = 2, color = c, label = text)
        ax[1].plot(gamma * diatoms.c, diatoms.z, linewidth = 2, color = c, label = text)
        ax[0].plot(I, diatoms.z, linewidth = 2, color = c, label = text)

for axis in ax:
    axis.grid(alpha = 0.5)
    axis.legend(frameon=False)
    axis.set_ylim(-H,0)
ax[0].set_ylabel('Depth (m)')


# fig.savefig('Algae Concentration with diurnal light.png')
# fig.savefig('Algae Concentration with constant light.png')
plt.show()
