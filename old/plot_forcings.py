import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import cmocean as cmo 
sys.path.append("model_code/")
import watercolumn_lib as lib



#********************** SPATIAL DOMAIN  ***************************
N = 80    # number of grid points
H = 10    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
dt = 10 #60   # (seconds) size of time step 
ten_days = int(2*24*3600/dt)
M  = ten_days# seventy_two_hrs# 1440*18*2 # 400  # number of time steps 

#********************** FIXED CONSTANTS  ***************************
rhoA = 1.23  # DENSITY OF AIR, kg / m^3
rhoW = 1000  # Density of Water
specific_heat_water = 4181 # J/kg-degC
specific_heat_air = 1007 # J/kg-degCxrh
c_d = 0.05   # Drag coefficient 

#********************** INITIAL CONDITION ***************************
# Initialize thermocline based on tanh curve 
base_temp = 22
dtemp = 1.5 
stretch = 0.25 

#********************** DEFINE HYDRODYNAMIC FORCINGS ***************************
# (1) PRESSURE 
Px0 = 2e-6  # gradient forcing
T_Px = 12 #$12 #12 # 12.0  # Period [hours] on pressure gradient forcing. Set to 0 for steady

# (2) Wind
# Wind = 0                        # u_star =m/s >> 0.05 is  drag coefficient, 10 is my wind speed 
# WIND = (c_d * Wind)**2 * rhoA  # this is rho * u*^2

# (3) Heat flux
STRATIFIED_INIT_TEMP = False 
flux_max = 0.005

# Water temperature  temperature(t, bottom_temp, top_temp, phase_shift = 12)
top_temp = 33
bottom_temp = 30
bottom_speed = 0 
top_speed = 3.5  

TIDAL_PHASE_SHIFT = 0
LIGHT_PHASE_SHIFT = 0
TEMP_PHASE_SHIFT = 0 
WIND_PHASE_SHIFT = 0


model = lib.WCModel(N=N,
                    H=H, 
                    dt=dt,
                    N_time_steps=M,
                    base_temp = base_temp) 
model.set_pressure_parameters(Px0, T_Px)


# Create a vector of time steps 
Times = model.get_time_steps() 

# Light 
Light =  model.diurnal_light(Times, 450, phase_shift=LIGHT_PHASE_SHIFT, diurnal=True)

# Temperature profile + Wind speed 

Wind = model.wind_speed(Times, bottom_speed, top_speed, phase_shift=WIND_PHASE_SHIFT)
wind = (c_d * Wind)**2 * rhoA 

# Pressure
time = Times/3600
period = (2*np.pi)/model.T_Px
Px0 = 2e-6
# Update pressure forcing term for the current timestep
if model.T_Px == 0.0:
    Px = 0 + Px0 # Steady and constant forcing for now
else: 
    Px =  Px0*np.cos(period * (time - TIDAL_PHASE_SHIFT)) 

Temp = model.temperature(Times, bottom_temp=bottom_temp, top_temp=top_temp, phase_shift =TEMP_PHASE_SHIFT)

# Px   = model.get_pressure_at_timestep(Times, phase_shift=TIDAL_PHASE_SHIFT)

hours = Times/3600

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 6))
axs = axs.ravel() 

axs[0].plot(hours, Temp, linewidth=3, alpha=0.8,  color="#af0522")
ax0 = axs[0].twinx()
axs[0].tick_params(axis='y', colors='#af0522')
axs[0].set_ylabel("Water temperature [deg C]")
ax0.plot(hours, Light, linewidth=3, alpha=0.8, color="#f3b151")
ax0.tick_params(axis='y', colors='#f3b151')
ax0.set_ylabel("Light [PAR]")

axs[0].set_ylim(30, 35)
axs[0].set_title("Surface temperature/Light")

axs[1].plot(hours, Wind, linewidth=3, alpha=0.8, color="#aed8fb")
axs[1].set_title("Wind Speed")
axs[1].set_ylabel("Wind Speed [m/s]")

axs[2].plot(hours, Px, linewidth=3, alpha=0.8, color="#0e2d74")
axs[2].set_ylabel("Tidal Pressure [Pa]")
axs[2].hlines(0, 0, hours[-1], color="black", alpha=0.25)
axs[2].set_title("Pressure")

for ax in axs:
    ax.grid(alpha=0.25)
    ax.set_xlim(0,hours[-1])

axs[2].set_xlabel("Time [hours]")
plt.tight_layout()

fig.savefig("Forcings.png")



fig, ax = plt.figure(figsize=(8,4)), plt.gca()



for i in range(1, 52):
    ind = i*150
    color = cmo.cm.haline(i/52)
    Temp_Profile = model.analytical_temperature_profile(Times[ind], bottom_temp, Temp[ind])
    if i%4==0:
        ax.plot(Temp_Profile, model.z, linewidth=3, alpha=0.85, color=color, label="t = %d hr" % (Times[ind]/3600))
    else:
        ax.plot(Temp_Profile, model.z, linewidth=3, alpha=0.85, color=color)


ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
ax.grid(alpha = 0.5)
ax.set_ylabel('Depth (m)')
ax.set_xlabel("Temperature [deg C]")

plt.tight_layout()
fig.savefig("Temperature_profile.png") 
# for i, ind in reversed(list(enumerate(plot_index))):
#     time = ds.time.values[ind]
#     if time == 0:
#             ax.plot(ds.isel(time=0), 
#             ds.z, '--',
#             color = 'k',
#             linewidth = 2, 
#             label='Initial condition')
#     else:
#         if i%legend_ind==0: 
#             ax.plot(ds.isel(time=ind), 
#                     self.z, ls,
#                     color = mpl.cm.viridis(i/len(plot_index)),
#                     linewidth = 2.5, 
#                     alpha = 0.6,
#                     label='t = %2.1f hr' % seconds2hours(time))
#         else:
#             ax.plot(ds.isel(time=ind), 
#                     self.z, ls,
#                     color = mpl.cm.viridis(i/len(plot_index)),
#                     linewidth = 2.5, 
#                     alpha = 0.6)

