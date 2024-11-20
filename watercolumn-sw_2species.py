'''
    Master code for 200B water column (vertical) cases
    Original MATLAB code from Lisa Lucas, modified by Tina Chow.
    Spring 2018: Mark Stacey and Michaella Chung.
    Adapted to Python by Alexandre Georges, Spring 2022.

###############################################################################
  Set initial conditions for u, q^2, and other turbulence quantities

# Q -- square root of turbulent kinetic energy [sqrt(Q^2)]
# sm, sh are stability parameters, calculated in the mellor-yamada closure. coefficient on turbulent mixing coefficients
# sm is for momentum , sh is for scalars -- 
# nu_t is turbulent viscosity 
# Kq is turbulent diffusion cofficient (for Q^2)
# Kz is turbulent diffusion coefficient for scalars
###########################################################################
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

sys.path.append("model_code/")
import watercolumn_lib as lib
from phytoplankton import Algae_Species, SelfShading
from phytoplankton import self_shading

import importlib
importlib.reload(lib) # Reload the module if it has changed


import sys, os, time 

t1 = time.time()  # Time our simluation 


'''
Because the indexing is a little confusing in Python vs. Matlab (0 is the bed, N-1 is 
the top of the water column), and then the point below the top is N-2, when indexing 
the top of the water column, I defined a variable called "top" to be N-1. This is just 
to make the code a little more readable, and hopefully less confusing. Hopefully the concept
of U[0] = bed velocity is a little more intuitive.

There's also a separate python file called watercolumn_lib.py that contains some functions
as well as a class definition called "SavedProfiles". This class is used to store the
profiles at each requested time step, and then plot them at the end. Mostly this was easier
than passing a bunch of arrays around between functions.
'''


#********************** SPATIAL DOMAIN  ***************************
N = 80    # number of grid points
H = 10    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
dt = 10 #60   # (seconds) size of time step 
ten_days = int(7*24*3600/dt)
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
Px0 = 2e-6 #2e-5  # gradient forcing
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
top_speed=3.5  

TIDAL_PHASE_SHIFT = 0 #float(os.environ["tidal_phase"])
TEMP_PHASE_SHIFT = 0
WIND_PHASE_SHIFT = 0 # float(os.environ["wind_phase"])
LIGHT_PHASE_SHIFT = TEMP_PHASE_SHIFT

#********************** DEFINE ALGAL FORCINGS ***************************

# (1) Light 
DIURNAL_LIGHT = False #  
background_turbidity =  0.6 #016
I_in = 350 
init = 1e-2

#********************** DEFINING OUTPUT ...***************************

out_fn = "TEST_heatflux_0T_1W_0D.nc"

#***************************************************************************
read_from_input=False
save_output=True
plot=True
output=False

#***************************************************************************
#   Algae parameters
#***************************************************************************

'''
Diatoms ws = -1.38e-5 m/s
Cyanobacteria = 1.38e-4 m/s
'''

Algae1 = Algae_Species(k = 0.7,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                pmax = 0.1,# 5,     # maximum specific growth rate [1/hour]
                ws = -1.4e-5,       # vertical velocity [m/s]
                Hi = 40,            # half-saturation of light-limited growth [mu mol photons * m^2/s]
                Li = 0.006,        # specific loss rate [1/hour]
                name = "Diatoms",
                self_shading=True,
                net=True)

Algae2 = Algae_Species(k = 0.034, 
                    pmax = 0.008, 
                    ws = 1.4e-5,
                    Hi = 40,
                    Li = 0.004,
                    name = "HAB",
                    self_shading=True,
                    net=True)



#***************************************************************************
#***************************************************************************
#   # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
isave = 2000 #200 #00

model = lib.WCModel(N=N,
                    H=H, 
                    dt=dt,
                    N_time_steps=M,
                    base_temp = base_temp) 

# model.read_forcings_from_file("forcing_data/Holt_CIMIS_data_processed_August2024_10s.csv")


#***************************************************************************

model.set_pressure_parameters(Px0, T_Px)

# RUN_INFO='pressure=%2.2e_pmax1=%2.2e_pmax2=%2.2e_Wind=%2.1e' % (Px0, Algae1.pmax, Algae2.pmax, WIND)
TITLE= "Px=%2.1e_tidalphasing=%d_windphasing=%d_Wind=%1.1f" % (Px0, TIDAL_PHASE_SHIFT, WIND_PHASE_SHIFT, top_speed) #Px=%2.1e_pmax1=%2.2e_pmax2=%2.2e" % (Px0, Algae1.pmax, Algae2.pmax)

# if T_Px>0:
#     RUN_INFO+="_TIDAL"
#     TITLE+=" (TIDAL)"
# if DIURNAL_LIGHT:
#     RUN_INFO+="_DIURNAL_LIGHT"
#     TITLE+= " (DIURNAL_LIGHT LIGHT)"

########################################################################################## 
output_csv = "temp.csv" # os.environ["output_csv"] # "PHASING_OUTPUT.csv"
def save_at_end(tidal_phasing, wind_phasing, diatom_biomass, hab_biomass):
    lib.save_output(output_csv, tidal_phasing, wind_phasing, diatom_biomass, hab_biomass, header=["tidal_phasing", "wind_phasing", "diatom_biomass", "hab_biomass"])

#***************************************************************************

'''
Initialize all profiles and closure parameters 
Call once before time loop
Sets all forcing: pressure gradients, stresses, etc.
Should be used to adjust initial temperature/salinity profiles
Velocity initialized to zero
Turbulence quantities initialized to "SMALL"; Lengthscale parabolic
'''
#***************************************************************************
#   Initialize arrays 
#***************************************************************************
# Initialize z vector --> bottom at z[0]; top at z[N-1] or z[top] 
z = model.z

# Create a vector of time steps 
Times = model.get_time_steps() 

#***************************************************************************
#   Initialize algae  
#***************************************************************************
Algae1.set_initial_concentration(N, init=init, opt='constant')
Algae1.set_vertical_grid(H, N, dz)

Algae2.set_initial_concentration(N, init=init, opt='constant')
Algae2.set_vertical_grid(H, N, dz)

algae1 = Algae1.c
algae2 = Algae2.c

SelfShade = SelfShading(z, N, background_turbidity, self_shading=True)
#***************************************************************************
#   Initialize temperature / strafication profile 
#***************************************************************************
C = model.temp_profile(dtemp=dtemp, stretch=stretch, STRATIFIED_INIT_TEMP=False)

# Calculate density profile 
rho = model.calculate_rho(C) # Single scalar, linear equation of state

# Calculate Brunt-Vaisala frequency profile
N_BV2 = model.calculate_brunt_vaisala(rho)

#***************************************************************************
#   Initialize velocity + turbulent parameters 
#***************************************************************************

# Initalize velocity
U = np.zeros(N) 

Q2, Q2L, Q, L, Gh, nu_t, Kq, Kz = model.initialize_turbulent_functions(N_BV2)

# Initialize based on initial condition
Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, U = model.check_initial_condition(Px0)
##########################################################################################   

def wc_advance(Up, Cp, Q2p, Q2Lp, rhop, Lp, nu_tp, Kzp, Kqp, N_BV2p, Ap1, Ap2, time_index):

    '''
    Time-advancing algorithm. Steps a single timestep for c, rho, q2, q2l, l, kz, nu_t, kq
    All diffusion/viscous terms handled implicitly
    '''

    #********************** TIME VARYING FORCINGS ***************************
    timestep = Times[time_index] 

    # Light # keep track of what light we're using ..
    Light =  model.diurnal_light(time_index, 1000, phase_shift=LIGHT_PHASE_SHIFT, diurnal=True)

    # Temperature profile + Wind speed 
    Temp = model.temperature(time_index, bottom_temp=bottom_temp, top_temp=top_temp, phase_shift =TEMP_PHASE_SHIFT)
    Temp_Profile = model.analytical_temperature_profile(time_index, bottom_temp, Temp)
    Wind = model.wind_speed(time_index, bottom_speed, top_speed, phase_shift=WIND_PHASE_SHIFT)
    wind = (c_d * Wind)**2 * rhoA 
    Px   = model.get_pressure_at_timestep(timestep, phase_shift=TIDAL_PHASE_SHIFT)

    # air_temp = model.air_temperature(time_index, max_temp=37, diurnal=True)
    # wind = model.wind(time_index) 
    

    # heat_in_air = air_temp * specific_heat_air * rhoA
    # heat_in_water = specific_heat_water * rhoW * Cp[-1]

    # if Cp[-1]>air_temp:
    #     heat_flux = -heat_in_air/(specific_heat_water * rhoW)
    # else:
    #     heat_flux = heat_in_air/(specific_heat_water * rhoW)


    rho0 = 1000
    Wstress= wind * dt/(dz*rho0) 
 
    Qp = np.sqrt(Q2p)

    # Update shear velocity at bottom boundary. Note explicit dependence on C_D
    ustar = model.calculate_ustar(Up[0]) 

    #***************************************************************************
    #   ADVANCE HYDRODYNAMIC VARIABLES
    #***************************************************************************
    #   Advance velocity (U,V)
    U = model.advance_velocity(Up, nu_tp, Px, W=Wstress)

    #   Advance TKE / Q2 
    Q2 = model.advance_Q2(Q2p, Lp, Kqp, nu_tp, Up, Kzp, N_BV2p, ustar) 

    #   Advance Q2 * L     
    Q2L = model.advance_Q2L(Q2p, Q2Lp, Lp, Kqp, nu_tp, Up, Kzp, N_BV2p, ustar)

    #***************************************************************************
    #   ADVANCE ALAGE 
    #***************************************************************************
    light = SelfShade.calc_self_shading([Algae1, Algae2], I_in=Light)
    
    #*********** Advance algae 1!  *********************************************
    gamma1 = Algae1.get_loss_and_growth(I_in = light, current_concentration = Ap1)
    algae1 = model.advance_algae(Algae1.ws, gamma1, Kzp, Ap1)
    Algae1.c = algae1

    #*********** Advance algae 2!  *********************************************
    gamma2 = Algae2.get_loss_and_growth(I_in = light, current_concentration = Ap2)
    algae2 = model.advance_algae(Algae2.ws, gamma2, Kzp, Ap2)
    Algae2.c = algae2

    #***************************************************************************
    #   Advance temperature
    #**************************************************************************
    # Thomas algorithm to solve for C
    # C = model.advance_scalar(Kzp, Cp, heat_flux)
    C = Temp_Profile # model.analytical_temperature_profile(t=timestep, bottom_temp=30, top_temp=33)
    #***************************************************************************
    # Update density and Brunt-Vaisala frequency
    rho = model.calculate_rho(C) 
    N_BV2 = model.calculate_brunt_vaisala(rho)

    # Prevent negative values from causing instabilities
    Q2 = model.add_noise_floor(Q2) 
    Q2L = model.add_noise_floor(Q2L) 

    #  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    Q = np.sqrt(Q2)

    # Get length scale by dividing Q2L by Q2 and checking for stability 
    L = model.calculate_lengthscale(Q2, Q2L, N_BV2)

    # Calculate turbulent diffusivities 
    nu_t, Kq, Kz = model.calculate_turbulent_functions(Gh, Q, L)

    if (time_index%isave) == 0:
        # Pack data into dictionary structure before saving 
        data2d = {'U': U, 'C': C, 'Q2': Q2, 
                'rho': rho, 'nu_t': nu_t, 'Kz': Kz, 'Kq': Kq,
                'N_BV2': N_BV2, 'algae1': algae1, 'algae2': algae2,
                'net_growth1': gamma1, 'net_growth2' : gamma2}
        photic_depth = model.calculate_photic_depth(Light, light)
        data1d = {'biomass1': sum(algae1), 'biomass2' : sum(algae2), "photic_depth": photic_depth}
        model.save_2d_data(Times[time_index], **data2d)
        model.save_1d_data(Times[time_index], **data1d)

    return [U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, algae1, algae2]


#***************************************************************************
#   LOOP THROUGH TIME ! 
#***************************************************************************
for m in range(1,M):
   
    # if m%1000 == 0:
    #     print('Time step = %d' % m)

    # Advance the model by one timestep
    output = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, algae1, algae2, time_index=m) 

    # Unpack output
    U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, algae1, algae2 = output


output=False

attributes = {"Px0": Px0, "T_Px": T_Px, "I_in": I_in, "Diurnal" : int(DIURNAL_LIGHT),
              "Wind": np.nan, "pmax1": Algae1.pmax, "pmax2": Algae2.pmax, "ws1": Algae1.ws, "ws2": Algae2.ws, "background_turbidity": background_turbidity}
model.save_run_info(**attributes) 
model.save_dataset(out_fn)

# if output:
#     change = diatoms.total_mass[-1] / diatoms.total_mass[0] 
#     depth_av_kz = np.mean(Kz)
#     print("Depth averaged U = ")

#               tidal_phasing, wind_phasing,  diatom_biomass, hab_biomass
save_at_end(TIDAL_PHASE_SHIFT, WIND_PHASE_SHIFT, sum(algae1), sum(algae2))

print("Total time = %f" % (time.time()  - t1))

# f0, a0 = model.plot_shear(skip=3,passed_string=TITLE, show=True)
# f0.savefig('figures/phase/%s_Shear.png' % TITLE)

# f0, a0 = model.plot_profiles('C', skip=3, passed_string='', show=True)
# f0.savefig('figures/phase/%s_Temp.png' % TITLE)

# print("Saving...", TITLE)
# f0, a0 = model.plot_profiles('U', skip=3, passed_string='', show=True)
# f0.savefig('figures/phase/%s_U.png' % TITLE)

# f0, a0 = model.plot_profiles('Kz', skip=3, passed_string='', show=True)
# f0.savefig('figures/phase/%s_KZ.png' % TITLE)

# f0, a0 = model.plot_profiles('C', skip=1, passed_string='', show=True)
# f0.savefig('figures/phase/%s_C.png' % TITLE)

# f0, a0 = model.plot_profiles('N_BV2', skip=3, passed_string='', show=True)
# f0.savefig('figures/phase/%s_N_BV2.png' % TITLE)

# f0, a0 = model.plot_profiles('algae1', skip=3, passed_string=TITLE, show=True)
# f0.savefig('figures/phase/%s_algae1.png' % TITLE)

# f0, a0 = model.plot_profiles('algae2', skip=3, passed_string=TITLE, show=True)
# f0.savefig('figures/phase/%s_algae2.png' % TITLE)

time = model.get_time_steps()
light = model.diurnal_light(time, 450, phase_shift=TEMP_PHASE_SHIFT, diurnal=True)
wind = model.wind_speed(time, bottom_speed, top_speed, phase_shift=WIND_PHASE_SHIFT)
time = Times/3600
period = (2*np.pi)/model.T_Px
# Update pressure forcing term for the current timestep
if model.T_Px == 0.0:
    Px = 0 + Px0 # Steady and constant forcing for now
else: 
    Px =  Px0*np.cos(period * (time - TIDAL_PHASE_SHIFT)) 

f0, a0 = model.plot_phasing([Algae1, Algae2], pressure=Px, light=light, wind=wind, passed_string=TITLE, skip=3, show=True)

plt.suptitle('%s' % TITLE)
f0.savefig('figures/phase/%s_PHASING.png' % TITLE)

#***************************************************************************


# f0.savefig('figures/two_species_wind1/%s_phasing.png' % RUN_INFO)



# f0, a0 = saved_profiles.plot_profiles('Kz', skip=5, passed_string=' ', show=True)
# f0.savefig('figures/two_species_wind1/%s_Kz.png' % RUN_INFO)

# f0, a0 = saved_profiles.plot_profiles('C', skip=5, passed_string=' ', show=False)
# f0.savefig('figures/two_species_wind/%s_Kz.png' % RUN_INFO)


# f0, a0 = saved_profiles.plot_biomass([Algae1, Algae2], ['biomass1', 'biomass2'], passed_string='', show=True)
# f0.savefig('figures/two_species/biomass-%s.png' % RUN_INFO)

# f0, a0 = saved_profiles.plot_concentration([Algae1, Algae2], ['algae1', 'algae2'], passed_string='', skip=4, show=True)
# f0.savefig('figures/two_species/concentration-%s.png' % RUN_INFO)

# f0, a0 = saved_profiles.plot_concentration([Algae1, Algae2], ['net_growth1', 'net_growth2'], passed_string='', skip=4, show=True)





# f0, a0 = saved_profiles.plot_profiles('net_growth2', skip=4, passed_string=RUN_INFO, show=True)

# f0, a0 = saved_profiles.plot_profiles('Kz', skip=4, passed_string=RUN_INFO, show=False)
# # f0.savefig('figures/two_species/Kz-%s.png' % RUN_INFO)

# f0, a0 = saved_profiles.plot_profiles('U', skip=4, passed_string=RUN_INFO, show=False)
# # f0.savefig('figures/two_species/U-%s.png' % RUN_INFO)

# # f0.savefig('figures/two_species/biomass-%s.png' % RUN_INFO)

# f0, a0 = saved_profiles.plot_profiles('algae1', skip=4, passed_string=RUN_INFO, show=False)
# f0.savefig('figures/two_species/algae-%s.png' % RUN_INFO)





# if output:
#     result= saved_profiles.does_biomass_increase()
#     save_at_end(result)

    

# # saved_profiles.plot_profiles('C', skip=2)

# # saved_profiles.plot_profiles('L', skip=2)
# f0, a0 = saved_profiles.plot_profiles('Kz', skip=2, passed_string=RUN_INFO)
# f0.savefig('output/kappa_z-%s.png' % RUN_INFO)
# print(diatoms.total_mass)
# plt.plot( diatoms.total_mass)
# plt.show()
