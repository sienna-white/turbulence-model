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


# Filename for your netCDF output
out_fn = "demo_for_lisa.nc"

#********************** SPATIAL DOMAIN  ***************************
N = 80      # number of grid points
H = 10      # depth (meters)
dz = H/N    # grid spacing - may need to adjust to reduce oscillations
dt = 10   # (seconds) size of time step 
M  = 600 # number of time steps 

# Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
isave = 20

#********************** FIXED CONSTANTS  ***************************
rhoA = 1.23                 # DENSITY OF AIR, kg / m^3
rhoW = 1000                 # Density of Water
specific_heat_water = 4181  # J/kg-degC
specific_heat_air = 1007    # J/kg-degCxrh
c_d = 0.05                  # Drag coefficient for surface

#************ TEMPERATURE INITIAL CONDITION & FORCING ***************************
# Initialize thermocline based on tanh curve 
base_temp = 22
dtemp = 1.5 
stretch = 0.25 
STRATIFIED_INIT_TEMP = False 


# Heat flux
# flux_max = 0.005    # ignore for now unless trying to implement heat flux .. 
# top_temp = 33       # ignore unless trying to implement time-varying temperature
# bottom_temp = 30    # ignore unless trying to implement time-varying temperature

#********************** DEFINE HYDRODYNAMIC FORCINGS ***************************
# (1) PRESSURE 
Px0 = 2e-4         # Barotropic gradient forcing --> TIDES! 
T_Px = 1 #12           # Period [hours] on pressure gradient forcing. Set to 0 for steady

#************ WIND FORCING ***************************
WIND = 8            # constant wind speed, m/s
# bottom_speed = 0  # ignore unless trying to implement time-varying wind
# top_speed=3.5     # ignore unless trying to implement time-varying wind

#************ PHASING ***************************
# leave as-is unless you want to experiment here. only relevant if 
# you are implementing time-varying forcings.
TIDAL_PHASE_SHIFT = 0 
TEMP_PHASE_SHIFT = 0
WIND_PHASE_SHIFT = 0 
LIGHT_PHASE_SHIFT = 0

#********************** DEFINE ALGAL FORCINGS ***************************
# (1) Light 
# DIURNAL_LIGHT = False #     
background_turbidity =  0.6         # belive this is 1/m
I_in = 350                          # irradiance in
   


#***************************************************************************
#   Algae parameters
#***************************************************************************
'''
Define each type of algae as a class... 
Diatoms ws = -1.38e-5 m/s
Cyanobacteria = 1.38e-4 m/s
'''
initial_phytoplankton = 1e-2          

Algae1 = Algae_Species(k = 0.7,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                pmax = 0.1,         # maximum specific growth rate [1/hour]
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
# Shouldn't need to change this section EXCEPT for modifying initial condition for temperature? 
# ---------------------------------------------------------------------------

model = lib.WCModel(N=N,
                    H=H, 
                    dt=dt,
                    N_time_steps=M,
                    base_temp = base_temp) 

# Initialize z vector --> bottom at z[0]; top at z[N-1] or z[top] 
z = model.z

# Create a vector of time steps 
Times = model.get_time_steps() 

#***************************************************************************
#   Initialize algae  
#***************************************************************************
Algae1.set_initial_concentration(N, init=initial_phytoplankton, opt='constant')
Algae1.set_vertical_grid(H, N, dz)

Algae2.set_initial_concentration(N, init=initial_phytoplankton, opt='constant')
Algae2.set_vertical_grid(H, N, dz)

algae1 = Algae1.c
algae2 = Algae2.c

SelfShade = SelfShading(z, N, background_turbidity, self_shading=True)

#***************************************************************************
#   Initialize temperature / strafication profile 
#***************************************************************************

#  * * * EDIT HERE IF YOU WANT TO CHANGE INITIAL TEMP PROFILE!!!! 
C = np.tanh(z * stretch)*dtemp + base_temp 

# Calculate density profile 
rho = model.calculate_rho(C)    # Single scalar, linear equation of state

# Calculate Brunt-Vaisala frequency profile
N_BV2 = model.calculate_brunt_vaisala(rho)

#***************************************************************************
#   Initialize velocity + turbulent parameters 
#***************************************************************************

# Initialize pressure forcing
model.set_pressure_parameters(Px0, T_Px)

# Initalize velocity
U = np.zeros(N) 

# Turbulence quantities initialized to "SMALL"; Lengthscale parabolic
Q2, Q2L, Q, L, Gh, nu_t, Kq, Kz = model.initialize_turbulent_functions(N_BV2)

# Initialize based on initial condition
# only relevant if you saved an initial condition previously
# Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, U = model.check_initial_condition(Px0)
##########################################################################################   
#***************** end 
#***************************************************************************

def wc_advance(Up, Cp, Q2p, Q2Lp, rhop, Lp, nu_tp, Kzp, Kqp, N_BV2p, Ap1, Ap2, time_index):

    '''
    Time-advancing algorithm. Steps a single timestep for c, rho, q2, q2l, l, kz, nu_t, kq
    All diffusion/viscous terms handled implicitly
    '''

    #********************** TIME VARYING FORCINGS ***************************
    timestep = Times[time_index] 

    # [1] Light for phytoplankton growth
    Light =  model.diurnal_light(time_index, 1000, phase_shift=LIGHT_PHASE_SHIFT, diurnal=True)

    #  [2] Wind speed
    #       option 1: constant wind speed
    Wind = WIND

    #       option 2: time-varying wind speed
    # Wind = model.wind_speed(time_index, bottom_speed, top_speed, phase_shift=WIND_PHASE_SHIFT)

    #  Calculates wind stress 
    wind = (c_d * Wind)**2 * rhoA 
    Wstress= wind * dt/(dz*rhoW) 
    #**************************************************************************


    #***************************************************************************
    #   ADVANCE HYDRODYNAMIC VARIABLES
    #***************************************************************************
    # Extracts pressure from our forcing time series 
    Px   = model.get_pressure_at_timestep(timestep, phase_shift=TIDAL_PHASE_SHIFT)

    # TKE 
    Qp = np.sqrt(Q2p)

    # Update shear velocity at bottom boundary. Note explicit dependence on C_D
    ustar = model.calculate_ustar(Up[0]) 

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
    # option (1): Use Thomas algorithm to solve for C
    # C = model.advance_scalar(Kzp, Cp, heat_flux)

    # option (2): prescribe temperature profile (let values evolve over the day)
    # C = Temp_Profile # model.analytical_temperature_profile(t=timestep, bottom_temp=30, top_temp=33)

    # option(3): let same as the inital condition, if so, all options here stay commented
    

    # #***************************************************************************
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
        # photic_depth = model.calculate_photic_depth(Light, light)
        model.save_2d_data(Times[time_index], **data2d)

        # 1D data -- commenting out for now.  
        # data1d = {'biomass1': sum(algae1), 'biomass2' : sum(algae2), "photic_depth": photic_depth}
        # model.save_1d_data(Times[time_index], **data1d)

    return [U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, algae1, algae2]


#***************************************************************************
#   LOOP THROUGH TIME ! 
#***************************************************************************
for m in range(1,M):
   
    if m%1000 == 0:
        print('Time step = %d' % m)

    # Advance the model by one timestep
    output = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, algae1, algae2, time_index=m) 

    # Unpack output
    U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV2, algae1, algae2 = output



attributes = {"Px0": Px0, "T_Px": T_Px, "I_in": I_in, "Diurnal" : True,
              "Wind": WIND, "pmax1": Algae1.pmax, "pmax2": Algae2.pmax, "ws1": Algae1.ws, "ws2": Algae2.ws, "background_turbidity": background_turbidity}
model.save_run_info(**attributes) 
model.save_dataset(out_fn)

print("Total time = %f" % (time.time()  - t1))

