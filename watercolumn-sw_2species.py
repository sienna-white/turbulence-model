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

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import watercolumn_lib as lib
import watercolumn_run_lib as wrl
from phytoplankton import Algae_Species, SelfShading
from phytoplankton import self_shading
import time
import sys
import os 
import pandas as pd 

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

# Spatial Parameters 
N = 80    # number of grid points
H = 10    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
dt = 10 #60   # (seconds) size of time step 
seventy_two_hrs=259200
M  = seventy_two_hrs# 1440*18*2 # 400  # number of time steps 

# Initialize thermocline based on tanh curve 
base_temp = 15 
dtemp = 1.5 
centered_z = 2*z + H # center z vector around zero 
stretch = 0.25 

model = lib.WCModel(N=N,
                    H=H, 
                    dt=dt,
                    N_time_steps=M,
                    base_temp = 15) 


read_from_input=False
save_output=True
plot=True
output=False

# output_csv=os.getenv("output_csv")
# ws = float(os.getenv("ws"))
# pmax = float(os.getenv("pmax"))
# pressure = float(os.getenv("pressure"))

# Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
isave = 600 #200 #00


#***************************************************************************
#   Algae parameters
#***************************************************************************
background_turbidity =  0.0016
I_in = 350 
DIURNAL = True #  

'''
Diatoms ws = -1.38e-5 m/s
Cyanobacteria = 1.38e-4 m/s
'''

Algae1 = Algae_Species(k = 0.07,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                pmax = 0.07,# 5,     # maximum specific growth rate [1/hour]
                ws = -1.4e-5,       # vertical velocity [m/s]
                Hi = 40,            # half-saturation of light-limited growth [mu mol photons * m^2/s]
                Li = 0.0001,        # specific loss rate [1/hour]
                name = "Diatoms",
                self_shading=True,
                net=True)

Algae2 = Algae_Species(k = 0.0034, 
                    pmax = 0.01, 
                    ws = 1.4e-5,
                    Hi = 40,
                    Li = 0.001,
                    name = "HAB",
                    self_shading=True,
                    net=True)


init = 1
#***************************************************************************

# Pressure Forcing -> Need to modify to allow for time variable Px.
Px0 = 2e-6 # 2e-6  # Magnitude on pressure gradient forcing
T_Px = 0#12 #12 # 12.0  # Period [hours] on pressure gradient forcing. Set to 0 for steady

rhoA = 1.23  # kg / m^3
Wind = 2 
# u_star =  # m/s >> 0.05 is  drag coefficient, 10 is my wind speed 
WIND = (0.05 * Wind)**2 * rhoA  # this is rho * u*^2
RUN_INFO='pressure=%2.2e_pmax1=%2.2e_pmax2=%2.2e_Wind=%2.1e' % (Px0, Algae1.pmax, Algae2.pmax, WIND)
TITLE= "Px=%2.1e, Wind=%2.1e STRAT OFF" % (Px0,  WIND)

if T_Px>0:
    RUN_INFO+="_TIDAL"
    TITLE+=" (TIDAL)"
if DIURNAL:
    RUN_INFO+="_DIURNAL"
    TITLE+= " (DIURNAL LIGHT)"

########################################################################################## 

# def save_at_end(depth_av_kz, ws, result):
#     print("Depth averaged turbulent dissipation = %f" % depth_av_kz)
#     lib.save_output(output_csv, depth_av_kz, ws, result, header=["depth_averaged_kz", "ws", "output"])


##########################################################################################
# Create shorthand beta for use in discretization 
beta = (dt/dz**2)
top = N-1

# Create a vector of time steps 
time = model.get_time_steps() 

# np.zeros((M))
# t[1:M] = dt * (np.arange(1,M) - 1)
##########################################################################################


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
# Initialize z vector --> bottom at z[0]; top at z[N-1] or z[top] .. or zzTop .. just kidding
z = model.get_z() # np.array([(-H + dz*(i + 0.5)) for i in range(N)]) 

# Initalize arrays for temperature, density, Brunt-Vaisala frequency, velocity
empty_arrays = [np.zeros(N) for i in range(5)]
C, rho, N_BV, U, V = empty_arrays

# Intialize an object for saving profiles throughout the model run
RUN_TEST= wrl.WCRun(N, z, save_output=save_output) 

#***************************************************************************
#   Initialize algae  
#***************************************************************************
Algae1.set_initial_concentration(N, init=init, opt='constant')
Algae1.save_total_mass()
Algae1.set_vertical_grid(H, N, dz)

Algae2.set_initial_concentration(N, init=init, opt='constant')
Algae2.save_total_mass()
Algae2.set_vertical_grid(H, N, dz)

algae1 = Algae1.c
algae2 = Algae2.c

SelfShade = SelfShading(z, N, background_turbidity, self_shading=True)
#***************************************************************************
#   Initialize temperature / strafication profile 
#***************************************************************************
C = model.temp_profile(dtemp=dtemp, stretch= stretch)

# Calculate density profile 
rho = model.calculate_rho(C) # Single scalar, linear equation of state

# Calculate Brunt-Vaisala frequency profile
N_BV = model.initialize_N_BV(N_BV, rho)

#***************************************************************************
#   Initialize velocity + turbulent parameters 
#***************************************************************************
# Q2, Q2L, L, Q, Sm, Sh, nu_t, Kq, Kz = model.initialize_arrays()

Q2, Q2L, Q, L, Gh, nu_t, Kq, Kz = model.initialize_turbulent_functions(N_BV)

##########################################################################################
# Initialize based on initial condition
Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, U = lib.check_initial_condition(Px0)
##########################################################################################   


def wc_advance(Up, Cp, Q2p, Q2Lp, rhop, Lp, nu_tp, Kzp, Kqp, N_BVp, Ap1, Ap2, time_index):

    '''
    Time-advancing algorithm. Steps a single timestep for c, rho, q2, q2l, l, kz, nu_t, kq
    All diffusion/viscous terms handled implicitly
    '''
    #********************** TIME VARYING FORCINGS ***************************
    Light =  lib.diurnal_light(time[time_index], 350, diurnal=DIURNAL)
    Px  = model.get_pressure_at_timestep(time[time_index]) 
    rho0 = 1000
    Wstress= WIND * dt/(dz*rho0) 
 
    Qp = np.sqrt(Q2p)

    # Update shear velocity at bottom boundary. Note explicit dependence on C_D
    ustar = model.calculate_ustar(Up[0]) 

    #***************************************************************************
    #   Update stability parameters 
    #***************************************************************************
    # Update Gh (stratification correction)
    Gh = model.calculate_Gh(N_BVp, Lp, Qp)

    # Update turbulent diffusivities
    nu_t, Kq, Kz = model.calculate_turbulent_functions(Gh, Qp, Lp)

    #***************************************************************************
    #   ADVANCE HYDRODYNAMIC VARIABLES
    #***************************************************************************
    #   Advance velocity (U,V)
    U = model.advance_velocity(Up, nu_tp, Px, W=None)

    #   Advance TKE / Q2 
    Q2 = model.advance_Q2(Q2p, Lp, Kqp, nu_tp, Up, Kzp, N_BVp, ustar) 

    #   Advance Q2 * L     
    Q2L = model.advance_Q2L(Q2p, Q2Lp, Lp, Kqp, nu_tp, Up, Kzp, N_BVp, ustar)

    #***************************************************************************
    #   ADVANCE ALAGE 
    #***************************************************************************
    light = SelfShade.calc_self_shading([Algae1, Algae2], I_in=Light)
    
    #  [1]  Advance algae 1!  ****************************************************
    ws    = Algae1.ws 
    wsdtdz  = abs(ws*dt)/dz

    gamma1 = Algae1.get_loss_and_growth(I_in = light, current_concentration = Ap1)

    algae1 = model.advance_algae(ws, wsdtdz, gamma1, Kzp, Ap1)

    Algae1.c = algae1

    #  [2]  Advance algae 2!  ****************************************************
    ws    = Algae2.ws 
    wsdtdz  = abs(ws*dt)/dz

    gamma2 = Algae2.get_loss_and_growth(I_in = light, current_concentration = Ap2)

    algae2 = model.advance_algae(ws, wsdtdz, gamma2, Kzp, Ap2)
    
    Algae2.c = algae2

    #***************************************************************************
    #   Advance scalars/density (C, rho) 
    #**************************************************************************
    aC, bC, cC, dC = model.advance_scalar(Kzp, Cp)

    # Thomas algorithm to solve for C
    C =  Cp # lib.TDMA(aC, bC, cC, dC, N)

    # Update density and Brunt-Vaisala frequency
    rho = model.calculate_rho(C) 
    N_BV = model.calculate_brunt_vaisala(rho, N_BV)

    # Prevent negative values from causing instabilities
    Q2 = model.add_noise_floor(Q2) 

    # Prevent negative values in Q2L 
    Q2L = model.add_noise_floor(Q2L) 

    #  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    Q = np.sqrt(Q2)

    # Get length scale by dividing Q2L by Q2 and checking for stability 
    L = model.calculate_lengthscale(Q2, Q2L, N_BV)

    # Calculate turbulent diffusivities 
    nu_t, Kq, Kz = model.calculate_turbulent_functions(Gh, np.sqrt(Q2), L)

    if (time%isave) == 0:
        # Pack data into dictionary structure before saving 
        data2d = {'U': U, 'C': C, 'Q2': Q2, 
                'rho': rho, 'nu_t': nu_t, 'Kz': Kz, 'Kq': Kq,
                'N_BV': N_BV, 'algae1': algae1, 'algae2': algae2,
                'net_growth1': gamma1, 'net_growth2' : gamma2}
        if Light<1:
            photic_depth = 0
        else:
            photic_depth = z[light>(0.1 * Light)][0]
        data1d = {'biomass1': sum(algae1), 'biomass2' : sum(algae2), "photic_depth": photic_depth}
        RUN_TEST.save_2d_data(time, **data2d)
        RUN_TEST.save_1d_data(time, **data1d)

    return [U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2]


#***************************************************************************
#   LOOP THROUGH TIME ! 
#***************************************************************************
for m in range(1,M):
   
    # if m%200 == 0:
    #     print('Time step = %d' % m)

    # Advance the model by one timestep
    output = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2, time_index=m) 

    # Unpack output
    U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2 = output


output=False

attributes = {"Px0": Px0, "T_Px": T_Px, "I_in": I_in, "Diurnal" : int(DIURNAL),
              "Wind": Wind, "pmax1": Algae1.pmax, "pmax2": Algae2.pmax, "ws1": Algae1.ws, "ws2": Algae2.ws, "background_turbidity": background_turbidity}
RUN_TEST.save_run_info(**attributes) 
RUN_TEST.save_dataset("stratified_model_Px=6_diurnal_2.nc")

if output:
    change = diatoms.total_mass[-1] / diatoms.total_mass[0] 
    depth_av_kz = np.mean(Kz)
    print("Depth averaged U = ")
    save_at_end(depth_av_kz, diatoms.ws, change)
    
print("Total time = %f" % (time.time()  - t1))

f0, a0 = RUN_TEST.plot_profiles('U', skip=1, passed_string=TITLE, show=True)
f0.savefig('figures/stratification/%s_U.png' % TITLE)

f0, a0 = RUN_TEST.plot_profiles('Kz', skip=1, passed_string=TITLE, show=True)
f0.savefig('figures/stratification/%s_KZ.png' % TITLE)

f0, a0 = RUN_TEST.plot_profiles('C', skip=1, passed_string=TITLE, show=True)
f0.savefig('figures/stratification/%s_C.png' % TITLE)

f0, a0 = RUN_TEST.plot_profiles('N_BV', skip=1, passed_string=TITLE, show=True)
f0.savefig('figures/stratification/%s_N_BV.png' % TITLE)

f0, a0 = RUN_TEST.plot_phasing([Algae1, Algae2], passed_string=TITLE, skip=3, show=True)
plt.suptitle('%s' % TITLE)
f0.savefig('figures/stratification/%s_PHASING.png' % TITLE)

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




assert(False)
f0, a0 = saved_profiles.plot_profiles('net_growth2', skip=4, passed_string=RUN_INFO, show=True)

f0, a0 = saved_profiles.plot_profiles('Kz', skip=4, passed_string=RUN_INFO, show=False)
# f0.savefig('figures/two_species/Kz-%s.png' % RUN_INFO)

f0, a0 = saved_profiles.plot_profiles('U', skip=4, passed_string=RUN_INFO, show=False)
# f0.savefig('figures/two_species/U-%s.png' % RUN_INFO)

# f0.savefig('figures/two_species/biomass-%s.png' % RUN_INFO)

f0, a0 = saved_profiles.plot_profiles('algae1', skip=4, passed_string=RUN_INFO, show=False)
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
