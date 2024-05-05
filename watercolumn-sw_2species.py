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
from phytoplankton import Algae_Species
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
M  = 1440*18*2 # 400  # number of time steps 

read_from_input=False
plot=True
output=True


# Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
isave = 300 #200 #00

# Algae parameters 
background_turbidity =  0.16
I_in = 350 

'''
Diatoms ws = -1.38e-5 m/s
Cyanobacteria = 1.38e-4 m/s
'''
# Show --> ws=1e-7

Algae2 = Algae_Species(k = 0.0034, 
                    pmax = 0.05, #08,
                    ws = 1.4e-5,
                    Hi = 40,
                    Li = 0.004,
                    name = "HAB",
                    self_shading=True,
                    net=True)

Algae1 = Algae_Species(k = 0.07,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                pmax = 0.05,# 5, #0.1, #0.05,     # maximum specific growth rate [1/hour]
                ws = -1.4e-6, #-1.4e-6, #1e-5, #-1e-6, #-1e-9,#-1e-3, #-200,       # vertical velocity [m/s]
                Hi = 40,         # half-saturation of light-limited growth [mu mol photons * m^2/s]
                Li = 0.006,      # specific loss rate [1/hour]
                name = "Diatoms",
                self_shading=True,
                net=True)

init = 5 # 20 #25 #50 #200 

# Pressure Forcing -> Need to modify to allow for time variable Px.
Px0 = 2e-4 # 2e-6  # Magnitude on pressure gradient forcing
T_Px = 12 #12 # 12.0  # Period [hours] on pressure gradient forcing. Set to 0 for steady

W = 0.0
# RUN_INFO='pressure=%2.2e_pmax=%2.2e_W=%f' % (Px0, Algae1.pmax, W)
RUN_INFO='NOTIDESpx=%1.2e_a1ws=%1.2e_a1pmax=%1.2e_a2ws=1.2%e_a2pmax=%1.2e_W=%1.1f' % (Px0, Algae1.ws, Algae1.pmax, Algae2.ws, Algae2.pmax, W)


########################################################################################## 

Algae1.set_initial_concentration(N, init=init, opt='constant')
Algae1.save_total_mass()
Algae1.set_vertical_grid(H, N, dz)

Algae2.set_initial_concentration(N, init=init, opt='constant')
Algae2.save_total_mass()
Algae2.set_vertical_grid(H, N, dz)


algae1 = Algae1.c
algae2 = Algae2.c

##########################################################################################



# Initial conditions for temperature profile
delC   = 5       # Change in temperature at initial themocline [deg C]; set to zero for Unstratified Case
zdelC  = -5      # Position of initial thermocline
dzdelC = 4       # Thickness of initial thermocline 
alpha  = 0.0     # Thermal expansivity, set to zero for passive scalar case
base_temp = 15   # Temperature of water column [deg C]

# Physical parameters 
z0 = 0.01         # Bottom roughness [m]
zb = 10*z0        # Bottom height [m]
g  = 9.81         # Gravity [m/s^2]
C_D = 0.0025      # Friction coefficient 
SMALL = 1e-6      # Noise floor for turbulence quantities [m/s?]
kappa = 0.4       # Von Karman constant
nu = 1e-6         # Kinematic viscosity [m^2/s]
rho0 = 1000       # Water density [kg/m^3]

##########################################################################################
# Mellor-Yamada closure parameters. All dimensionless --> no need to change  
A1=0.92 # [-]
A2=0.74 # [-]
B1=16.6 # [-]
B2=10.1 # [-]
C1=0.08 # [-]
E1=1.8  # [-]
E2=1.33 # [-]
E3=0.25 # [-]
Sq=0.2  # [-]

# Create shorthand beta for use in discretization 
beta = (dt/dz**2)
top = N-1

# Create a vector of time steps 
t = np.zeros((M))
t[1:M] = dt * (np.arange(1,M) - 1)
##########################################################################################

#***************************************************************************
#   Define supporting functions
#***************************************************************************
def calculate_sm(gh):
    num = B1**(-1/3) - A1*A2*gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))
    dem = (1-3*A2*gh*(B2+6*A1))*(1-9*A1*A2*gh)
    Sm  = num/dem
    return Sm

def calculate_sh(gh):
    return A2*(1-6*A1/B1)/(1-3*A2*gh*(B2+6*A1))

def calculate_brunt_vaisala(rho_, N_BV):
    for i in range(0,top):
        dpdz = (rho_[i+1] - rho_[i])/dz     # Density gradient 
        N_BV[i] = math.sqrt(abs((-g/rho0)* dpdz))
    N_BV[top] = math.sqrt(abs((-g/rho0)*(rho[top] - rho[top-1])/(dz)))
    return N_BV
#***************************************************************************

'''
Initialize all profiles and closure parameters 
Call once before time loop
Sets all forcing: pressure gradients, stresses, etc.
Should be used to adjust initial temperature/salinity profiles
Velocity initialized to zero
Turbulence quantities initialized to "SMALL"; Lengthscale parabolic
'''
# Initialize z vector --> bottom at z[0]; top at z[N-1] or z[top] .. or zzTop .. just kidding
z = [(-H + dz*(i + 0.5)) for i in range(N)]
z = np.array(z)

# Initalize arrays for temperature, density, Brunt-Vaisala frequency, velocity
empty_arrays = [np.zeros(N) for i in range(5)]
C, rho, N_BV, U, V = empty_arrays


#***************************************************************************
#   Initialize temperature / strafication profile 
#***************************************************************************
# Thermocline will be half above midpoint, half below 
half_height_thermocline = 0.5*dzdelC
C = base_temp + delC*(z - zdelC + 0.5*dzdelC)/dzdelC
#   For values of z BELOW the thermocline, set C = base_temp
C[(z <= (zdelC - half_height_thermocline))] = base_temp
#   For values of z ABOVE the thermocline, set C = base_temp + delC
C[(z > (zdelC + half_height_thermocline))] = base_temp + delC

# Calculate density profile 
rho = rho0*(1 - alpha*(C - base_temp))  # Single scalar, linear equation of state

# Calculate Brunt-Vaisala frequency profile
dpdz = (rho[1:top+1]-rho[0:top])/dz 
N_BV[0:top]  = np.sqrt(abs((-g/rho0)*dpdz))
N_BV[top] = np.sqrt(abs((-g/rho0)*(rho[top] - rho[top-1])/(dz)))

#***************************************************************************
#   Initialize velocity + turbulent parameters 
#***************************************************************************
Q2  = SMALL*np.ones(N)   # "seed" the turbulent field with small values, then let it evolve
Q2L = SMALL*np.ones(N)

# Initial length scale 
L = -kappa*H*(z/H)*(1-(z/H)) # Q2L(n,1)/Q2(n,1) = 1 at initialization

# Initialize empty arrays with size N
empty_arrays = [np.zeros(N) for i in range(6)]
Q, Sm, Sh, nu_t, Kq, Kz = empty_arrays

# Initial U 
U = U*0 # 0.1*(z + 2)
Q = np.sqrt(Q2)

# Initialize Gh (stratification correction)
Gh = -((N_BV*L)/(Q + SMALL))**2
Gh = np.clip(Gh, -0.28, 0.0233)

# Calculate Sm, Sh, nu_t, Kq, Kz
Sm = calculate_sm(Gh) 
Sh = calculate_sh(Gh) 
nu_t = (Sm * Q * L) + nu    # Turbulent diffusivity for Q2
Kq = (Sq * Q * L) + nu      # Turbulent viscosity
Kz = (Sh * Q * L) + nu

##########################################################################################

# Initialize based on initial condition
Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, U = lib.check_initial_condition(Px0)

##########################################################################################   

# Intialize an object for saving profiles throughout the model run
n_profiles = int(M/isave)
variables_to_save = ['U', 'C', 'Q2', 'Q2L', 'rho', 'L', 'nu_t', 'Kz', 'Kq', 'N_BV',
                    'algae1', 'algae2', 'biomass1', 'biomass2', 'net_growth1', 'net_growth2']
saved_profiles = lib.SavedProfiles(n_profiles, variables_to_save, N, isave)  
saved_profiles.Px0 = Px0
saved_profiles.T_Px = T_Px
saved_profiles.I_in = I_in


# Save initial condition (first profile at time zero)
data = {'U': U, 'C': C, 'Q2': Q2, 
        'Q2L': Q2L, 'rho': rho, 'L': L,
        'nu_t': nu_t, 'Kz': Kz, 'Kq': Kq,
        'N_BV': N_BV, 'algae1': Algae1.c, 'algae2' : Algae2.c, 'biomass1': sum(Algae1.c), 
        'biomass2':  sum(Algae2.c), 'net_growth1': Algae1.c*0, 'net_growth2' : Algae2.c*0}
saved_profiles.save_profile_at_timestep(0, 0, **data)

# Store z in our object so we can plot the profiles later 
saved_profiles.store_z(z)


def wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2, time):

    '''
    Time-advancing algorithm. Steps a single timestep for c, rho, q2, q2l, l, kz, nu_t, kq
    All diffusion/viscous terms handled implicitly
    '''
    #***************************************************************************
    #  Initialize Tridiagonal Arrays 
    #***************************************************************************
    # Initialize tridiagonal arrays for C/temperature. dC is the RHS vector
    aC, bC, cC, dC = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for turbulent kinetic energy
    aQ2, bQ2, cQ2, dQ2 = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for Q^2 * L (turbulent kinetic energy times a lengthscale)
    aQ2L, bQ2L, cQ2L, dQ2L = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for velocity
    aU, bU, cU, dU = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for algae1 
    aA1, bA1, cA1, dA1 = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for algae 
    aA2, bA2, cA2, dA2 = lib.initialize_abcd(N)

    Px = np.zeros(N)
    Q = np.sqrt(Q2)

    # Update pressure forcing term for the current timestep
    if T_Px == 0.0:
        Px = Px + Px0 # Steady and constant forcing for now
    else: 
        Px = Px + Px0*math.cos((2*math.pi*t[m] - 3600*3)/(3600*T_Px)) 

    # Update shear velocity at bottom boundary. Note explicit dependence on C_D
    ustar = abs(U[0])*math.sqrt(C_D); 

    #***************************************************************************
    #   Update stability parameters 
    #***************************************************************************
        
    # Update Gh (stratification correction)
    Gh = -((N_BV*L)/(Q + SMALL))**2
    Gh = np.clip(Gh, -0.28, 0.0233)

    # Calculate Sm, Sh, nu_t, Kq, Kz
    Sm = calculate_sm(Gh) 
    Sh = calculate_sh(Gh) 
    nu_t = (Sm * Q * L) + nu # Turbulent diffusivity for Q2
    Kq = (Sq * Q * L) + nu   # Turbulent viscosity
    Kz = (Sh * Q * L) + nu
    Kz = Kz.clip(SMALL,)     # Set floor on Kz so it's never zero 

    #***************************************************************************
    #   Store last time's step variables (f --> fp, q2 --> q2p, etc)
    #***************************************************************************
    Ap1 = algae1 
    Ap2 = algae2
    Cp = C
    Q2p,Q2Lp  = Q2, Q2L
    Lp, Kzp, Kqp, nu_tp = L, Kz, Kq, nu_t
    N_BVp = N_BV
    Up, Vp = U, V
    
    #***************************************************************************
    #   Advance velocity (U,V)
    #***************************************************************************

    # Moved function itself into watercolumn_lib.py
    aU, bU, cU, dU = lib.advance_velocity(Up, N, top, beta, nu_tp, Px, C_D, kappa, dt, W)
    
    # Use Thomas algorithm to solve for U
    U = lib.TDMA(aU, bU, cU, dU, N)

    #********************** LIGHT FOR THIS TIME STEP ***************************
    #***************************************************************************
    light =  lib.diurnal_light(t[m], 350, diurnal=True)
    background_turbidity = 0.0 #16
    light, photic_depth = self_shading([Algae1, Algae2], I_in=light, turbidity=background_turbidity, self_shading=True)
    #***************************************************************************
    
    #***************************************************************************
    #  [1]  Advance algae 1! 
    #***************************************************************************
    ws    = Algae1.ws 
    wsdtdz  = abs(ws*dt)/dz

    gamma1 = Algae1.get_loss_and_growth(I_in = light, current_concentration = Ap1)

    aA1, bA1, cA1, dA1 = lib.advance_algae(ws, wsdtdz, gamma1, beta, Kzp, Ap1, N, top, dt)

    # Thomas algorithm to solve for C
    algae1 = lib.TDMA(aA1, bA1, cA1, dA1, N)  
    Algae1.c = algae1

    #***************************************************************************
    #   ADVANCE SECOND ALGAL SPECIES
    #***************************************************************************
    ws    = Algae2.ws 
    wsdtdz  = abs(ws*dt)/dz

    gamma2 = Algae2.get_loss_and_growth(I_in = light, current_concentration = Ap2)
    
    aA2, bA2, cA2, dA2 = lib.advance_algae(ws, wsdtdz, gamma2, beta, Kzp, Ap2, N, top, dt)

    # Thomas algorithm to solve for C
    algae2 = lib.TDMA(aA2, bA2, cA2, dA2, N)  
    Algae2.c = algae2


    #***************************************************************************
    #   Advance scalars/density (C, rho) 
    #***************************************************************************
    aC[1:top] = -0.5*beta*(Kzp[1:top] + Kzp[0:top-1])
    bC[1:top] = 1 + 0.5*beta*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1])
    cC[1:top] = -0.5*beta*(Kzp[1:top] + Kzp[2:top+1])
    dC[1:top] = Cp[1:top]

    # Bottom-Boundary: no flux for scalars
    bC[0] = 1+0.5*beta*(Kzp[1] + Kzp[0])
    cC[0] = -0.5*beta*(Kzp[1] + Kzp[0])
    dC[0] =  Cp[0]

    # Top-Boundary: no flux for scalars
    aC[top] = -0.5*beta*(Kzp[top] + Kzp[top-1])
    bC[top] = 1+0.5*beta*(Kzp[top] + Kzp[top-1])
    dC[top] = Cp[-1]

    # Thomas algorithm to solve for C
    C = lib.TDMA(aC, bC, cC, dC, N)

    # Update density and Brunt-Vaisala frequency
    rho = rho0*(1-alpha*(C - 15))  
    N_BV = calculate_brunt_vaisala(rho, N_BV)

    #***************************************************************************
    #   Advance TKE / Q2  
    #***************************************************************************

    # Dissipation is a size (N-2) vector used for the non-boundary terms in the Q2 equation 
    diss = (2 * dt *(Q2p[1:top]**0.5))/(B1*Lp[1:top])
    aQ2[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[0:top-1])
    bQ2[1:top] = 1 + 0.5*beta*(Kqp[2:top+1] + 2*Kqp[1:top] + Kqp[0:top-1]) + diss 
    cQ2[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[2:top+1])
    dQ2[1:top] = Q2p[1:top] + 0.25*beta*nu_tp[1:top]*(Up[2:top+1]-Up[0:top-1])**2 - dt*Kzp[1:top]*(N_BVp[1:top]**2)

    # Bottom-Boundary Condition 
    Q2bot = B1**(2/3) * ustar**2
    bdryterm = 0.5*beta*Kqp[0]*Q2bot
    dissipation = 2 * dt *((Q2p[0]**0.5)/(B1*Lp[0]))
    bQ2[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + dissipation
    cQ2[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2[0] = Q2p[0] + dt*((ustar**4)/nu_tp[0]) - dt*Kzp[0]*(N_BVp[0]**2) + bdryterm

    # Top boundary condition
    dissipation =  2 * dt *((Q2p[top]**0.5)/(B1*Lp[top]))
    aQ2[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
    bQ2[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kq[top-1]) + dissipation
    dQ2[top] = Q2p[top] + 0.25*beta*nu_tp[top]*((Up[top] - Up[top-1])**2) -4*dt*Kzp[top]*(N_BVp[top]**2)

    # TDMA to solve for q2
    Q2 = lib.TDMA(aQ2, bQ2, cQ2, dQ2, N)

    # Prevent negative values from causing instabilities
    Q2[Q2 < 0] = SMALL

    #***************************************************************************
    #   Advance Q2 * L 
    #***************************************************************************
    diss = 2*dt*((Q2p[1:top]**0.5) / (B1*Lp[1:top]))*(1+E2*(Lp[1:top]/(kappa*abs(-H-z[1:top])))**2 \
                                                      + E3*(Lp[1:top]/(kappa*abs(z[1:top])))**2)

    aQ2L[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[0:top-1])
    bQ2L[1:top] = 1 + 0.5*beta*(Kqp[2:top+1] + 2*Kqp[1:top] + Kqp[0:top-1]) + diss
    cQ2L[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[2:top+1]) 
    dQ2L[1:top] = Q2Lp[1:top] + 0.25*beta*nu_tp[1:top]*E1*Lp[1:top] * (Up[2:top+1]-Up[0:top-1])**2 \
                                - 2*dt*Lp[1:top]*E1*Kzp[1:top]*(N_BVp[1:top]**2)

    # Bottom boundary Condition
    q2lbot = B1**(2/3) * (ustar**2) * kappa * zb
    bdryterm = 0.5*beta*Kqp[0]*q2lbot
    diss =  2 * dt *(Q2p[0]**0.5)/(B1*Lp[0])*(1+E2*(Lp[0]/(kappa*abs(-H-z[0])))**2 + E3*(Lp[0]/(kappa*abs(z[0])))**2)
    bQ2L[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + diss
    cQ2L[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2L[0] = Q2Lp[0] + dt*((ustar**4)/nu_tp[0])*E1*Lp[0] - dt*Lp[0]*E1*Kzp[0]*(N_BVp[0]**2) + bdryterm

    # Top boundary condition
    dissipation =  2 * dt *(Q2p[top]**0.5)/(B1*Lp[top])*(1+E2*(Lp[top]/(kappa*abs(-H-z[top])))**2 \
                                                  + E3*(Lp[top]/(kappa*abs(z[top])))**2)
    aQ2L[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
    bQ2L[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kqp[top-1]) + dissipation # Are we using kq or kqp here?
    dQ2L[top] = Q2Lp[top] + 0.25*beta*nu_tp[top]*E1*Lp[top]*(Up[top]-Up[top-1])**2 - 2*dt*Lp[-1]*E1*Kzp[top]*(N_BVp[top]**2)
    
    # TDMA to solve for q2
    Q2L = lib.TDMA(aQ2L, bQ2L, cQ2L, dQ2L, N)

    # Prevent negative values in Q2L 
    Q2L[Q2L < 0] = SMALL

    #  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    Q = np.sqrt(Q2)
    L = Q2L/(Q2 + SMALL)

    # Check length scale 
    ind = ((L**2)*(N_BV**2)) > (0.281*Q2) # Vectorized if-statement 
    if sum(ind) > 0: 
        Q2L[ind] = Q2[ind]*math.sqrt(0.281*Q2[ind]/(N_BV[ind]**2 + SMALL))
        L[ind] = Q2L[ind] / Q2[ind]
    L[abs(L) <= zb] = zb

    Kq = Sq*Q*L + nu
    nu_t = Sm*Q*L + nu
    Kz = Sh*Q*L + nu   

    return [U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2, gamma1, gamma2]


#***************************************************************************
#   LOOP THROUGH TIME ! 
#***************************************************************************
for m in range(1,M):
   
    # if m%200 == 0:
    #     print('Time step = %d' % m)

    # Advance the model by one timestep
    output = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2, m) 

    # Unpack output
    U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae1, algae2, gamma1, gamma2 = output

    if m%isave == 0:
        # Algae1.save_total_mass()
        # Algae2.save_total_mass()
        # Pack data into dictionary structure before saving 
        data = {'U': U, 'C': C, 'Q2': Q2, 
                'Q2L': Q2L, 'rho': rho, 'L': L,
                'nu_t': nu_t, 'Kz': Kz, 'Kq': Kq,
                'N_BV': N_BV, 'algae1': algae1, 'algae2': algae2,
                'biomass1': sum(algae1), 'biomass2' : sum(algae2),
                'net_growth1': gamma1, 'net_growth2' : gamma2}
        saved_profiles.save_profile_at_timestep(m, t[m], **data)

print(time.time()  - t1) 
# saved_profiles.output_final_to_csv("initial_condition-%s.csv" % RUN_INFO)

#***************************************************************************


f0, a0 = saved_profiles.plot_phasing([Algae1, Algae2], passed_string='', skip=3, show=True)
# f0.savefig('figures/two_species_wind/%s_phasing.png' % RUN_INFO)

f0, a0 = saved_profiles.plot_profiles('U', skip=4, passed_string='', show=True)
# f0.savefig('figures/two_species_wind/%s_U.png' % RUN_INFO)


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
