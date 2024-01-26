'''
    Master code for 200B water column (vertical) cases
    Original MATLAB code from Lisa Lucas, modified by Tina Chow.
    Spring 2018: Mark Stacey and Michaella Chung.
    Adapted to Python by Alexandre Georges, Spring 2022.

'''

import math
import numpy as np
import matplotlib.pyplot as plt
from turbulence_model import Turbulence_Model as TM 
import matplotlib as mpl
import watercolumn_lib as lib



PLOT=True 
isave=100     # increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
savecount=1

# Spatial Parameters 
N=80    # number of grid points
H=20    # depth (meters)
dz=H/N  # grid spacing - may need to adjust to reduce oscillations
dt=1 #60   # (seconds) size of time step 
M=2000  # number of time steps 

# Physical parameters 
z0 = 0.01         # Bottom roughness [m]
zb = 10*z0        # Bottom height [m]
g  = 9.81         # Gravity [m/s^2]
C_D = 0.0025      # Friction coefficient 
SMALL = 1e-6      # Noise floor for turbulence quantities
kappa = 0.4       # Von Karman constant
nu = 1e-6         # Kinematic viscosity [m^2/s]
rho0 = 1000       # Water density [kg/m^3]

# Initial conditions for temperature profile
delC   = 5       # Change in temperature at initial themocline [deg C]; set to zero for Unstratified Case
zdelC  = -5      # Position of initial thermocline
dzdelC = 4       # Thickness of initial thermocline 
alpha  = 0.0     # Thermal expansivity, set to zero for passive scalar case
base_temp = 15   # Temperature of water column [deg C]

# Pressure Forcing -> Need to modify to allow for time variable Px.
Px0 = 0.001  # Magnitude on pressure gradient forcing
T_Px = 12.0 # Period [hours] on pressure gradient forcing. Set to 0 for steady

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

beta = (dt/dz**2)

top = N-1


'''
Initialize all profiles and closure parameters 
Call once before time loop
Sets all forcing: pressure gradients, stresses, etc.
Should be used to adjust initial temperature/salinity profiles
Velocity initialized to zero
Turbulence quantities initialized to "SMALL"; Lengthscale parabolic

'''

# Initialize z vector --> bottom at z=-H, free surface at 0
z = [(-H + dz*(i + 0.5)) for i in range(N)]
z = np.array(z)

# Initalize arrays for temperature, density, Brunt-Vaisala frequency, velocity
empty_arrays = [np.zeros(N) for i in range(5)]
C, rho, N_BV, U, V = empty_arrays

####################################
## Initialize temperature profile ## 
####################################
# Thermocline will be half above midpoint, half below 
half_height_thermocline = 0.5*dzdelC
C = base_temp + delC*(z - zdelC + 0.5*dzdelC)/dzdelC

#   For values of z BELOW the thermocline, set C = base_temp
C[(z <= (zdelC - half_height_thermocline))] = base_temp

#   For values of z ABOVE the thermocline, set C = base_temp + delC
C[(z > (zdelC + half_height_thermocline))] = base_temp + delC

# Calculate density profile 
rho = rho0*(1 - alpha*(C - base_temp))  # Single scalar, linear equation of state

# Caclulate Brunt-Vaisala Frequency 
# N_BV = np.sqrt(abs((-g/rho0)*dpdz)) 
# N_BV[0] = math.sqrt(abs((-g/rho0)*(rho[1]-rho[0])/(dz)))
# for i in range(0,N-1):
#   dpdz = (rho[i+1] - rho[i])/dz     # Density gradient 
#   N_BV[i] = math.sqrt(abs((-g/rho0)* dpdz))

dpdz = (rho[1:top+1]-rho[0:top])/dz 
print(dpdz.shape)

N_BV[0:top]  = np.sqrt(abs((-g/rho0)*dpdz))
N_BV[top] = np.sqrt(abs((-g/rho0)*(rho[N-1] - rho[N-2])/(dz)))



def calculate_sm(gh):
    num = B1**(-1/3) - A1*A2*gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))
    dem = (1-3*A2*gh*(B2+6*A1))*(1-9*A1*A2*gh)
    Sm  = num/dem
    return Sm

def calculate_sh(gh):
    return A2*(1-6*A1/B1)/(1-3*A2*gh*(B2+6*A1))

def calculate_brunt_vaisala(rho_, N_BV):
    for i in range(0,N-1):
        dpdz = (rho_[i+1] - rho_[i])/dz     # Density gradient 
        N_BV[i] = math.sqrt(abs((-g/rho0)* dpdz))
    N_BV[N-1] = math.sqrt(abs((-g/rho0)*(rho[N-1] - rho[N-2])/(dz)))
    return N_BV
###########################################################################
####  Set initial conditions for u, q^2, and other turbulence quantities

# Q -- square root of turbulent kinetic energy [sqrt(Q^2)]
# sm, sh are stability parameters, calculated in the mellor-yamada closure. coefficient on turbulent mixing coefficients
# sm is for momentum , sh is for scalars -- 
# nu_t is turbulent viscosity 
# Kq is turbulent diffusion cofficient (for Q^2)
# Kz is turbulent diffusion coefficient for scalars
###########################################################################
Q2  = SMALL*np.ones(N)   # "seed" the turbulent field with small values, then let it evolve
Q2L = SMALL*np.ones(N)

L = -kappa*H*(z/H)*(1-(z/H)) # Q2L(n,1)/Q2(n,1) = 1 at initialization

# Initialize empty arrays with size N
empty_arrays = [np.zeros(N) for i in range(6)]
Q, Sm, Sh, nu_t, Kq, Kz = empty_arrays

# Initialize U 
U = 0.1*(z + 2)
U[z<=-2] = 0
Q = np.sqrt(Q2)

# Gh is a stratification correction 
Gh = -((N_BV*L)/(Q + SMALL))**2
Gh = np.clip(Gh, -0.28, 0.0233)

# Calculate Sm, Sh, nu_t, Kq, Kz
Sm = calculate_sm(Gh) 
Sh = calculate_sh(Gh) 
nu_t = (Sm * Q * L) + nu # Turbulent diffusivity for Q2
Kq = (Sq * Q * L) + nu # Turbulent viscosity
Kz = (Sh * Q * L) + nu


n_profiles = int(M/isave)

variables_to_save = ['U', 'C', 'Q2', 'Q2L', 'rho', 'L', 'nu_t', 'Kz', 'Kq', 'N_BV']

# Intialize an object for saving profiles throughout the model run
saved_profiles = lib.SavedProfiles(n_profiles, variables_to_save, N, isave)  

# Save initial condition (first profile at time zero)
saved_profiles.save_profile_at_timestep(0, 0, U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV)

# Store z in our object so we can plot the profiles later 
saved_profiles.store_z(z)


'''
 Time-advancing algorithm
    Steps a single timestep for c, rho, q2, q2l, l, kz, nu_t, kq
    All diffusion/viscous terms handled implicitly
'''

def TDMA(aX, bX, cX, dX, N):
    # Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
    # a = Lower Diag, b = Main Diag, c = Upper Diag, d = solution vector
    x = np.zeros(N)
    for i in range(1, N):
        bX[i] = bX[i] - aX[i]/bX[i-1]*cX[i-1]
        dX[i] = dX[i] - aX[i]/bX[i-1]*dX[i-1]
    x[-1] = dX[-1]/bX[-1]
    for i in range(N-2, -1, -1):
        x[i] = (1/bX[i])*(dX[i] - cX[i]*x[i+1])
    return x



def wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV):

    # Pre-define Tridiagonal Arrays 
    # Initialize tridiagonal arrays for C/temperature and dC is the RHS vector
    aC, bC, cC, dC = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for turbulent kinetic energy
    aQ2, bQ2, cQ2, dQ2 = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for Q^2 * L (turbulent kinetic energy times a lengthscale)
    aQ2L, bQ2L, cQ2L, dQ2L = lib.initialize_abcd(N)

    # Initialize tridiagonal arrays for velocity
    aU, bU, cU, dU = lib.initialize_abcd(N)

    Px = np.zeros(N)
    Q = np.sqrt(Q2)

    # Update pressure forcing term for the current timestep
    if T_Px == 0.0:
        Px = Px + Px0 # Steady and constant forcing for now
    else: 
        Px = Px + Px0*math.cos(2*math.pi*t[m]/(3600*T_Px))

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

#***************************************************************************
#   Store last time's step variables (f --> fp, q2 --> q2p, etc)
#***************************************************************************
    Cp = C
    Q2p = Q2
    Q2Lp = Q2L
    Lp = L
    Kzp = Kz
    nu_tp = nu_t
    Kqp = Kq
    N_BVp = N_BV
    Up = U
    Vp = V

#***************************************************************************
#   Advance velocity (U,V)
#***************************************************************************
    aU[1:N-1] = -0.5*beta*(nu_tp[1:N-1] + nu_tp[0:N-2])
    bU[1:N-1] = 1 + 0.5*beta*(nu_tp[2:N] + 2*nu_tp[1:N-1] + nu_tp[0:N-2])
    cU[1:N-1] = -0.5*beta*(nu_tp[1:N-1] + nu_tp[2:N])
    dU[1:N-1] = Up[1:N-1] - dt*Px[1:N-1]

    # Bottom boundary: log-law
    bU[0] = 1 + 0.5*beta*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
    cU[0] = -0.5*beta*(nu_tp[1] + nu_tp[0])
    dU[0] = Up[0] - dt*Px[0]

    # Top boundary: no stress
    aU[top] = -0.5*beta*(nu_tp[top]+nu_tp[top-1])
    bU[top] = 1+0.5*beta*(nu_tp[top]+nu_tp[top-1])
    dU[top] = Up[top] - dt*Px[top]

    # Use Thomas algorithm to solve for U
    U = lib.TDMA(aU, bU, cU, dU, N)

#***************************************************************************
#   Advance scalars/density (C, rho) 
#***************************************************************************
    aC[1:N-1] = -0.5*beta*(Kzp[1:N-1] + Kzp[0:N-2])
    bC[1:N-1] = 1 + 0.5*beta*(Kzp[2:N] + 2*Kzp[1:N-1] + Kzp[0:N-2])
    cC[1:N-1] = -0.5*beta*(Kzp[1:N-1] + Kzp[2:N])
    dC[1:N-1] = Cp[1:N-1]

    # Bottom-Boundary: no flux for scalars
    bC[0] = 1+0.5*beta*(Kzp[1] + Kzp[0])
    cC[0] = -0.5*beta*(Kzp[1] + Kzp[0])
    dC[0] =  Cp[0]

    # Top-Boundary: no flux for scalars
    aC[top] = -0.5*beta*(Kzp[top] + Kzp[N-2])
    bC[top] = 1+0.5*beta*(Kzp[top] + Kzp[N-2])
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
    diss = (2 * dt *(Q2p[1:N-1]**0.5))/(B1*Lp[1:N-1])
    aQ2[1:N-1] = -0.5*beta*(Kqp[1:N-1] + Kqp[0:N-2])
    bQ2[1:N-1] = 1 + 0.5*beta*(Kqp[2:N] + 2*Kqp[1:N-1] + Kqp[0:N-2]) + diss 
    cQ2[1:N-1] = -0.5*beta*(Kqp[1:N-1] + Kqp[2:N])
    dQ2[1:N-1] = Q2p[1:N-1] + 0.25*beta*nu_tp[1:N-1]*(Up[2:N]-Up[0:N-2])**2 - dt*Kzp[1:N-1]*(N_BVp[1:N-1]**2)

    # Bottom-Boundary Condition 
    Q2bot = B1**(2/3) * ustar**2
    bdryterm = 0.5*beta*Kqp[0]*Q2bot
    dissipation = 2 * dt *((Q2p[0]**0.5)/(B1*Lp[0]))
    bQ2[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + dissipation
    cQ2[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2[0] = Q2p[0] + dt*((ustar**4)/nu_tp[0]) - dt*Kzp[0]*(N_BVp[0]**2) + bdryterm

    # Top-Boundary Condition
    dissipation =  2 * dt *((Q2p[top]**0.5)/(B1*Lp[top]))
    aQ2[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
    bQ2[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kq[top-1]) + dissipation
    dQ2[top] = Q2p[top] + 0.25*beta*nu_tp[top]*((Up[top] - Up[top-1])**2) -4*dt*Kzp[top]*(N_BVp[top]**2)

    # TDMA to solve for q2
    Q2 = TDMA(aQ2, bQ2, cQ2, dQ2, N)

    # Kluge to prevent negative values from causing instabilities
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
    dQ2L[top] = Q2Lp[top] + 0.25*beta*nu_tp[top]*E1*Lp[top]*(Up[top]-Up[N-2])**2 - 2*dt*Lp[-1]*E1*Kzp[top]*(N_BVp[top]**2)
    
    # TDMA to solve for q2
    Q2L = TDMA(aQ2L, bQ2L, cQ2L, dQ2L, N)

    # Prevent negative values in Q2L 
    Q2L[Q2L < 0] = SMALL

    #  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    Q = np.sqrt(Q2)
    L = Q2L/(Q2 + SMALL)

    # Check length scale 
    ind = (L**2)*(N_BV**2) > (0.281*Q2) # Vectorized if-statement 
    if sum(ind) > 0:
        Q2L[ind] = Q2[ind]*math.sqrt(0.281*Q2[ind]/(N_BV[ind]**2 + SMALL))
        L[ind] = Q2L[ind] / Q2[ind]
    L[abs(L) <= zb] = zb

    Kq = Sq*Q*L + nu
    nu_t = Sm*Q*L + nu
    Kz = Sh*Q*L + nu   

    return U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV

####################### 
## Loop through time ##
########################

t = np.zeros((M))
t[1:M] = dt * (np.arange(1,M) - 1)

for m in range(1,M):
   
    if m%100 == 0:
        print('Time step = %d' % m)

    # Because of how Python handles variables compared to MATLAB, we pass the variables as arguments and get them returned
    # Uses BGO/Mellor-Yamada 2-equation closure
    U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV) 

    if m%isave == 0:
        saved_profiles.save_profile_at_timestep(m, t[m], U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV)

saved_profiles.plot_profiles('U', skip=2)
saved_profiles.plot_profiles('Q2', skip=2)
saved_profiles.plot_profiles('C', skip=2)
saved_profiles.plot_profiles('nu_t', skip=2)