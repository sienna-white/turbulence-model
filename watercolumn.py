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

turb_model = TM()

N=80 # number of grid points
H=20 # depth (meters)
dz=H/N # grid spacing - may need to adjust to reduce oscillations
dt=60 # (seconds) size of time step 
M=2000 # number of time steps 

beta = (dt/dz**2)

isave=1 #increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
savecount=1

z = [(-H + dz*(i + 0.5)) for i in range(N)]
z = np.array(z)
# Store model information in class object
turb_model.set_grid_points(N)
turb_model.set_grid_depth(H)
turb_model.set_grid_spacing()
turb_model.set_time_step(dt)
turb_model.set_num_time_steps(M)

turb_model.initilize_grid() 

'''
for i in range(N): # Initialize grid
   z[i]=-H+dz*(i + 0.5) # bottom at z=-H, free surface at 0
   # Because of the difference in indentation between Python and MATLAB, we use + 0.5 here instead of - 0.5
'''


'''
Initialize all profiles and closure parameters 
Call once before time loop
Sets all forcing: pressure gradients, stresses, etc.
Should be used to adjust initial temperature/salinity profiles
Velocity initialized to zero
Turbulence quantities initialized to "SMALL"; Lengthscale parabolic

'''



# Physical parameters 
z0 = 0.01         # Bottom roughness [m]
zb = 10*z0        # Bottom height [m]
g  = 9.81         # Gravity [m/s^2]
C_D = 0.0025      # Friction coefficient 
SMALL = 1e-6      # not sure?
kappa = 0.4       # Von Karman constant
nu = 1e-6         # Kinematic viscosity [m^2/s]
rho0 = 1000       # Water density [kg/m^3]


# Forcing parameters - Need to modify to allow for time variable Px.
Px0 = .001  # Magnitude on pressure gradient forcing
T_Px = 12.0 # Period [hours] on pressure gradient forcing. Set to 0 for steady

# Turbulence closure parameters 
# all dimensionless for mellor-yamada / no need to change unless we modify model
A1=0.92
A2=0.74
B1=16.6
B2=10.1
C1=0.08
E1=1.8
E2=1.33
E3=0.25
Sq=0.2 

# Setup initial conditions for scalar and density
delC   = 5       # Change in temperature at initial themocline [deg C]; set to zero for Unstratified Case
zdelC  = -5      # Position of initial thermocline
dzdelC = 4 #-2   # Width of initial thermocline
alpha  = 0.0     # Thermal expansivity, set to zero for passive scalar case
base_temp = 15   # Temperature of water column [deg C]

# Initalize arrays 
empty_arrays = [np.zeros(N) for i in range(5)]
C, rho, N_BV, U, V = empty_arrays


#for i in range(N):
 # if z[i] >= -5:
  #  C[i] = 20
  #else:
  #  C[i] = 15 

half_height_thermocline = 0.5*dzdelC

C = base_temp + delC*(z - zdelC + 0.5*dzdelC)/dzdelC
#   For values of z BELOW the thermocline, set C = base_temp
C[(z <= (zdelC - half_height_thermocline))] = base_temp
#   For values of z ABOVE the thermocline, set C = base_temp + delC
C[(z >  (zdelC + half_height_thermocline))] = base_temp + delC

print(z)
# import matplotlib.pyplot as plt
# plt.plot(C,z)
# ax = plt.gca()
# ax.grid(True)
# plt.show()

# Calculate density profile 
rho = rho0*(1 - alpha*(C - base_temp))  # Single scalar, linear equation of state


# Caclulate Brunt-Vaisala Frequency 
dpdz = turb_model.calculate_gradient(rho) # Calculate vertical density gradient 
# N_BV = np.sqrt(abs((-g/rho0)*dpdz)) 
# N_BV[0] = math.sqrt(abs((-g/rho0)*(rho[1]-rho[0])/(dz)))
for i in range(0,N-1):
  dpdz = (rho[i+1] - rho[i])/dz     # Density gradient 
  N_BV[i] = math.sqrt(abs((-g/rho0)* dpdz))
N_BV[N-1] = math.sqrt(abs((-g/rho0)*(rho[N-1] - rho[N-2])/(dz)))



###########################################################################
####  Set initial conditions for u, q^2, and other turbulence quantities
###########################################################################
Q2  = SMALL*np.ones(N)   # "seed" the turbulent field with small values, then let it evolve
Q2L = SMALL*np.ones(N)

L = -kappa*H*(z/H)*(1-(z/H)) # Q2L(n,1)/Q2(n,1) = 1 at initialization

# Initialize empty arrays with size N
empty_arrays = [np.zeros(N) for i in range(6)]
Q, Sm, Sh, nu_t, Kq, Kz = empty_arrays
# Q -- square root of turbulent kinetic energy [sqrt(Q^2)]
# sm, sh are stability parameters, calculated in the mellor-yamada closure. coefficient on turbulent mixing coefficients
# sm is for momentum , sh is for scalars -- 
# nu_t is turbulent viscosity 
# Kq is turbulent diffusion cofficient (for Q^2)
# Kz is turbulent diffusion coefficient for scalars

# Gh is a stratification correction 

for i in range(N):
  if z[i] <= -2:
    U[i] = 0
  else:
    U[i] = 0.1*(z[i] + 2)
  Q[i] = math.sqrt(Q2[i])


  Gh = -((N_BV[i]*L[i])/(Q[i] + SMALL))**2
  Gh = min(Gh, 0.0233)
  Gh = max(Gh, -0.28)   

  # 
#   y1 = (1/3) - (2*A1/B1)
#   y2 = (B2/B1) + (6*A1/B1)    
#   sm_siena = (A1/A2)*(B1*(y1 - C1) - B1*(y1-C1) + 6*(A1 + ))   
  num= B1**(-1/3) - A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))
  dem= (1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh)
  Sm[i] = num/dem
  Sh[i] = A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1))
  nu_t[i] = Sm[i] * Q[i] * L[i] + nu # Turbulent diffusivity for Q2
  Kq[i] = Sq*Q[i]*L[i] + nu # Turbulent viscosity
  Kz[i] = Sh[i] * Q[i] * L[i] + nu # Turbulent scalar diffusivity




#  Save initial conditions as first columns in saved matrix
n_profiles = int(M/isave)
Um = np.zeros((N,4))
Cm = np.zeros((N,4))
Q2m = np.zeros((N,4))
Q2Lm = np.zeros((N,4))
rhom = np.zeros((N,4))
Lm = np.zeros((N,4))
nu_tm = np.zeros((N,4))
Kzm = np.zeros((N,4))
Kqm = np.zeros((N,4))
N_BVm = np.zeros((N,4))

# Initialize first column 
Um[:,0] = U
Cm[:,0] = C
Q2m[:,0] = Q2
Q2Lm[:,0] = Q2L
rhom[:,0] = rho
Lm[:,0] = L
nu_tm[:,0] = nu_t
Kzm[:,0] = Kz
Kqm[:,0] = Kq
N_BVm[:,0] = N_BV


'''
 Time-advancing algorithm
    Steps a single timestep for c, rho, q2, q2l, l, kz, nu_t, kq
    All diffusion/viscous terms handled implicitly
'''

def TDMA(aX, bX, cX, dX, N):
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
    # Pre-define Tridiagonal Arrays - Just in case
    aC = np.zeros(N)
    bC = np.zeros(N)
    cC = np.zeros(N)
    dC = np.zeros(N)
    # ^ these vectors represent the diagonals of the tridiagonal matrix and dC is the RHS vector
    aQ2 = np.zeros(N)
    bQ2 = np.zeros(N)
    cQ2 = np.zeros(N)
    dQ2 = np.zeros(N)
    # ^ these vectors represent the diagonals of the tridiagonal matrix for turbulent kinetic energy
    aQ2L = np.zeros(N)
    bQ2L = np.zeros(N)
    cQ2L = np.zeros(N)
    dQ2L = np.zeros(N)
    # ^ these vectors represent the diagonals of the tridiagonal matrix for Q^2 * L (turbulent kinetic energy times a lengthscale)



    Px = np.zeros(N)

    #  Update pressure forcing term for the current timestep
    for i in range(N):
        if T_Px == 0.0:
            Px[i] = Px0  # Steady and constant forcing for now
        else:
            Px[i] = Px0*math.cos(2*math.pi*t[m]/(3600*T_Px))

    #Update shear velocity at bottom boundary for use later
    ustar = abs(U[0])*math.sqrt(C_D); # Explicit dependence on C_D


    # Update parameters for the model, Sm and Sh
    for i in range(N):		
        Gh=-(N_BV[i]*L[i]/(Q[i]+SMALL))**2 
        
        # set LIMITER for Gh 
        Gh=min(Gh, 0.0233)
        Gh=max(Gh, -0.28)

        num=B1**(-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))
        dem=(1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh)
        Sm[i]=num/dem
        Sh[i]=A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1)) 

    #  Place previous variable f into fp (i.e. q2 into q2p, etc)
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


    # Advance scalars/density (C, rho) 
    for i in range(1, N-1):
        aC[i] = -0.5*beta*(Kzp[i] + Kzp[i-1])
        bC[i] = 1+0.5*beta*(Kzp[i+1] + 2*Kzp[i] + Kzp[i-1])
        cC[i] = -0.5*beta*(Kzp[i] + Kzp[i+1])
        dC[i] = Cp[i]
    # Bottom-Boundary: no flux for scalars
    bC[0] = 1+0.5*beta*(Kzp[1] + Kzp[0])
    cC[0] = -0.5*beta*(Kzp[1] + Kzp[0])
    dC[0] =  Cp[0]
    # Top-Boundary: no flux for scalars
    aC[-1] = -0.5*beta*(Kzp[-1] + Kzp[N-2])
    bC[-1] = 1+0.5*beta*(Kzp[-1] + Kzp[N-2])
    dC[-1] = Cp[-1]

    # Thomas algorithm to solve for C
    C = TDMA(aC, bC, cC, dC, N)

    #update density and Brunt-Vaisala frequency
    for i in range(N):
        rho[i] = rho0*(1-alpha*(C[i] - 15))
    N_BV[0] = math.sqrt((((-g/rho0)*(rho[1]-rho[0])/dz)))
    for i in range(1, N-1):
        N_BV[i] = math.sqrt(((-g/rho0)*(rho[i+1] - rho[i-1])/(2*dz)))
    N_BV[-1] = math.sqrt(((-g/rho0)*(rho[-1] - rho[N-2])/dz))


    #  Advance turbulence parameters (q2, q2l - q2 first, then q2l)
    for i in range(1, N-1):
            diss = (2 * dt *(Q2p[i]**0.5))/(B1*Lp[i]) # Coefficient for linearized term
            aQ2[i] = -0.5*beta*(Kqp[i] + Kqp[i-1])
            bQ2[i] = 1+0.5*beta*(Kqp[i+1] + 2*Kqp[i] + Kqp[i-1]) + diss
            cQ2[i] = -0.5*beta*(Kqp[i] + Kqp[i+1])
            dQ2[i] = Q2p[i] + 0.25*beta*nu_tp[i]*(Up[i+1]-Up[i-1])**2 -dt*Kzp[i]*(N_BVp[i]**2)

    # Bottom-Boundary Condition 
    Q2bot = B1**(2/3) * ustar**2
    bdryterm = 0.5*beta*Kqp[0]*Q2bot
    diss =  2 * dt *((Q2p[0]**0.5)/(B1*Lp[0]))
    bQ2[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + diss
    cQ2[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2[0] = Q2p[0] + dt*((ustar**4)/nu_tp[0]) - dt*Kzp[0]*(N_BVp[0]**2) + bdryterm

    # Top-Boundary Condition
    diss =  2 * dt *((Q2p[-1]**0.5)/(B1*Lp[-1]))
    aQ2[-1] = -0.5*beta*(Kqp[-1] + Kqp[N-2])
    bQ2[-1] = 1+0.5*beta*(Kqp[-1] + 2*Kqp[-1] + Kq[N-2]) + diss
    dQ2[-1] = Q2p[-1] + 0.25*beta*nu_tp[-1]*((Up[-1] - Up[N-2])**2) -4*dt*Kzp[-1]*(N_BVp[-1]**2)

    # TDMA to solve for q2
    Q2 = TDMA(aQ2, bQ2, cQ2, dQ2, N)
    # Kluge to prevent negative values from causing instabilities
    for i in range(N):
        if Q2[i] < 0:
            Q2[i] = SMALL
    
    

    for i in range(1, N-1):
            diss = 2*dt*((Q2p[i]**0.5) / (B1*Lp[i]))*(1+E2*(Lp[i]/(kappa*abs(-H-z[i])))**2 + E3*(Lp[i]/(kappa*abs(z[i])))**2)
            aQ2L[i] = -0.5*beta*(Kqp[i] + Kqp[i-1])
            bQ2L[i] = 1+0.5*beta*(Kqp[i+1] + 2*Kqp[i] + Kqp[i-1]) + diss
            cQ2L[i] = -0.5*beta*(Kqp[i] + Kqp[i+1])
            dQ2L[i] = Q2Lp[i] + 0.25*beta*nu_tp[i]*E1*Lp[i]*(Up[i+1]-Up[i-1])**2 - 2*dt*Lp[i]*E1*Kzp[i]*(N_BVp[i]**2)
    # Bottom-Boundary Condition
    q2lbot = B1**(2/3) * (ustar**2) * kappa * zb
    bdryterm = 0.5*beta*Kqp[0]*q2lbot
    diss =  2 * dt *(Q2p[0]**0.5)/(B1*Lp[0])*(1+E2*(Lp[0]/(kappa*abs(-H-z[0])))**2 + E3*(Lp[0]/(kappa*abs(z[0])))**2)
    bQ2L[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + diss
    cQ2L[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2L[0] = Q2Lp[0] + dt*((ustar**4)/nu_tp[0])*E1*Lp[0] - dt*Lp[0]*E1*Kzp[0]*(N_BVp[0]**2) + bdryterm

    # Top-Boundary Condition
    diss =  2 * dt *(Q2p[-1]**0.5)/(B1*Lp[-1])*(1+E2*(Lp[-1]/(kappa*abs(-H-z[-1])))**2 + E3*(Lp[-1]/(kappa*abs(z[-1])))**2)
    aQ2L[-1] = -0.5*beta*(Kqp[-1] + Kqp[N-2])
    bQ2L[-1] = 1+0.5*beta*(Kqp[-1] + 2*Kqp[-1] + Kqp[N-2]) + diss # Are we using kq or kqp here?
    dQ2L[-1] = Q2Lp[-1] + 0.25*beta*nu_tp[-1]*E1*Lp[-1]*(Up[-1]-Up[N-2])**2 - 2*dt*Lp[-1]*E1*Kzp[-1]*(N_BVp[-1]**2)
    # TDMA to solve for q2
    Q2L = TDMA(aQ2L, bQ2L, cQ2L, dQ2L, N)
    # Making sure to prevent negative values
    for i in range(N):
        if Q2L[i] < 0:
            Q2L[i] = SMALL
    
    #  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    #     Works will all updated values 
    for i in range(N):
            Q[i] = math.sqrt(Q2[i])
            L[i] = Q2L[i]/(Q2[i] + SMALL)
            
            # Limit due to stable stratification
            if (L[i]**2)*(N_BV[i]**2) > 0.281*Q2[i]:
                # Adjust Q2L as well as L
                Q2L[i] = Q2[i]*math.sqrt(0.281*Q2[i]/(N_BV[i]**2 + SMALL))
                L[i] = Q2L[i] / Q2[i]  
            # Keep L from becoming zero -- zb=bottom roughness parameter
            
            if abs(L[i]) <= zb:
                L[i] = zb
            # Update diffusivities 
            Kq[i] = Sq*Q[i]*L[i] + nu
            nu_t[i] = Sm[i]*Q[i]*L[i] + nu
            Kz[i] = Sh[i]*Q[i]*L[i] + nu      
    #return C, Q2, Q2L
    return U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV


t = np.zeros(M)
for m in range(1,M):
   
    if m%100 == 0:
        print('m = %d' % m)
    t[m]=dt*(m-1) #define time
    # Because of how Python handles variables compared to MATLAB, we pass the variables as arguments and get them returned
    [U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV] = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV) #uses BGO/Mellor-Yamada 2-equation closure
    if m == 1 or m == 10 or m == 140:
        Cm[:,savecount] = C
        Q2m[:,savecount] = Q2
        Q2Lm[:,savecount] = Q2L
        rhom[:,savecount] = rho
        Lm[:,savecount] = L
        nu_tm[:,savecount] = nu_t
        Kzm[:,savecount] = Kz
        Kqm[:,savecount] = Kq
        N_BVm[:,savecount] = N_BV
        Um[:,savecount] = U
        savecount += 1 # Because of how Python handles identation compared to Matlab, 
    #savecount starts at 1 (instead of 0) and only gets incremented at the end of the step
savecount = 1



print('Plotting...')
plt.plot(Um[:,1], z)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12,12))

plt.show()