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
import time
import sys



def run_watercolumn(ws, Px0, output_csv):

    t1 = time.time()  # Time our simluation 
    print("Running simulation w/ ws = %e and Px0 = %e" % (ws, Px0))

    def save_at_end(result):
        lib.save_output(output_csv, Px0, ws, result, header=["pressure", "ws", "output"])

    # Spatial Parameters 
    N = 80    # number of grid points
    H = 10    # depth (meters)
    dz = H/N  # grid spacing - may need to adjust to reduce oscillations
    dt = 10 #60   # (seconds) size of time step 
    M  =  1440*3*6 # 400  # number of time steps 

    read_from_input=False
    plot=False
    output=True

    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 1

    # Algae parameters 
    background_turbidity =  0.16
    I_in = 350 

    '''

    Diatoms ws = -1.38e-5 m/s
    Cyanobacteria = 1.38e-4 m/s
    '''
    # Show --> ws=1e-7
    diatoms = Algae_Species(k = 0.07,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                    pmax = 0.05,# 5, #0.1, #0.05,     # maximum specific growth rate [1/hour]
                    ws = ws, #-1.4e-6, #1e-5, #-1e-6, #-1e-9,#-1e-3, #-200,       # vertical velocity [m/s]
                    Hi = 40,         # half-saturation of light-limited growth [mu mol photons * m^2/s]
                    Li = 0.006,      # specific loss rate [1/hour]
                    name = "Diatoms",
                    self_shading=False,
                    net=True)
    init = 200 

    # Pressure Forcing -> Need to modify to allow for time variable Px.
    # Px0 = 2e-6 # 2e-6  # Magnitude on pressure gradient forcing
    T_Px = 0 # 12.0  # Period [hours] on pressure gradient forcing. Set to 0 for steady


    diatoms.set_initial_concentration(N, init=init, opt='linear')
    diatoms.save_total_mass()
    diatoms.set_vertical_grid(H, N, dz)
    algae = diatoms.c
    courant = abs(diatoms.ws * dt)
    assert(courant<dz)
    print(diatoms.pmax)
    RUN_INFO=' growth= %2.2e' % diatoms.pmax
    # RUN_INFO=' ws= %2.2e' % diatoms.ws


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

    # Intialize an object for saving profiles throughout the model run
    n_profiles = int(M/isave)
    variables_to_save = ['U', 'C', 'Q2', 'Q2L', 'rho', 'L', 'nu_t', 'Kz', 'Kq', 'N_BV','algae', 'biomass', 'net_growth']
    saved_profiles = lib.SavedProfiles(n_profiles, variables_to_save, N, isave)  

    # Save initial condition (first profile at time zero)
    data = {'U': U, 'C': C, 'Q2': Q2, 
            'Q2L': Q2L, 'rho': rho, 'L': L,
            'nu_t': nu_t, 'Kz': Kz, 'Kq': Kq,
            'N_BV': N_BV, 'algae': algae, 'biomass': sum(algae), 'net_growth': algae*0}
    saved_profiles.save_profile_at_timestep(0, 0, **data)

    # Store z in our object so we can plot the profiles later 
    saved_profiles.store_z(z)

    def wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae, time):

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

        # Initialize tridiagonal arrays for algae 
        aA, bA, cA, dA = lib.initialize_abcd(N)

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
        Kz = Kz.clip(SMALL,)     # Set floor on Kz so it's never = zero 

        #***************************************************************************
        #   Store last time's step variables (f --> fp, q2 --> q2p, etc)
        #***************************************************************************
        Ap = algae 
        Cp = C
        Q2p,Q2Lp  = Q2, Q2L
        Lp, Kzp, Kqp, nu_tp = L, Kz, Kq, nu_t
        N_BVp = N_BV
        Up, Vp = U, V
        
        #***************************************************************************
        #   Advance velocity (U,V)
        #***************************************************************************
        aU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[0:top-1])
        bU[1:top] = 1 + beta/2*(nu_tp[2:top+1] + 2*nu_tp[1:top] + nu_tp[0:top-1])
        cU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[2:top+1])
        dU[1:top] = Up[1:top] - dt*Px[1:top]

        # Bottom boundary: log-law
        bU[0] = 1 + beta/2*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
        cU[0] = -beta/2*(nu_tp[1] + nu_tp[0])
        dU[0] = Up[0] - dt*Px[0]

        # Top boundary: no stress
        aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
        bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
        dU[top] = Up[top] - dt*Px[top]

        # Use Thomas algorithm to solve for U
        U = lib.TDMA(aU, bU, cU, dU, N)

        #***************************************************************************
        #   Advance algae! // need to figure out to include settling velocity / source + sink
        #***************************************************************************
        ws    = diatoms.ws 
        wsdtdz  = abs(ws*dt)/dz
        light =  350 # lib.diurnal_light(time, 350)
        gamma = diatoms.get_loss_and_growth(I_in = light, current_concentration = Ap)

        # If settling speed is UPWARD (swimming!)
        if ws>0:
            aA[1:top]  = -wsdtdz - beta/2 * (Kzp[0:top-1]+ Kzp[1:top]) 
            bA[1:top]  = 1 + wsdtdz - gamma[1:top]*dt  + beta/2*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1]) 
            cA[1:top]  = -beta/2 * (Kzp[1:top] + Kzp[2:top+1])
            dA = Ap

            # Bottom-Boundary: no flux for scalars
            bA[0] =  1 - (gamma[0]*dt) + beta/2*(Kzp[1] + Kzp[0]) + wsdtdz
            cA[0] = -beta/2 * (Kzp[1] + Kzp[0])
            dA[0] =  Ap[0]

            # Top-Boundary: no flux for scalars
            aA[top] = -wsdtdz - beta/2 * (Kzp[top] + Kzp[top-1])
            bA[top] = 1  - gamma[top]*dt + beta/2 * (Kzp[top] + Kzp[top-1]) # removed + wsdtdz 
            dA[top] = Ap[top]

        # If settling speed is DOWNWARD (sinking!)
        if ws<=0:
            aA[1:top]  = -beta/2 * (Kzp[0:top-1] + Kzp[1:top]) 
            bA[1:top]  = 1 + wsdtdz - gamma[1:top]*dt  + beta/2*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1]) 
            cA[1:top]  = -wsdtdz - beta/2 * (Kzp[1:top] + Kzp[2:top+1])
            dA = Ap

            # Bottom-Boundary: no flux for scalars
            bA[0] =  1 + wsdtdz - (gamma[0]*dt) + beta/2*(Kzp[1] + Kzp[0]) 
            cA[0] =  -wsdtdz -beta/2 * (Kzp[1] + Kzp[0])
            dA[0] =  Ap[0]

            # Top-Boundary: no flux for scalars
            aA[top] =  -beta/2 * (Kzp[top] + Kzp[top-1])
            bA[top] =  1 - gamma[top]*dt + beta/2 * (Kzp[top] + Kzp[top-1]) + wsdtdz # okay adding this here 
            dA[top] = Ap[top]

        # Thomas algorithm to solve for C
        algae = lib.TDMA(aA, bA, cA, dA, N)  
        diatoms.c = algae
        # if ws>0:
        #     algae[top] = algae[top] + Ap[top]*wsdtdz

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

        return U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae, gamma


    #***************************************************************************
    #   LOOP THROUGH TIME ! 
    #***************************************************************************
    for m in range(1,M):

        U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae, gamma = wc_advance(U, C, Q2, Q2L, rho, L, nu_t, Kz, Kq, N_BV, algae, m) 

        diatoms.save_total_mass()
        if m%isave == 0:
            diatoms.save_total_mass()
            # Pack data into dictionary structure before saving 
            data = {'U': U, 'C': C, 'Q2': Q2, 
                    'Q2L': Q2L, 'rho': rho, 'L': L,
                    'nu_t': nu_t, 'Kz': Kz, 'Kq': Kq,
                    'N_BV': N_BV, 'algae': algae, 'biomass': sum(algae),
                    'net_growth' : gamma}
            saved_profiles.save_profile_at_timestep(m, t[m], **data)

    print(time.time()  - t1) 

    if output:
        epsilon = 1e-3
        change = diatoms.total_mass[0] - diatoms.total_mass[-1]
        if change > epsilon:
            save_at_end(False)
            return 
        elif change < epsilon:
            save_at_end(True)
            return 
        if abs(change) < epsilon: 
            save_at_end(np.nan)
            return

