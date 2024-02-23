
# Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
isave = 30

# Spatial Parameters 
N = 80    # number of grid points
H = 20    # depth (meters)
dz = H/N  # grid spacing - may need to adjust to reduce oscillations
dt = 60   # (seconds) size of time step 
M = 1800 # 400  # number of time steps 

# Physical parameters 
z0 = 0.01         # Bottom roughness [m]
zb = 10*z0        # Bottom height [m]
g  = 9.81         # Gravity [m/s^2]
C_D = 0.0025      # Friction coefficient 
SMALL = 1e-6      # Noise floor for turbulence quantities [m/s?]
kappa = 0.4       # Von Karman constant
nu = 1e-6         # Kinematic viscosity [m^2/s]
rho0 = 1000       # Water density [kg/m^3]

# Algae parameters 
background_turbidity =  0.16
I_in = 350 

diatoms = Algae_Species(k = 0, #0.7,    # specific light attenuation coefficient [cm^2 / 10^6 cells]
                   pmax = 0.05,     # maximum specific growth rate [1/hour]
                   ws = 200,          # vertical velocity [cm/hour]
                   Hi = 40,         # half-saturation of light-limited growth [mu mol photons * m^2/s]
                   Li = 0.006,      # specific loss rate [1/hour]
                   name = "Diatoms")

diatoms.set_initial_concentration(N, init=100, opt='constant')
diatoms.set_vertical_grid(H, N, dz)
algae = diatoms.c
courant = abs(diatoms.ws * dt)/dz
assert(courant < 1)

# Initial conditions for temperature profile
delC   = 5       # Change in temperature at initial themocline [deg C]; set to zero for Unstratified Case
zdelC  = -5      # Position of initial thermocline
dzdelC = 4       # Thickness of initial thermocline 
alpha  = 0.0     # Thermal expansivity, set to zero for passive scalar case
base_temp = 15   # Temperature of water column [deg C]

# Pressure Forcing -> Need to modify to allow for time variable Px.
Px0 = 0.0001/100  # Magnitude on pressure gradient forcing
T_Px = 0 # 12.0  # Period [hours] on pressure gradient forcing. Set to 0 for steady