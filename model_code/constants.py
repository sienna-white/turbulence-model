variable2name = {} 
variable2name['U'] = 'velocity'
variable2name['C'] = 'temperature'
variable2name['Q2'] = 'turbulent kinetic energy'
variable2name['Q2L'] = 'turbulent kinetic energy times a length scale'
variable2name['rho'] = 'density'
variable2name['L'] = 'length scale'
variable2name['nu_t'] = r'turbulent viscosity ($\nu_t$)'
variable2name['Kz'] = r'turbulent diffusivity ($\kappa_z$)'
variable2name['Kq'] = r'turbulent diffusivity ($\kappa_q$)'
variable2name['N_BV2'] = r'Brunt-Vaisala frequency ($N_{BV}^2$)'
variable2name['algae'] = r'Algae concentration'
variable2name['biomass'] = r'Total algae biomass'
variable2name['net_growth'] = r'Net Algae Growth (loss + growth)'

variable2name['algae1'] = r'Algae concentration'
variable2name['biomass1'] = r'Total algae biomass'
variable2name['net_growth1'] = r'Net Algae Growth (loss + growth)'
variable2name['algae2'] = r'Algae concentration'
variable2name['biomass2'] = r'Total algae biomass'
variable2name['net_growth2'] = r'Net Algae Growth (loss + growth)'

variable2units = {}
variable2units['U'] = 'm/s'
variable2units['C'] = 'deg C'
variable2units['Q2'] = r'tke'
variable2units['Q2L'] = r'tke*L'
variable2units['rho'] = r'kg/m$^3$'
variable2units['L'] = r'm'
variable2units['nu_t'] = r'm$^2$/s'
variable2units['Kz'] = r'm$^2$/s'
variable2units['Kq'] = r'm$^2$/s'
variable2units['N_BV2'] = r'1/s$^2$'
variable2units['algae'] = r'10$^6$ cells/cm$^3$' #'10$^6$ cells/m$^2$'
variable2units['biomass'] = r'10$^6$ cells'
variable2units['net_growth'] = r'hour$^{-1}$'
variable2units['algae1'] = r'10$^6$ cells/cm$^3$'
variable2units['algae2'] = r'10$^6$ cells/cm$^3$'
variable2units['net_growth1'] = r's$^{-1}$'
variable2units['net_growth2'] = r's$^{-1}$'


variable2convert = {}
for key in variable2units.keys():
    variable2convert[key] = 1.0
variable2convert["algae1"] = 1e6 * (4* 1e-6) * 1000  # 10^6 cells/mL --> ug/L 
variable2convert["algae2"] = 1e6 * (4* 1e-6) * 1000  
variable2convert["algae"] = 1e6 * (4* 1e-6) * 1000  


# Physical parameters 
z0 = 0.01         # Bottom roughness [m]
zb = 10*z0        # Bottom height [m]
g  = 9.81         # Gravity [m/s^2]
C_D = 0.0025      # Friction coefficient 
SMALL = 1e-6      # Noise floor for turbulence quantities [m/s?]
kappa = 0.4       # Von Karman constant
nu = 1e-6         # Kinematic viscosity [m^2/s]
rho0 = 1000       # Water density [kg/m^3]
alpha  =  2.1e-4      # Thermal expansivity, set to zero for passive scalar case
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
# Lisa thinks there's a chnace one of these parameters is off ... 

#********************** FIXED CONSTANTS  ***************************
rhoA = 1.23  # DENSITY OF AIR, kg / m^3
rhoW = 1000  # Density of Water
specific_heat_water = 4181 # J/kg-degC
specific_heat_air = 1007 # J/kg-degCxrh
c_d = 0.05   # Drag coefficient 