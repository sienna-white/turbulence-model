import numpy as np
import math
import pandas as pd 
from constants import * 
# def generate_forcings(self, pressure=0, wind=0, temp=0, light=0, fn=None):
#     if fn is not None:
#         print("File provided with external forcings!")
#         self.read_forcings_from_file(self, fn)



def get_pressure_at_timestep(self, time):
    Px = np.zeros(self.N)
    # return Px 
    # Update pressure forcing term for the current timestep
    if self.T_Px == 0.0:
        Px = Px + self.Px0 # Steady and constant forcing for now
    else: 
        Px = Px + self.Px0*math.cos((2*math.pi*time - 3600*3)/(3600*self.T_Px)) 
    return Px 

def read_forcings_from_file(self, filename):
    # Read in the forcings from a file 
    data = pd.read_csv(filename)
    print("Reading time-varying forcings from file %s" % filename)  
    forcings = {}
    forcings['I_in'] = data['Sol Rad (PAR)'].values
    forcings['T_air'] = data['Air Temp (C)'].values
    forcings['Wind_speed'] = data['Wind Speed (m/s)'].values
    self.forcings = forcings

def air_temperature(self, t, max_temp, diurnal):
    '''Function to generate estimate of heat flux according to diurnal cycle.
    Inputs: t (in seconds); flux_max (maximum heat intensity/ light at noon)'''
    return self.forcings['T_air'][int(t)]
    # if diurnal:
    #     hour = t/3600
    #     period = (2*np.pi)/24
    #     phase_shift = 12 
    #     temp = max_temp * np.cos(period * (hour - phase_shift))
    #     temp = temp.clip(min=18) # During night, light is zero
    #     return temp
    # else:
    #     return temp
    

def diurnal_light(self, t, I_max, diurnal):
    '''Function to generate estimate of light according to diurnal cycle.
    Inputs: t (in seconds); I_max (maximum light intensity/ light at noon)'''
    return self.forcings['I_in'][int(t)]

    # if diurnal:
    #     hour = t/3600
    #     period = (2*np.pi)/24
    #     phase_shift = 12 
    #     light = I_max * np.cos(period * (hour - phase_shift))
    #     light = light.clip(min=0) # During night, light is zero
    #     return light
    # else:
    #     return I_max

def wind(self, t):
    w = self.forcings["Wind_speed"][int(t)] 
    wind = (c_d * w)**2 * rhoA 
    return wind