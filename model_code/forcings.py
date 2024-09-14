import numpy as np
import math
import pandas as pd 
from constants import * 


def generate_forcings(self, pressure=0, wind=0, temp=0, light=0, fn=None):
    if fn is not None:
        print("File provided with external forcings!")
        self.read_forcings_from_file(self, fn)
    # if fn is None: 
    #     self.light = self.diurnal_light(self.time, light)


def diurnal_light(self, t, I_max, phase_shift =12, diurnal=True):
    '''Function to generate estimate of light according to diurnal cycle.
    Inputs: t (in seconds); I_max (maximum light intensity/ light at noon)'''
    # return self.forcings['I_in'][int(t)]
    if diurnal:
        hour = t/3600
        period = (2*np.pi)/24
        light = I_max * np.cos(period * (hour - phase_shift))
        light = light.clip(min=0) # During night, light is zero
        return light
    else:
        return I_max

def get_pressure_at_timestep(self, time, phase_shift=3):
    Px = np.zeros(self.N)
    time = time/3600
    period = (2*np.pi)/self.T_Px
    Px0 = self.Px0
    # Update pressure forcing term for the current timestep
    if self.T_Px == 0.0:
        Px = Px + Px0 # Steady and constant forcing for now
    else: 
        Px = Px + Px0*np.cos(period * (time - phase_shift)) 
    return Px 


def wind_speed(self, t, bottom_speed, top_speed, phase_shift = 12):
    hour = t/3600
    period = (2*np.pi)/24
    difference = (top_speed - bottom_speed) #/2
    wind_speed = np.cos(period * (hour - phase_shift))*difference  + (top_speed - difference)
    wind_speed = wind_speed.clip(min=bottom_speed) # During night, light is zero
    return wind_speed

def temperature(self, t, bottom_temp, top_temp, phase_shift = 12):
    hour = t/3600
    period = (2*np.pi)/24
    difference = (top_temp - bottom_temp) #/2
    temp = np.cos(period * (hour - phase_shift))*difference  + (top_temp - difference)
    temp = temp.clip(min=bottom_temp) # During night, light is zero
    return temp


def analytical_temperature_profile(self, t, bottom_temp, top_temp):
    hour = t/3600
    period = (2*np.pi)/24
    phase_shift = 12 
    temp = top_temp * np.cos(period * (hour - phase_shift))
    temp = temp.clip(min=18) # During night, light is zero

    x = np.linspace(-0.99, 0.99, self.N) 
    y = np.arctanh(x)
    difference = (top_temp - bottom_temp)/4 
    y = difference*y + (top_temp + bottom_temp)/2 
    y[0:17] = y[17]
    return y



# def get_pressure_at_timestep(self, time):
#     Px = np.zeros(self.N)
#     # return Px 
#     # Update pressure forcing term for the current timestep
#     if self.T_Px == 0.0:
#         Px = Px + self.Px0 # Steady and constant forcing for now
#     else: 
#         Px = Px + self.Px0*math.cos((2*math.pi*time - 3600*3)/(3600*self.T_Px)) 
#     return Px 

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
    # return self.forcings['T_air'][int(t)]
    # if diurnal:
    #     hour = t/3600
    #     period = (2*np.pi)/24
    #     phase_shift = 12 
    #     temp = max_temp * np.cos(period * (hour - phase_shift))
    #     temp = temp.clip(min=18) # During night, light is zero
    #     return temp
    # else:
    #     return temp
    


def wind(self, t):
    w = self.forcings["Wind_speed"][int(t)] 
    wind = (c_d * w)**2 * rhoA 
    return wind