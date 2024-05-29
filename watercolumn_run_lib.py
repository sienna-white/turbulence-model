import pandas as pd 
import xarray as xr 

class Run:
    def __init__(self, N, z):
        self.N = N
        self.z = z
    def store_z(self, z):
        self.z = z
    
    def add_1d_data(self, time, **kwargs):
        ''' 1D data inclues variables like biomass, forcings & photic depth
            that varies only with time (and not with depth)'''
        df = pd.DataFrame(index=[time]) 
        df.index.name = 'time'
        for key, value in kwargs.items():
            df[key] = value
        dfnc = df.to_xarray() 
        print(dfnc)
    
    
    def create_dataframe_from_saved_profiles(self, time, **kwargs):
        time_index = [time] * self.N
        index = pd.MultiIndex.from_arrays([self.z, time_index], names=["z", "time"])
        df = pd.DataFrame(index=index) 
        for key, value in kwargs.items():
            df[key] = value

        dfnc = df.to_xarray()
        print(df.U)
        print(dfnc)