import xarray as xr 
import pandas as pd 

def save_run_info(self, **kwargs):
    # Add to attributes
    if self.save_output:
        self.dataset = self.dataset.assign_attrs(**kwargs)

def save_1d_data(self, time, **kwargs):
    ''' 1D data inclues variables like biomass, forcings & photic depth
        that varies only with time (and not with depth)'''
    if self.save_output:
        df = pd.DataFrame(index=[time]) 
        df.index.name = 'time'
        for key, value in kwargs.items():
            df[key] = value
        dfnc = df.to_xarray() 
        self.dataset = xr.concat([dfnc, self.dataset], dim='time')
    else: 
        return 

def save_2d_data(self, time, **kwargs):
    if self.save_output:
        time_index = [time] * self.N      
        index = pd.MultiIndex.from_arrays([self.z, time_index], names=["z", "t"])
        df = pd.DataFrame(index=index) 
        for key, value in kwargs.items():
            df[key] = value
        dfnc = df.to_xarray()
        # dfnc = (df.unstack('z').to_xarray().to_array('time'))

        if self.first_data: 
            self.dataset= dfnc 
            self.first_data = False
        else:
            self.dataset = xr.concat([dfnc, self.dataset], dim='t')
    else:
        return

def save_dataset(self, filename):
    if self.save_output:
        self.dataset = self.dataset.sortby('t')
        self.dataset.to_netcdf(filename)
    return
