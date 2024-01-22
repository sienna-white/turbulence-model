import numpy as np 

class Turbulence_Model:
    def __init__(self):
        print('test')
    # Functions for setting model parameters 
    def set_grid_points(self, N):
        # Set the number of vertical grid points for the model
        self.N = N
    def set_grid_depth(self, H):
        # Set the depth of the water column
        self.H = H  
    def set_grid_spacing(self):
        # Set the spacing between grid points
        self.dz = self.H / self.N   
    def set_time_step(self, dt):
        # Set the time step
        self.dt = dt
    def set_num_time_steps(self, M):    
        # How long should the model run for?
        self.M = M
        print('Model will run for %d seconds' % (self.M * self.dt))


    def initilize_grid(self):
        self.z = [(-self.H + self.dz*(i + 0.5)) for i in range(self.N)]

    def calculate_gradient(self,vector):
        gradient = (vector[1:-1] - vector[0:-2])/self.dz
        return gradient 
    

        # vector[]

