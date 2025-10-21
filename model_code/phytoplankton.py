

import importlib  
import numpy as np 
import matplotlib.pyplot as plt

# N = 80    # number of grid points
# H = 20    # depth (meters)
# dz = H/N  # grid spacing - may need to adjust to reduce oscillations
# z = [(-H + dz*(i + 0.5)) for i in range(N)]
# z = np.array(z)

'''
Huisman paper

Diatoms ws = -1.38e-5 m/s
Cyanobacteria = 1.38e-4 m/s
'''

class Algae_Species:
    '''
    k     : specific light attenuation coefficient [cm^2 / 10^6 cells]
    pmax  : maximum specific growth rate [1/hour] --> 
    ws    : vertical velocity [cm/hour]
    Hi    : half-saturation of light-limited growth [mu mol photons * m^2/s]
    Li    : specific loss rate [1/hour]
    Im    : Incident light intensity [mu mol photons/m^2/s]
    K_bkg : background turbidity [1/m]
    z     : depth of water column layer [m]
    '''

    def __init__(self, k=0, pmax = 0, ws=0, Hi = 0, Li = 0, name = None, self_shading=False, net=False):
        ''' Initialize species with provided values'''
        self.k = k * (1e-2)**2 # Convert to m^2/ 10^6 cells? #  0.1 *  0.0001 # k is specific light attenuation coefficient [cm^2 / 10^6 cells]-->m2
        self.pmax = self.hour2second(pmax)      # 1/hour to 1/second 
        self.ws =   ws #self.hour2second(ws)/100    # 100 cm/m ; 3600 seconds/hour 
        self.Hi = Hi
        self.Li = self.hour2second(Li)
        self.name = name
        self.total_mass = [] 
        self.net_growth = []
        self.self_shading=self_shading
        self.net=net

    def set_vertical_grid(self, H, N, dz):
        self.N = N
        self.dz = dz
        z = [(-H + dz*(i + 0.5)) for i in range(N)]
        self.z = np.array(z)

    def set_initial_concentration(self, N, init, opt='constant'):
        ''' Set initial concentration of species
            opt can be set to linear if you'd like to initialize 
            the algae with some sort of gradient. It's pretty 
            hacky as an option right now. '''
        self.c = np.zeros(N) + init 
        if opt == 'linear':
            self.c = np.linspace(0,init,N)
        else: 
            print("Setting constant concentration for %s @ %d" % (self.name, init))

    def hour2second(self, input_rate):
        ''' Convert per-hour rate to per-second rate'''
        return input_rate * (1/3600)

    def monod_growth_rate(self, I): 
        '''Monod growth rate 
        Inputs: 
            pmax : maximum specific growth rate [1/hour]
            I    : light intensity [mu mol photons/m^2/s]
            Hi   : half-saturation of light-limited growth [mu mol photons * m^2/s]
        Output:
            pi   : specific growth rate [1/hour converted to 1/second]
        '''
        pi = self.pmax * I/(self.Hi + I)
        return pi 
    
    def get_light_intensity(self, I_in, other_species=None):
        background_turbidity =  0.26
        if other_species is None:
            I, photic_depth = self_shading([self, other_species], I_in=I_in, turbidity=background_turbidity, self_shading=True)
        else: 
            I, photic_depth = self_shading([self], I_in=I_in, turbidity=background_turbidity, self_shading=True)

        return I, photic_depth # photic_depth = self_shading([self], I_in=I_in, turbidity=background_turbidity)
    
    def get_loss_and_growth(self, I_in, current_concentration):
        ''' Get loss and growth rates for a given
        n time step'''
        if self.net:
            self.c = current_concentration
            # I, photic_depth = self.get_light_intensity(I_in)
            growth = self.monod_growth_rate(I_in) 
            loss = self.Li
            net  = growth - loss
        else:
            net = np.zeros(self.N,)
        return net 
    
    def save_total_mass(self):
        ''' Save total mass of species at each time step'''
        self.total_mass.append(np.sum(self.c)) 

    def save_net_growth(self, net):
        ''' Save net growth of species at each time step'''
        self.net_growth.append(net*self.c)


# Lives outside the class since we need to calculate light intensity for all species 

class SelfShading():
    def __init__(self, z, N, turbidity, self_shading=True):
        self.turbidity = turbidity
        self.self_shading = self_shading 
        self.z  = z
        self.N = N
        self.kxC = np.zeros(self.N)
        self.corrected_z = (-1) * self.z 
        self.add_txz = self.turbidity * self.corrected_z
        self.I = np.zeros(self.N)

    def calc_self_shading(self, ListofAlgae, I_in):
        if self.self_shading:
            for species in ListofAlgae:
                self.kxC += species.k * species.c
        vector = self.kxC*self.corrected_z +  self.add_txz  # This makes sense to be plus to me 
        self.kxC =  self.kxC*0 
        # I = np.zeros(self.N)
        # for i in (range(self.N)):
        #     I[i] = np.sum(vector[i:-1]) 
        # I = I * self.z
  
        I = np.zeros(self.N)
        cumsum_vector = np.cumsum(vector[::-1])[::-1]
        I = cumsum_vector * self.z
        I = I_in * np.exp(I)
        

        # photic_depth =  self.z[I<(0.1 * I_in)][-1] # Get first element where light is less than 90% of I_in
        # try: 
        #     photic_depth =  self.z[I<(0.1 * I_in)][-1] # Get first element where light is less than 90% of I_in
        # except:
        #     photic_depth = 0
        return I #photic_depth


def self_shading(ListofAlgae, I_in, turbidity, self_shading=True):
    ''' Use Lambert-Beer's Law to calculate light intensity at each depth '''

    z  = ListofAlgae[0].z       # m 
    dz = ListofAlgae[0].dz      # m 
    N = ListofAlgae[0].N        # [-]
    # We want water depth [cm] to be zero at the top, 20 at the bottom (pos. numbers)
    corrected_z = (-1) * z 
    kxC = np.zeros(N)
    if self_shading:
        for species in ListofAlgae:
            kxC += species.k * species.c

    I = np.zeros(N)
    vector = kxC + (turbidity*corrected_z) # This makes sense to be plus to me 
    for i in (range(N)):
        I[i] = np.sum(vector[i:-1]) * z[i]

    I0 = I_in * np.exp(I)
    print("When c1 = %f and c2=%f, light=" % (ListofAlgae[0].c[N], ListofAlgae[1].c[N]))
    print(I0)
    assert(False)
    try: 
        photic_depth =  z[I<(0.9 * I_in)][-1] # Get first element where light is less than 90% of I_in
    except:
        photic_depth = 0
    return I0, photic_depth #photic_depth  

    # print("Photic depth is %2.2f meters" % photic_depth)

    # fig = plt.figure()
    # ax = plt.gca()
    # plt.plot(I,z, linewidth = 4, color = 'skyblue')
    # ax.set_xlabel(r'Light Intensity ($\mu$mol photons/m$^2$/s)')
    # ax.set_ylabel('Depth (m)')
    # ax.hlines(photic_depth, 0, I_in, color = 'black', linestyle = '--', label = 'Photic Depth')
    # ax.set_ylim(-20,0)
    # ax.grid(alpha = 0.5)
    # plt.show()

    


  


# k     : specific light attenuation coefficient [cm^2 / 10^6 cells]
# pmax  : maximum specific growth rate [1/hour]
# ws    : vertical velocity [cm/hour]
# Hi    : half-saturation of light-limited growth [mu mol photons * m^2/s]
# Li    : specific loss rate [1/hour]
# Im    : Incident light intensity [mu mol photons/m^2/s]
# K_bkg : background turbidity [1/m]
# z     : depth of water column layer [m]

# diatoms 
# diatoms = Algae_Species(k = 0.7, 
#                    pmax = 0.05,
#                    ws = -5,
#                    Hi = 40,
#                    Li = 0.006,
#                    name = "Diatoms")

# # microcystis 
# hab = Algae_Species(k = 0.034, 
#                     pmax = 0.008,
#                     ws = 50,
#                     Hi = 40,
#                     Li = 0.004,
#                     name = "HAB")

# hab.set_initial_concentration(N, init=0.05)
# diatoms.set_initial_concentration(N, init=10)

# background_turbidity =  0.16
# I_in = 350 


# I, photic_depth = self_shading([diatoms], I_in=I_in, turbidity=background_turbidity)


# fig = plt.figure()
# ax = plt.gca()
# plt.plot(I,z, linewidth = 4, color = 'skyblue')
# ax.set_xlabel(r'Light Intensity ($\mu$mol photons/m$^2$/s)')
# ax.set_ylabel('Depth (m)')
# ax.hlines(photic_depth, 0, I_in, color = 'black', linestyle = '--', label = 'Photic Depth')
# ax.set_ylim(-H,0)
# ax.grid(alpha = 0.5)
# plt.show()
   


