
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

sys.path.append("model_code/")
import watercolumn_lib as lib
from phytoplankton import Algae_Species, SelfShading
from phytoplankton import self_shading
from scipy.interpolate import CubicSpline

import importlib
importlib.reload(lib) # Reload the module if it has changed


import sys, os, time 

t1 = time.time()  # Time our simluation 


import xarray as xr


data = xr.open_dataset("TEST_heatflux_0T_1W_0D.nc")

z = data.z.values
N = len(z)

########################################
dt = 10
Nsteps = 800 
Nsave = 50
Nparticles = 450 
shift = 8
variance=1e-1
print("Total time = %d minutes" % (dt*Nsteps/60))
########################################

# Pull out turbulent dissipation
kz = data.Kz.values
kz = kz[~np.isnan(kz)]
kz = kz.reshape(N, len(kz)//N)

data.close()

Kzs = [10, 20, 28]

# Initialize figure 
fig, axs = plt.subplots(nrows=3, ncols =2, figsize=(6,11), sharey=True)
axs = axs.ravel()

# Initialize time vector 
time = np.arange(0, Nsteps*dt, dt) 

for k,kzN in enumerate(Kzs):
    print('k=',k)

    kz1 = kz[:, kzN]

    # Smooth Kz + fit cubic spline 
    N = 18              # Kernel size for Kz smoothing  ( bigger=more smooth; smaller=less smooth)
    kz1 = np.convolve(kz1, np.ones(N)/N, mode='same')
    cs = CubicSpline(z, kz1)

    # Calculate gradient of kz1 
    dkdz = np.gradient(cs(z))

    axs[k*2].plot(kz1, z, color="skyblue", linewidth=5, label="smoothed modeled $K_z$")
    axs[k*2].plot(cs(z), z, '--', label="spline fit", color="crimson")
    axs[k*2].grid(alpha=0.5)
    axs[k*2].set_ylim(-10,0)
    axs[k*2].set_xlim(0,2e-3)
    axs[k*2].legend()
    axs[k*2].set_title("diffusivity")
    axs[k*2].ticklabel_format(style='sci',scilimits=(-3,4),axis='both')



    # wss = [-1.4e-4, -1.4e-5, 0]
    wss = [0] #, 0 , 0] #[1.4e-4, 1.4e-5, 0]
    shift = 321

    # Initialize particles 
    # init = np.zeros((Nparticles,)) - 5 #shift
    init = np.linspace(-10, 0, Nparticles)
    time4plot = time = np.arange(0, Nsteps*dt, dt*Nsave) 
    colors = ['#268a6d', '#bc7e88', '#d3dd32']
    ax0 = axs[2*k + 1] 

    for w,ws in enumerate(wss):
        zp = init

        # generate random matrix 
        r = np.random.normal(0, variance, size=(Nparticles, Nsteps))

        # generate output matrix 
        path = np.zeros((Nparticles, Nsteps//Nsave))

        path[:,0] = zp 


        for i in range(1,Nsteps):
            # Calculate diffusivity values @ each particle's point 
            kz1 = cs(zp)

            # Calculate dk/dz @ each particle's point 
            dk = cs(zp, 1) 

            # K_corrected term 
            k_corrected = cs(zp + 1/2*dk*dt)

            # Diffusive walk eqn 
            zn = zp + r[:,i]*np.sqrt(2*dt*k_corrected/variance) + ws*dt  + dk*dt
            
            # Implement boundary conditions 
            zn[zn>0] = 0 
            zn[zn<-10] = np.NaN
            zp = zn     
            if i%Nsave==0:
                path[:,i//Nsave] = zn



        h = ax0.plot([], [], color=colors[w], linewidth=4, label="ws=%1.1e" % ws)

        for i in range(Nparticles):
            ax0.plot(time/3600, path[i,:], linewidth=2, alpha=0.05, color=colors[w])#, label="z0=%2.1f" % shift)
    ax0.set_xlabel("hr")
    
    ax0.grid(alpha=0.5)
    ax0.legend()
    # ax0.set_ylim(-2.5,-1.5)
    # plt.legend()
figname = "COMBINED_ini_swimmingt=%1.1f.png" % shift
fig.savefig(figname)
print('saved %s' % figname)

fig, ax = plt.figure(), plt.gca()
ax.hist(path[:,-1], bins=30, density=True, orientation='horizontal')
# ax.hlines(-shift, 0, 5, color='k')
ax.set_ylim(-10, 1)
ax.grid(alpha=0.3)
fig.savefig("hist_%2.1f.png" % shift)





# final = np.zeros((Nparticles,))
# for n,z0 in enumerate(fives): #[-5, -5, -5, -5, -5]: #np.linspace(0,-10, 6): #: #, -7, -9]:
    
#     # path = np.zeros((Nsteps*dt//Nsave,))
#     # path[0] = z0 

#     zp = z0 
#     variance=1e-3
#     r = np.random.normal(0, variance, size=(Nsteps,))
    
#     for i in range(1,Nsteps):
        
#         ind = ind4z(zp)
#         kz = kz1[ind]
#         dk = dkdz[ind]


#         kz2 = kz1[ind4z(zp + 1/2*kz*dt)]
#         zn = zp + (dk*dt) + r[i]*np.sqrt(2*dt*kz2/variance) + ws*dt 
#         if zn>0:
#             zn = 0 
#         if zn<-10:
#             zn = np.NaN

#         # if i%Nsave==0:
#         #     path[i//Nsave] = zn
#         zp = zn     



#     # ax.plot(time/3600, path, linewidth=1, label="z0=%2.1f" % z0)
#     final[n] = zp #path[-1]
# ax.set_xlabel("hr")

# ax.grid(alpha=0.5)
# ax.set_ylim(-10,0)
# # plt.legend()
# fig.savefig("path.png")

# fig, ax = plt.figure(), plt.gca()
# ax.hist(final, bins=50)
# fig.savefig("hist.png")

# for n,z0 in enumerate(fives): #[-5, -5, -5, -5, -5]: #np.linspace(0,-10, 6): #: #, -7, -9]:
    