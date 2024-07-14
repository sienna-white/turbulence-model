
import xarray as xr 
import sys
sys.path.append("../../")
import importlib 
import numpy as np 
import watercolumn_run_lib as wrl 
importlib.reload(wrl)

import matplotlib.pyplot as plt

import cmocean as cmo

# cmdict = cmo.tools.get_dict(cmo.cm.phase, N=11)
# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=cmdict)

# stem = "Px=6"


def filter_nan(arr):
    return arr[:, ~np.isnan(arr).any(axis=0)]
stem = "0T_0W_0D"
title ="no tidal, no wind, no diurnal"

zlevel = 0
strat = wrl.WCRun(None, None)
strat.read_from_file('../../output/stratified_model_%s.nc' % stem)

unstrat = wrl.WCRun(None, None) 
unstrat.read_from_file('../../output/unstratified_model_%s.nc' % stem)


def interpret_dataset_string(label):
    run_type = ""
    # if "unstrat" in label:
    #     run_type += "Unstratified, "
    # else:
    #     run_type += "Stratified, "
    if "0T" in label:
        run_type += "constant pressure, "
    else:
        run_type += "tidal, "
    if "0W" in label:
        run_type += "no wind, "
    else:
        run_type += "wind, "
    # if "0D" in label:
    #     run_type += "constant light,"
    # else:
    #     run_type += "diurnal light,"
    return run_type



def get_xhi_from_dataset(ds): 
    dataset = ds.dataset.isel(time=slice(0, -24))

    diffusivity = dataset.Kq.dropna(dim="time")
    z = diffusivity.z.values
    time = diffusivity.time
    diffusivity = diffusivity.values # filter_nan(strat.dataset.Kq.values)

    array_of_z = np.zeros(diffusivity.shape) + z[:, None]

 
    photic_depth = dataset.photic_depth.values
    photic_depth = photic_depth[~np.isnan(photic_depth)]
    photic_depth = photic_depth[0:len(diffusivity.T)]

    # What cells are in the photic zone
    is_in_photic_zone= array_of_z> -2

    masked_diffusivity = np.ma.array(diffusivity, mask=is_in_photic_zone)

    mean_photic_diffusivity = np.ma.mean(masked_diffusivity, axis=0)
    
    xhi1 = mean_photic_diffusivity / (photic_depth * ds.dataset.attrs['ws1'])
    xhi2 = mean_photic_diffusivity / (photic_depth * ds.dataset.attrs['ws2'])

    settling_time_scale = abs(photic_depth / ds.dataset.attrs['ws1']) 
    diffusive_time_scale = photic_depth**2 / mean_photic_diffusivity


    ind = np.isnan(settling_time_scale)
    time = time[~ind]
    settling_time_scale = settling_time_scale[~ind]
    diffusive_time_scale = diffusive_time_scale[~np.isnan(diffusive_time_scale)]
    xhi1 = xhi1[~ind]
    xhi2 = xhi2[~ind]

    return time, abs(xhi1), abs(xhi2), settling_time_scale, diffusive_time_scale


if 0: 

    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True,sharey = "row",figsize=(15, 10))
    ax = ax.ravel()

    labels = ["0T_0W_0D" ,"0T_1W_0D", "1T_1W_0D", "1T_0W_0D"] #,  "0T_0W_1D",  "1T_1W_1D"]
    colors = ['blue', 'red', 'green', 'purple', 'orange']

    for i,label in enumerate(labels):
        s1 = "stratified_model_%s.nc"   % label 
        s2 = "unstratified_model_%s.nc" % label
        strat.read_from_file('../../output/' + s1)
        unstrat.read_from_file('../../output/' + s2)

        time, xhi1, xhi2, st, dt = get_xhi_from_dataset(strat)
        run_label = interpret_dataset_string(s1)
        
        ax[0].plot(time/3600, xhi1, '--', label=r"Diatoms %s" % run_label, linewidth = 2, color = colors[i])

        algae1 = strat.dataset.biomass1 #.dropna(dim="time")
        grad_al = np.gradient(np.flip(algae1.values))
        # ax[2].plot(strat.dataset.time/3600, algae1, '--' ,label=r"$\nabla$ biomass %s" % run_label, linewidth = 2, color = colors[i])

        algae1 = strat.dataset.biomass2 #.dropna(dim="time")
        grad_al = np.gradient(np.flip(algae1.values))
        # ax[2].plot(time/3600, algae1, '.',label=r"HABS", color = colors[i])

        time, xhi1, xhi2, st, dt = get_xhi_from_dataset(unstrat)
        run_label = interpret_dataset_string(s2)
        ax[1].plot(time/3600, xhi1, label=r"Diatoms %s" % run_label, linewidth = 2, color = colors[i])

        algae1 = unstrat.dataset.biomass1.dropna(dim="time")
        grad_al = np.gradient(np.flip(algae1.values))
        # ax[3].plot(time/3600, algae1 ,label=r"$\nabla$ biomass %s" % run_label, linewidth = 2, color = colors[i])

        algae1 = unstrat.dataset.biomass2.dropna(dim="time")
        grad_al = np.gradient(np.flip(algae1.values))
        # ax[3].plot(time/3600, algae1, '-.',label=r"HABS", linewidth = 2, color = colors[i])

    # ax.plot(time, xhi2, label=r"habs, unstratified $\chi$")

    # Set log scale for y axis
    # ax[0].set_yscale('log')
    # ax[0].set_yscale('log')

    ax[0].set_ylabel(r"$\chi$ ")
    ax[1].set_ylabel(r"$\chi$ ")

    ax[2].set_ylabel(r"$\nabla$ biomass")
    ax[3].set_ylabel(r"$\nabla$ biomass")

    for a in ax:
        a.grid(alpha=0.3)
        a.legend()  
    fig.savefig("test.png")


########################################


fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True,sharey = "row",figsize=(15, 6))
ax = ax.ravel()

labels = ["0T_0W_0D" ,"0T_1W_0D", "1T_1W_0D", "1T_0W_0D"]#,  "0T_0W_1D",  "1T_1W_1D"]
colors = ['blue', 'red', 'green', 'purple', 'orange']

for i,label in enumerate(labels):
    s1 = "stratified_model_%s.nc"   % label 
    s2 = "unstratified_model_%s.nc" % label
    strat.read_from_file('../../output/' + s1)
    unstrat.read_from_file('../../output/' + s2)

    time, xhi1, xhi2, st, dt = get_xhi_from_dataset(strat)
    run_label = interpret_dataset_string(s1)

    ax[0].plot(time/3600, st, '--', label=r"Settling $t$ [%s]" % run_label, linewidth = 2, color = colors[i])
    ax[0].plot(time/3600, dt, '-', label=r"Diffusive $t$ [%s]" % run_label, linewidth = 3, color = colors[i])
    print(time/3600)
    ax[0].set_title("Stratified")


    algae1 = strat.dataset.biomass1.dropna(dim="time")
    grad_al = np.gradient(np.flip(algae1.values))
    # ax[2].plot(time/3600, algae1, '--' ,label=r"$\nabla$ biomass %s" % run_label, linewidth = 2, color = colors[i])

    algae1 = strat.dataset.biomass2.dropna(dim="time")
    grad_al = np.gradient(np.flip(algae1.values))
    # ax[2].plot(time/3600, algae1, '.',label=r"HABS", color = colors[i])

    time, xhi1, xhi2, st, dt = get_xhi_from_dataset(unstrat)
    run_label = interpret_dataset_string(s2)
    ax[1].plot(time/3600, st, '--', label=r"Settling $t$ [%s]" % run_label, linewidth = 2, color = colors[i])
    ax[1].plot(time/3600, dt, '-', label=r"Diffusive $t$ [%s]" % run_label, linewidth = 3, color = colors[i])
    ax[1].set_title("Unstratified")   

    algae1 = unstrat.dataset.biomass1.dropna(dim="time")
    grad_al = np.gradient(np.flip(algae1.values))
    # ax[3].plot(time/3600, algae1 ,label=r"$\nabla$ biomass %s" % run_label, linewidth = 2, color = colors[i])

    algae1 = unstrat.dataset.biomass2.dropna(dim="time")
    grad_al = np.gradient(np.flip(algae1.values))
    # ax[3].plot(time/3600, algae1, '-.',label=r"HABS", linewidth = 2, color = colors[i])

# ax.plot(time, xhi2, label=r"habs, unstratified $\chi$")

# Set log scale for y axis
ax[0].set_yscale('log')
ax[0].set_yscale('log')

ax[0].set_ylabel(r"Time scale ")
ax[1].set_ylabel(r"Time scale ")

# ax[2].set_ylabel(r"$\nabla$ biomass")
# ax[3].set_ylabel(r"$\nabla$ biomass")

for a in ax:
    a.grid(alpha=0.3)
    a.legend(loc='lower center', bbox_to_anchor=(0.5, -0.55),
          ncol=2, fancybox=False, shadow=False)
plt.tight_layout()
fig.savefig("time_scales.png")





fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True,sharey = "row",figsize=(15, 6))
ax = ax.ravel()

labels = ["0T_0W_0D" ,"0T_1W_0D", "1T_1W_0D", "1T_0W_0D"]#,  "0T_0W_1D",  "1T_1W_1D"]
colors = ['blue', 'red', 'green', 'purple', 'orange']

for i,label in enumerate(labels):
    s1 = "stratified_model_%s.nc"   % label 
    s2 = "unstratified_model_%s.nc" % label
    strat.read_from_file('../../output/' + s1)
    unstrat.read_from_file('../../output/' + s2)

    time, xhi1, xhi2, st, dt = get_xhi_from_dataset(strat)
    run_label = interpret_dataset_string(s1)

    ax[0].plot(time/3600, st/dt, '-', label=r"$\chi$ [%s]" % run_label, linewidth = 2, color = colors[i])
    print(time/3600)
    ax[0].set_title("Stratified")

    algae1 = strat.dataset.biomass1.dropna(dim="time")
    # grad_al = np.gradient(np.flip(algae1.values))
    ax[2].plot(algae1.time/3600, algae1, '--' ,label=r"Algae biomass [%s]" % run_label, linewidth = 2,alpha = 0.75, color = colors[i])

    algae1 = strat.dataset.biomass2.dropna(dim="time")
    # grad_al = np.gradient(np.flip(algae1.values))
    ax[4].plot(algae1.time/3600, algae1, '-',label=r"HABS biomass [%s]" % run_label , color = colors[i])

    time, xhi1, xhi2, st, dt = get_xhi_from_dataset(unstrat)
    run_label = interpret_dataset_string(s2)
    ax[1].plot(time/3600, st/dt, '-', label=r"$\chi$ [%s]" % run_label, linewidth = 2, color = colors[i])
    ax[1].set_title("Unstratified")   

    algae1 = unstrat.dataset.biomass1.dropna(dim="time")
    # grad_al = np.gradient(np.flip(algae1.values))
    ax[3].plot(algae1.time/3600, algae1 , '--', label=r"Algae biomass [%s]" % run_label, linewidth = 2, alpha = 0.75, color = colors[i])

    algae1 = unstrat.dataset.biomass2.dropna(dim="time")
    grad_al = np.gradient(np.flip(algae1.values))
    ax[5].plot(algae1.time/3600, algae1, '-',label=r"HABS biomass [%s]" % run_label, linewidth = 2, color = colors[i])

# ax.plot(time, xhi2, label=r"habs, unstratified $\chi$")

# Set log scale for y axis
ax[0].set_yscale('log')
ax[4].set_yscale('log')

ax[2].set_yscale('log')
ax[2].set_ylim(0, 1e3)
ax[4].set_ylim(0, 1e3)

ax[0].set_ylabel(r"$\chi$")
ax[1].set_ylabel(r"$\chi$ ")

ax[2].set_ylabel(r"$\nabla$ biomass")
ax[3].set_ylabel(r"$\nabla$ biomass")

for a in ax:
    a.grid(alpha=0.3)
    a.legend(loc='lower center', bbox_to_anchor=(0.5, -0.45),
          ncol=2, fancybox=False, shadow=False)
plt.tight_layout()
fig.savefig("chi2.png")


assert(False)
def plot1d(ds, variable):
    ds = ds.dataset[variable].dropna(dim="time")
    return ds.time, ds.values
########################################
#           Plot photic depth
########################################
fig = plt.figure(figsize=(8,4))
ax = plt.gca()
ax.grid(alpha=0.3)

x,y = plot1d(strat, "photic_depth")
plt.plot(x,y, '-', linewidth=3, alpha=0.6,  label="Stratified")
x,y = plot1d(unstrat, "photic_depth")
plt.plot(x,y, '--', linewidth=3, alpha=0.6, label="Unstratified")

ax.legend() 
ax.set_ylim()
ax.set_title(title)
fig.tight_layout()
fig.savefig("photic_depth_%s.png" % title)


########################################
#           Plot biomass
########################################
fig = plt.figure(figsize=(8,4))
ax = plt.gca()
ax.grid(alpha=0.3)

x,y = plot1d(strat, "biomass1")
plt.plot(x,y, '-', color='skyblue', linewidth=3, label="Algae biomass (stratified)")
x,y = plot1d(unstrat, "biomass1")
plt.plot(x,y, '--', color='skyblue', linewidth=3, label="Algae biomass (unstratified)")

x,y = plot1d(strat, "biomass2")
plt.plot(x,y, '-', color='crimson', linewidth=3, label="HAB biomass (stratified)")
x,y = plot1d(unstrat, "biomass2")
plt.plot(x,y, '--', color='crimson', linewidth=3, label="HAB biomass (unstratified)")

ax.legend() 
ax.set_ylim()
ax.set_title(title)
fig.tight_layout()
fig.savefig("biomass_%s.png" % title)

########################################
# fig = plt.figure(figsize=(8,4))
# ax = plt.gca()
# ax.grid(alpha=0.3)

# plt.plot(strat.dataset["Kz"].dropna(dim="time").isel(z=zlevel), '-o', label="Stratified")
# plt.plot(unstrat.dataset["Kz"].dropna(dim="time").isel(z=zlevel), '--', label="Untratified")
# ax.legend() 
# ax.set_ylim()


plot_var = [ "C","N_BV", "U", "Kz"] 

fig, axs = plt.subplots(nrows=4, ncols=2, sharex='row', figsize=(10, 20))
axs = axs.ravel()

for i,variable in enumerate(plot_var): 
    fig, _ = strat.plot_profiles(variable, skip=1, passed_string='(stratified)', show=False, ax=axs[i*2])
    fig, _ = unstrat.plot_profiles(variable, skip=1, passed_string='(unstratified)', show=False , ax=axs[i*2+1])

fig.suptitle(title)
fig.savefig("Composite_%s.png" % title)

plot_var = [ "algae1","algae2"] #, , , "algae1"]

fig, axs = plt.subplots(nrows=2, ncols=2, sharex='row', figsize=(10, 10))
axs = axs.ravel()
for i,variable in enumerate(plot_var): 
    fig, ax = strat.plot_profiles(variable, skip=1, passed_string='(stratified)', show=False, ax=axs[i*2])
    fig, ax = unstrat.plot_profiles(variable, skip=1, passed_string='(unstratified)', show=False , ax=axs[i*2+1])
fig.suptitle(title)
fig.tight_layout()
fig.savefig("Algae_%s.png" % title)


# Add biomass 
# axs[0].set_title("Biomass over time")
# label1 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[0].name, ListOfSpecies[0].pmax, ListOfSpecies[0].ws)
# label2 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[1].name, ListOfSpecies[1].pmax, ListOfSpecies[1].ws)
# axs[0] = self.add_time_series_to_axis('biomass1', label1, axs[0])
# axs[0] = self.add_time_series_to_axis('biomass2', label2, axs[0])
# axs[0].legend()




# unstrat.add_profile_to_axis("C", plot_index, legend_ind, ax) # label="Stratified")

