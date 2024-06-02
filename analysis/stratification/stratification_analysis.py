
import xarray as xr 
import sys
sys.path.append("../../")
import importlib 

import watercolumn_run_lib as wrl 
importlib.reload(wrl)

import matplotlib.pyplot as plt


zlevel = 0
strat = wrl.WCRun(None, None)
strat.read_from_file('../../output/stratified_model_nowind.nc')

unstrat = wrl.WCRun(None, None)
unstrat.read_from_file('../../unstratified_model_nowind.nc')



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
plt.plot(x,y, '-', label="Stratified")
x,y = plot1d(unstrat, "photic_depth")
plt.plot(x,y, '-', label="Unstratified")

ax.legend() 
ax.set_ylim()
fig.savefig("photic_depth.png")


########################################
fig = plt.figure(figsize=(8,4))
ax = plt.gca()
ax.grid(alpha=0.3)

plt.plot(strat.dataset["Kz"].dropna(dim="time").isel(z=zlevel), '-o', label="Stratified")
plt.plot(unstrat.dataset["Kz"].dropna(dim="time").isel(z=zlevel), '--', label="Untratified")
ax.legend() 
ax.set_ylim()



plot_var = [ "C","N_BV", "U", "Kz"] #, , , "algae1"]

fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(10, 20))
axs = axs.ravel()
for i,variable in enumerate(plot_var): 
    fig, ax = strat.plot_profiles(variable, skip=11, passed_string='(stratified)', show=False, ax=axs[i*2])
    fig, ax = unstrat.plot_profiles(variable, skip=11, passed_string='(unstratified)', show=False , ax=axs[i*2+1])
fig.savefig("Composite_nowind.png")

assert(False)


# unstrat.add_profile_to_axis("C", plot_index, legend_ind, ax) # label="Stratified")

plt.show()
# # Add biomass 
# axs[0].set_title("Biomass over time")
# label1 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[0].name, ListOfSpecies[0].pmax, ListOfSpecies[0].ws)
# label2 = "Biomass of %s (pmax = %1.1e, ws = %1.1e)" % (ListOfSpecies[1].name, ListOfSpecies[1].pmax, ListOfSpecies[1].ws)
# axs[0] = self.add_time_series_to_axis('biomass1', label1, axs[0])
# axs[0] = self.add_time_series_to_axis('biomass2', label2, axs[0])
# axs[0].legend()