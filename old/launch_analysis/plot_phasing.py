
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 

fn="../tidalphase_vs_windphase_pressure=2e-6_wind=3.5.csv"
title = fn


data = pd.read_csv(fn)

tidal = data.tidal_phasing
wind = data.wind_phasing


base_case = data.loc[(data.tidal_phasing==0) & (data.wind_phasing==0)]

bc_diatoms = base_case.diatom_biomass.values
bc_habs = base_case.hab_biomass.values


diatoms = 100*(data.diatom_biomass.values - bc_diatoms)/bc_diatoms
habs = 100*(data.hab_biomass.values - bc_habs)/bc_habs


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,8), sharex=True, sharey=True)

# fig = plt.figure(figsize=(4,8))
scale = -1
i = axs[0].scatter(tidal, wind, c=diatoms, marker= "s", 
                cmap = mpl.cm.seismic, vmin=-10**(-scale), vmax=10**(-scale), s=62)
plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))
# cbar = plt.colorbar(i, shrink = 0.5, orientation="vertical" , ax=axs[0])#, label = units)
# cbar.set_label(units)

j = axs[1].scatter(tidal, wind, c=habs, marker= "s", 
                cmap = mpl.cm.seismic, vmin=-10**(-scale), vmax=10**(-scale), s=62)
axs[0].set_title("Diatoms")
axs[1].set_title("HABs")
# ticks0 = [0, 0.25 , 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
# cbar.set_ticks(ticks0)
# cbar.set_ticklabels(["All diatoms died", "25% left",  "50% left", "75% left", "No population change",
#                         "125% growth", "150% growth", "175% growth", "Population doubles"])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.35, 0.03, 0.3])
cbar= fig.colorbar(i, shrink = 0.25, cax=cbar_ax)
cbar.set_label("% Shift")
# cbar = fig.colorbar(i, shrink = 0.5, ax=axs[1], orientation="horizontal",pad=0.01)

for ax in axs: 
    ax.set_xlabel("tidal phasing")
    ax.set_ylabel("wind phasing")
    # ax.set_title(title)
    ax.grid(False)
    # plt.tight_layout()
    ax.set_ylim(0,23)
    ax.set_xlim(0,)
# plt.show()
# fn = "test"
print(fn)
fig.savefig("%s.png" % fn)