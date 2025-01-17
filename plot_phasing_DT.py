
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import cmocean as cmo

fn="temp_vs_windphase_pressure=2e-6_wind=0.csv"

title = fn


data = pd.read_csv(fn)

temp = data.del_temp
wind = data.wind_phasing


# base_case = data.loc[(data.tidal_phasing==0.0) & (data.wind_phasing==0.0)]
# print(base_case)
# bc_diatoms = base_case.diatom_biomass.values
# bc_habs = base_case.hab_biomass.values

bc_diatoms = 0.4
bc_habs = 0.4


# print("bc_diatoms", bc_diatoms)
# if len(bc_diatoms) > 1:
#     bc_diatoms = bc_diatoms[0]
#     bc_habs = bc_habs[0]


diatoms = 100*(data.diatom_biomass.values - bc_diatoms)/bc_diatoms
habs = 100*(data.hab_biomass.values - bc_habs)/bc_habs

print(diatoms[0:15])
print(sum(diatoms>0))
print(habs[0:15])
print(sum(habs>0))

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,8), sharex=True, sharey=True)

vlim = 20 
print(np.mean(diatoms))
i = axs[0].scatter(temp, wind, c=diatoms, marker= "s", 
                cmap = cmo.cm.balance, vmin=-vlim, vmax=vlim, s=66)
# plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))

j = axs[1].scatter(temp, wind, c=habs, marker= "s", 
                cmap = cmo.cm.balance, vmin=-vlim, vmax=vlim, s=66)


# i = axs[0].scatter(tidal[diatoms>0], wind[diatoms>0], c="crimson", marker= "s", s=60)
# j = axs[0].scatter(tidal[diatoms<0], wind[diatoms<0], c="skyblue", marker= "s", s=60)


# j = axs[1].scatter(tidal[habs>0], wind[habs>0], c="crimson", marker= "s", s=60)
# j = axs[1].scatter(tidal[habs<0], wind[habs<0], c="skyblue", marker= "s", s=60)



# scale = 1
# i = axs[0].scatter(tidal, wind, c=diatoms, marker= "s", 
#                 cmap = cmo.cm.haline, vmin=-73, vmax=-68, s=62)
# plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))


# j = axs[1].scatter(tidal, wind, c=habs, marker= "s", 
#                 cmap = cmo.cm.haline, vmin=-47, vmax=-43, s=62)




axs[0].set_title("Diatoms")
axs[1].set_title("HABs")

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.35, 0.03, 0.3])
cbar= fig.colorbar(i, shrink = 0.25, cax=cbar_ax)
cbar.set_label("% Shift")
# cbar = fig.colorbar(i, shrink = 0.5, ax=axs[1], orientation="horizontal",pad=0.01)

for ax in axs: 
    ax.set_xlabel("Temperature difference (top/bottom)")
    ax.set_ylabel("Wind phasing")
    # ax.set_title(title)
    ax.grid(False)
    # plt.tight_layout()
    ax.set_ylim(0,23)
    ax.set_xlim(0,)
# plt.show()
# fn = "test"
print(fn)
fig.savefig("%s.png" % fn)