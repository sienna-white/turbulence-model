
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 

fn="ws_vs_pmax_pressure=2.2e-5"
title = "Comparing $w_s$ to pmax (pressure=2.2e-4)"
csv = "../output/%s.csv.csv" % fn 
data = pd.read_csv(csv)


# x = np.reshape(data.pressure, size)
# y = np.reshape(data.ws, size)
# z = np.reshape(data.output, size)


# true = data[data.output==1]
# false = data[data.output==0]
# fig = plt.figure(figsize=(5,5))

# plt.plot(true.pmax, true.ws, "s", color="skyblue", label= "Growing", markersize=9)
# plt.plot(false.pmax, false.ws, "s", color="crimson", label= "Dying", markersize=9)

# ax = plt.gca() 
# ax.set_xlabel("pmax (growth rate [1/hr])")
# ax.set_ylabel("$w_s$")
# # plt.hlines(0, xmax = max(data.pmax), xmin=min(data.pmax))
# ax.set_title("Comparing $w_s$ to growth rate  (pressure=2.2e-5)")
# plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
# ax.grid(False) #alpha = 0) #0.3)
# ax.legend()

# plt.tight_layout()
# fig.savefig("ws_vs_pmax_px=2.2e-.png")


fig = plt.figure(figsize=(7,5))

output = data.output


i = plt.scatter(data.ws, data.pmax, c=output, marker= "s", 
                cmap = mpl.cm.seismic, vmin=0, vmax=2, s=90)

ax = plt.gca() 
ax.set_xlabel("$w_s$")
ax.set_ylabel("pmax")


ax.set_title(title)
plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))

cbar = plt.colorbar(i, shrink = 0.9, orientation="vertical" )#, label = units)
# cbar.set_label(units)


ticks0 = [0, 0.25 , 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
cbar.set_ticks(ticks0)
cbar.set_ticklabels(["All diatoms died", "25% loss",  "50% loss", "75% loss", "No population change",
                        "125% loss", "150% growth", "175% loss", "Population doubles"])


ax.grid(False)
plt.tight_layout()
# plt.show()
fig.savefig("%s.png" % fn)
