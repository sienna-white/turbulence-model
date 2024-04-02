
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 

fn="pressure_vs_ws_pmax=0.01"
csv = "../output/%s.csv" % fn 
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



true = data[data.output==1]
false = data[data.output==0]
fig = plt.figure(figsize=(5,5))

output = 1 - data.output


i = plt.scatter(abs(data.depth_averaged_u), data.ws, c=output, marker= "s", 
                cmap = mpl.cm.seismic, vmin=0, vmax=1)
# plt.plot(true.pressure, true.ws, "s", color="skyblue", label= "Growing", markersize=9)
# plt.plot(false.pressure, false.ws, "s", color="crimson", label= "Dying", markersize=9)

ax = plt.gca() 
ax.set_xlabel("pressure")
ax.set_ylabel("$w_s$")
# ax.set_xscale('log')
# plt.hlines(0, xmax = max(data.pmax), xmin=min(data.pmax))
ax.set_title("Comparing pressure to growth rate  (pmax=0.03)")

units = "percent change"
cbar = plt.colorbar(i, shrink = 0.4, orientation="horizontal" , label = units)
cbar.set_label(units)

# ax.set_title("Comparing $w_s$ to growth rate  (pressure=2.2e-5)")
plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
ax.grid(False) #alpha = 0) #0.3)
ax.set_xlim(0,10)
# ax.legend()

plt.tight_layout()
# plt.show()
fig.savefig("test.png") #%s.png" % fn)
