
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 


fn="pressure_vs_ws_pmax=0.08"
title = fn

csv = "../output/%s.csv" % fn 
data = pd.read_csv(csv)

fig = plt.figure(figsize=(7,5))

output = data.output

i = plt.scatter(data.ws, abs(data.depth_averaged_u), c=output, marker= "s", 
                cmap = mpl.cm.seismic, vmin=0, vmax=2, s=50)

ax = plt.gca() 
ax.set_xlabel("$w_s$ [m/s]")
ax.set_ylabel("depth-averaged velocity [m/s]")


ax.set_title(title)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))

cbar = plt.colorbar(i, shrink = 0.9, orientation="vertical" )#, label = units)
# cbar.set_label(units)


ticks0 = [0, 0.25 , 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
cbar.set_ticks(ticks0)
cbar.set_ticklabels(["All diatoms died", "25% left",  "50% left", "75% left", "No population change",
                        "125% growth", "150% growth", "175% growth", "Population doubles"])


ax.grid(False)
plt.tight_layout()
# plt.show()
# fn = "test"
fig.savefig("%s.png" % fn)