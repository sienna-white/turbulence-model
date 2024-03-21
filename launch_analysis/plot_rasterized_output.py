

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
csv = "../ws_vs_pmax_px2.2e-6.csv"
data = pd.read_csv(csv)

print(len(data.pmax))
print(len(data.ws))
print(len(data.output))

# x = np.reshape(data.pressure, size)
# y = np.reshape(data.ws, size)
# z = np.reshape(data.output, size)
plt.scatter(data.pmax, data.ws, c=data.output) #, cmap=‘viridis’)  # cmap specifies the color map
# plt.colorbar(label=‘Z values’)  # Add color bar with label
ax = plt.gca() 
ax.set_xlabel("pmax")
ax.set_ylabel("ws")
# plt.title(‘Raster-Style Plot with Gridded Data)
fig = plt.gcf()
fig.savefig('test.png')
