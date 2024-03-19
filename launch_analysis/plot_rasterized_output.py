

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
csv = "../pressure_vs_ws_pmax0.05old.csv"
data = pd.read_csv(csv)

print(len(data.pressure))
print(len(data.ws))
print(len(data.output))

size = (35,35)
# x = np.reshape(data.pressure, size)
# y = np.reshape(data.ws, size)
# z = np.reshape(data.output, size)
plt.scatter(data.pressure, data.ws, c=data.output) #, cmap=‘viridis’)  # cmap specifies the color map
# plt.colorbar(label=‘Z values’)  # Add color bar with label
ax = plt.gca() 
ax.set_xlabel("Pressure")
ax.set_ylabel("ws")
# plt.title(‘Raster-Style Plot with Gridded Data)
fig = plt.gcf()
fig.savefig('test.png')
