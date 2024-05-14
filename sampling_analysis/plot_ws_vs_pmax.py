
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import os

# Get list of csv files in a folder
folder = "../launch_analysis"
all_files = os.listdir(folder)    
files = list(filter(lambda f: f.endswith('.csv'), all_files))

for file in files: 
    if not "="  in file:
        continue 

    title = file.replace(".csv", '')
    csv = "../launch_analysis/%s" % file 
    data = pd.read_csv(csv)

    fig = plt.figure(figsize=(7,5))

    output = data.output

    i = plt.scatter(data.ws, data.pmax, c=output, marker= "s", 
                    cmap = mpl.cm.seismic, vmin=0, vmax=2, s=80)

    ax = plt.gca() 
    ax.set_xlabel("$w_s$")
    ax.set_ylabel("pmax")


    ax.set_title(title.replace('_', ' '))
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
    fig.savefig("%s.png" % title)
