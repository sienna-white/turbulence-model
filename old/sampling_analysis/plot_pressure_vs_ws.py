
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 


files = ['pressure_vs_ws_pmax=0.1_1SS_diurnallight_tidal.csv',
'pressure_vs_ws_pmax=0.01_1SS_constantlight.csv' ,  
'pressure_vs_ws_pmax=0.1_1SS_constantlight.csv' ,  
'pressure_vs_ws_pmax=0.01_1SS_diurnallight.csv'  ,       
'pressure_vs_ws_pmax=0.01_1SS_diurnallight_tidal.csv'  ,
'pressure_vs_ws_pmax=0.01_1SS_tidal.csv'     ,           
'pressure_vs_ws_pmax=0.05_1SS_constantlight.csv' ,
'pressure_vs_ws_pmax=0.05_1SS_diurnallight.csv',
'pressure_vs_ws_pmax=0.05_1SS_diurnallight_tidal.csv',
'pressure_vs_ws_pmax=0.1_1SS_tidal.csv']


# fn="pressure_vs_ws_pmax=0.1_1SS_diurnallight"

for fn in files: 
    fn = fn.replace(".csv", "")
    title = fn.replace('_', ' ')

    csv = "../launch_analysis/%s.csv" % fn 
    data = pd.read_csv(csv)

    fig = plt.figure(figsize=(7,5))


    output = data.output

    # i = plt.plot(data.ws, data.pressure, 'o')
    i = plt.scatter(data.ws, data.depth_averaged_kz, c=output, marker= "s", 
                    cmap = mpl.cm.seismic, s=50, vmin=0, vmax=2) # norm=mpl.colors.LogNorm()

    black = abs(abs(data.output)-1) < 1e-10
    # plt.scatter(data.ws[black], data.pressure[black], color="black", marker= "s", s=40)

    ax = plt.gca() 
    ax.set_xlabel("$w_s$ [m/s]")
    # ax.set_ylim(0, 0.025)
    ax.set_ylabel("depth-averaged dissipation [m$^2$/s]")
    # ax.set_yscale('log')

    ax.set_title(title)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
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
    print(fn)