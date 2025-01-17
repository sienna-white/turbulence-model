import numpy as np 
import os 
import pandas as pd
import itertools 


# Set range for the two variables of interest 
# points = 35
# tidal_phasing = np.arange(0,6,step=0.5)
wind_phasing = np.arange(0, 24,step=0.5)
# options = list(itertools.product(tidal_phasing, wind_phasing))
# tidal_phase, wind_phase = list(map(list, zip(*options)))
# data = {"tidal_phase": tidal_phase, "wind_phase": wind_phase}


temp_strat = np.arange(0, 3,step=0.25)
options = list(itertools.product(temp_strat, wind_phasing))
temp_strat, wind_phase = list(map(list, zip(*options)))
data = {"temp_strat": temp_strat, "wind_phase": wind_phase}

df = pd.DataFrame(data)
df.to_csv("temp_vs_wind_options.csv", index=False)

assert(False)




output_csv = "pressure_vs_ws.csv"


for p0 in pressure:

    for ws0 in ws:
        # print("Pressure = %f, ws = %f" % (p0, ws0))

        data = {"pressure": [p0], "ws": [ws0]}
        print(data)
        df = pd.DataFrame(data)

        # Append the DataFrame to the CSV file
        file_exists = os.path.isfile(output_csv)
        if file_exists:
            mode = 'a'
            header=False
        else:
            mode = 'w'
            header=True

        df.to_csv(output_csv, mode=mode, index=False, header=header)