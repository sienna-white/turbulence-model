import numpy as np 
import os 
import pandas as pd

# Set range for the two variables of interest 
points = 35
pressure=np.linspace(2e-7,1e-5, num=points)
ws = np.linspace(-1e-4, 1e-4, num=points)

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