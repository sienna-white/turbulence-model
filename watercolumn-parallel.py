
import dask 
from dask.distributed import Client, LocalCluster
from watercolumn_fun_pmax import run_watercolumn
import numpy as np 


if __name__ == '__main__':

    # Set up the dask cluster 
    cluster = LocalCluster(n_workers = 64)
    c = Client(cluster)
    tasks = []

    # Set range for the two variables of interest 
    points = 50
    # pressure=np.logspace(1e-7,1e-1, num=points)
    pmax=np.linspace(0.0005,1, num=points)
    ws = np.logspace(-1e-7, 1e-7, num=points)

    output_csv = "ws_vs_pmax_px2.2e-6.csv"
    # output_csv = "pressure_vs_ws_pmax0.05.csv"

    print("There are %d tasks" % len(ws)) 

    for p0 in pmax:
        for ws0 in ws:
            # tasks.append(dask.delayed(run_watercolumn)(p0, ws0, output_csv))

            tasks.append(dask.delayed(run_watercolumn)(ws0, p0, output_csv))

    results = dask.compute(tasks)

    cluster.close()
