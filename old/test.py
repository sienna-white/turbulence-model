

import matplotlib.pyplot as plt 
import pandas as pd 


# from sklearn.linear_model import LinearRegression

stime="2018-01-01"
etime="2018-04-01"

stime = os.getenv("stime")
etime = os.getenv("etime")

syear=int(stime.split("-")[0])
smonth=int(stime.split("-")[1])
sday=int(stime.split("-")[2])

eyear=int(etime.split("-")[0])
emonth=int(etime.split("-")[1])
eday=int(etime.split("-")[2])


print(syear, smonth, sday)
assert(False)
csv="/Users/siennawhite/research/turbulence-model/output/ws_vs_pmax_pressure=2.2e-6.csv"

data0 = pd.read_csv(csv)
data = data0.loc[abs(data0.output)<0.1]

x = data.ws.values
y = data.pmax.values
regression = skl.linear_model.QuadraticRegression() 
regression.fit(x, y)

# extract the data within a threshold +/- 0.5 

''' 
DataFrame.where(cond, other=nan, *, inplace=False, axis=None, level=None)
cond would be range : btwn -0.1, 0.1

--> i get the indices of the data
--> pull out the x, y values 

# do a regression on the x,y values --> sklearn.linear_model.LinearRegression 

# or PolynomialRegression(degree=2) or QuadraticRegression()
regression = linear_model.LinearRegression(degree=2) 
regression.fit(x, y)


experiment with seaborn library for contour plots ... i should look into this 


'''