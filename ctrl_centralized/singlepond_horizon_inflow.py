from gurobipy import *
import numpy  as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")

# Create a model 
m = Model("sp_in")

# Create intial volume 
InitVolume = 500
inflow = 20

# Create Variable
horizon = [i for i in range(0,1000)]
outflow = m.addVars(horizon, lb=0, vtype=GRB.CONTINUOUS, name="valve")
volume = m.addVars(horizon, lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="volume")

# Define the constrat for updating volume
m.addConstr(volume[0] == InitVolume, "InitialConditions")
m.addConstrs(volume[time]  == volume[time-1] - outflow[time] + 6.0 for time in horizon[1:])

# Define the constraint for maintaning water under threshold 
m.addConstrs(outflow[time] <= 10 for time in horizon)

# Define the constraint for maintaning water under threshold 
m.addConstrs(outflow[time] <= (volume[time-1]*6.14647148e-02 + 3.81)*0.447 for time in horizon[1:])

m.setObjective(volume.sum(), GRB.MINIMIZE)
m.optimize()

volume = []
outflow = []
for v in m.getVars():
    if v.varName.split("[")[0] == 'volume':
        volume.append(v.x)
    else:
        outflow.append(v.x)

print("Obj:", m.objVal)

plt.subplot(1,2,1)
plt.plot(outflow)
plt.ylabel("Outflow")

plt.subplot(1,2,2)
plt.plot(volume)
plt.ylabel("Volume")

plt.show()
