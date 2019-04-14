from gurobipy import *
import numpy  as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")

# Create a model 
m = Model("twoponds_timeoftravel")

# Create intial volume 
InitialVolume = 500
# Basins beign controlled 
basins = [0, 1]

# Time horizon for solving
horizon = [i for i in range(0,200)]

# Create volumes for each basin 
volumes = m.addVars(horizon, basins, ub=500, name="volume")

# Create Valves for each basin 
outflows = m.addVars(horizon, basins, name="outflows")


# Create Constraints 
## 1. Max volume capacity constraint already enforced

## 2. Flow limit in channels 
m.addConstrs((outflows[time, 0] <= 10 for time in horizon), name="Upstream Channel")
m.addConstrs((outflows[time, 1] <= 5 for time in horizon), name="Downstream Channel")

## 3. Flow limit based on the volume in the ponds
m.addConstrs((outflows[time, 0] <= (volumes[time-1, 0]*6.14647148e-02 + 3.81)*0.447 for time in horizon[1:]), name="Upstream Pond")
m.addConstrs((outflows[time, 1] <= (volumes[time-1, 1]*6.14647148e-02 + 3.81)*0.447 for time in horizon[1:]), name="Downstream Pond")

## 4. Add some initial volume 
m.addConstr(volumes[0, 0] == InitialVolume, "Initial Volume Pond 1")
m.addConstr(volumes[0, 1] == InitialVolume, "Initial Volume Pond 2")

## 5. Volume constraints 
m.addConstrs((volumes[time, 0]  == volumes[time-1, 0] - outflows[time, 0] for time in horizon[1:]), name="Mass Balance pond:0")
m.addConstrs(((volumes[time, 1]  == volumes[time-1, 1] - outflows[time, 1] if time <=10 else volumes[time, 1]  == volumes[time-1, 1] - outflows[time, 1] + outflows[time-3, 0]) for time in horizon[1:]), name="Mass Balance pond:1")

## 5. Minimize the water in the network ponds
m.setObjective(volumes.sum(), GRB.MINIMIZE)
m.optimize()

# Parse using regex 
data = {}
data["volume_0"] = []
data["volume_1"] = []

data["outflows_0"] = []
data["outflows_1"] = []

for v in m.getVars():
    temp = re.split('\[|,|\]', v.varName)
    data[temp[0] + "_" + temp[2]].append(v.x)

plt.subplot(2,2,1)
plt.plot(data["volume_0"])
plt.ylabel("Volume")
plt.title("Upstream")

plt.subplot(2,2,3)
plt.plot(data["outflows_0"])
plt.ylabel("Outflows")

plt.subplot(2,2,2)
plt.plot(data["volume_1"])
plt.title("Downstream")

plt.subplot(2,2,4)
plt.plot(data["outflows_1"])

plt.show()



