from gurobipy import *
import numpy  as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

# Create a model 
m = Model("threeponds_timeoftravel")

# Create intial volume 
InitialVolume = 500
# Basins beign controlled 
basins = [0, 1, 2]

# Time horizon for solving
horizon = [i for i in range(0,600)]

#inflows 
inflows = np.zeros((600))
inflows[:200] = np.sin(np.linspace(0,1.0,200)*np.pi) * 5.0
#inflows[:200] = np.ones(200)*1.0

print(2*inflows.sum())
# Create volumes for each basin 
volumes = m.addVars(horizon, basins, ub=1000, name="volume")

# Create Valves for each basin 
outflows = m.addVars(horizon, basins, name="outflows")

# Create Constraints 
## 1. Max volume capacity constraint already enforced

## 2. Flow limit in channels 
m.addConstrs((outflows[time, 0] <= 10 for time in horizon), name="Upstream Channel1")
m.addConstrs((outflows[time, 1] <= 10 for time in horizon), name="Upstream Channel2")
m.addConstrs((outflows[time, 2] <= 5 for time in horizon), name="Downstream Channel")

## 3. Flow limit based on the volume in the ponds
m.addConstrs((outflows[time, 0] <= (volumes[time-1, 0]*1.13047751*0.01 +  3.65949363) for time in horizon[1:]), name="Upstream Pond1")
m.addConstrs((outflows[time, 1] <= (volumes[time-1, 1]*1.13047751*0.01  + 3.65949363) for time in horizon[1:]), name="Upstream Pond2")
m.addConstrs((outflows[time, 2] <= (volumes[time-1, 2]*1.13047751*0.01  + 3.65949363) for time in horizon[1:]), name="Downstream Pond")

## 4. Add some initial volume 
m.addConstr(volumes[0, 0] == 0.6*InitialVolume, "Initial Volume Pond 1")
m.addConstr(volumes[0, 1] == 0.4*InitialVolume, "Initial Volume Pond 2")
m.addConstr(volumes[0, 2] == InitialVolume, "Initial Volume Pond 3")

## 5. Volume constraints 
m.addConstrs((volumes[time, 0]  == volumes[time-1, 0] + inflows[time] - outflows[time, 0]  for time in horizon[1:]), name="Mass Balance pond:1")
m.addConstrs((volumes[time, 1]  == volumes[time-1, 1] + inflows[time] - outflows[time, 1]  for time in horizon[1:]), name="Mass Balance pond:2")
m.addConstrs(((volumes[time, 2]  == volumes[time-1, 2]  - outflows[time, 2] if time <=20 else volumes[time, 2]  == volumes[time-1, 2] - outflows[time, 2] + outflows[time-20, 0] + outflows[time-20, 1]) for time in horizon[1:]), name="Mass Balance pond:3")

## 5. Minimize the water in the network ponds
m.setObjective(volumes.sum(), GRB.MINIMIZE)
m.optimize()

# Parse using regex 
data = {}
data["volume_0"] = []
data["volume_1"] = []
data["volume_2"] = []

data["outflows_0"] = []
data["outflows_1"] = []
data["outflows_2"] = []

for v in m.getVars():
    temp = re.split('\[|,|\]', v.varName)
    data[temp[0] + "_" + temp[2]].append(v.x)

plt.figure(1)
plt.subplot(3,3,1)
plt.plot(inflows)
plt.ylabel("Inflow")
plt.title("P1")
plt.ylim([0, 20])

plt.subplot(3,3,2)
plt.plot(inflows)
plt.title("P2")
plt.ylim([0, 20])


plt.subplot(3,3,3)
plt.plot(np.asarray(data["outflows_1"]) + np.asarray(data["outflows_0"]))
plt.title("P3")
plt.ylim([0, 20])

plt.subplot(3,3,4)
plt.plot(data["volume_0"])
plt.ylabel("Volume")
plt.ylim([0, 1100])

plt.subplot(3,3,5)
plt.plot(data["volume_1"])
plt.ylim([0, 1100])


plt.subplot(3,3,6)
plt.plot(data["volume_2"])
plt.ylim([0, 1100])

plt.subplot(3,3,7)
plt.plot(data["outflows_0"])
plt.ylabel("Outflows")
plt.ylim([0, 15])

plt.subplot(3,3,8)
plt.plot(data["outflows_1"])
plt.ylim([0, 15])

plt.subplot(3,3,9)
plt.plot(data["outflows_2"])
plt.ylim([0, 15])

plt.show()
