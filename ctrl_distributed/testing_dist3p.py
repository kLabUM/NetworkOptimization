## Plotting Libraries 
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
import copy 
import numpy as np
from gurobipy import *

def to_numpy(gurobi_variable):
    temp = [gurobi_variable[i].x for i in gurobi_variable.keys()]
    return np.asarray(temp)

def MasterProblem(horizon=20, init_vol=[1000, 1000, 1000], inflow_vol=[0, 0, 0]):
    # Create Model 
    m = Model("MP")

    upper_volume = 1000 # Volume limits.

    basins  = [i for i in range(0, 3)]
    max_threshold  = [10*horizon, 10*horizon, 5*horizon]
    # Create Basin Volume
    volume_basins = m.addVars(basins, lb=0, ub=upper_volume, vtype=GRB.CONTINUOUS, name="Volume_Basin")
    # Create Volume Exchange 
    volume_links  = m.addVars(basins, lb=0, ub=max_threshold, vtype=GRB.CONTINUOUS, name="Volume_Link")
    # Constraints 
    ## 1. Link and basin upper bound constraint set in variable. 
    ## 2. Mass Balance constraint 
    m.addConstr(volume_basins[0] == init_vol[0] - volume_links[0] + inflow_vol[0], "P0")
    m.addConstr(volume_basins[1] == init_vol[1] - volume_links[1] + inflow_vol[1], "P1")
    m.addConstr(volume_basins[2] == init_vol[2] + volume_links[0] + inflow_vol[2] + volume_links[1] - volume_links[2], "P1")
    m.setObjective(volume_basins.sum(), GRB.MINIMIZE)
    m.optimize()
    return to_numpy(volume_links), to_numpy(volume_basins)

def Subproblem(Volume_Links, init_vol, inflows, inflow_p12, horizon_time=20):
    # Common Containers and variables
    obj = GRB.MAXIMIZE
    horizon  = [time for time in range(0, horizon_time)]  # Time steps

    # Pond1 -> Generate outflows
    ## Variables 
    m_p1       = Model("Pond1")
    volume_p1  = m_p1.addVars(horizon, lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="Volume")
    outflow_p1 = m_p1.addVars(horizon, lb=0, ub=10 , vtype=GRB.CONTINUOUS, name="Outflows")

    # Constarints 
    m_p1.addConstr(volume_p1[0]      == init_vol[0] - outflow_p1[0] + inflow_p12[0], "Inital Volume")
    m_p1.addConstr(outflow_p1.sum()  <= Volume_Links[0], "Total Mass Movement")
    m_p1.addConstrs((volume_p1[time] == volume_p1[time-1] - outflow_p1[time] + inflow_p12[time] for time in horizon[1:]), name="Mass balance")
    m_p1.addConstrs(outflow_p1[time] <= (volume_p1[time-1]*1.13047751*0.01 +  3.65949363) for time in horizon[1:])

    # Objective 
    m_p1.setObjective(outflow_p1.sum(), obj)
    m_p1.optimize()


    # Pond2 -> Generate Outflows 
    ## Variables 
    m_p2       = Model("Pond2")
    volume_p2  = m_p2.addVars(horizon, lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="Volume")
    outflow_p2 = m_p2.addVars(horizon, lb=0, ub=10 , vtype=GRB.CONTINUOUS, name="Outflows")

    # Constarints 
    m_p2.addConstr(volume_p2[0]      == init_vol[1] - outflow_p2[0]  + inflow_p12[0] , "Inital Volume")
    m_p2.addConstr(outflow_p2.sum()  <= Volume_Links[1], "Total Mass Movement")
    m_p2.addConstrs((volume_p2[time] == volume_p2[time-1] - outflow_p2[time]  + inflow_p12[time] for time in horizon[1:]), name="Mass balance")
    m_p2.addConstrs(outflow_p2[time] <= (volume_p2[time-1]*1.13047751*0.01 +  3.65949363)  for time in horizon[1:])

    # Objective 
    m_p2.setObjective(outflow_p2.sum(), obj)
    m_p2.optimize()


    # Pond3 -> Take in the outflows and plan how you want to release them.
    ## Generate the outflows from the above into the downstream 
    travel_time = 10
    inflow_p1    = to_numpy(outflow_p1)
    inflow_p2    = to_numpy(outflow_p2)
    inflows[travel_time:] = inflows[travel_time:] + inflow_p1 + inflow_p2  # Travel time inserted here. 

    m_p3       = Model("Pond3")
    volume_p3  = m_p3.addVars(horizon, lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="Volume")
    outflow_p3 = m_p3.addVars(horizon, lb=0, ub=5 , vtype=GRB.CONTINUOUS, name="Outflows")
    # Constarints
    m_p3.addConstr(volume_p3[0]      == init_vol[2] + inflows[0] - outflow_p3[0], "Inital Volume")
    m_p3.addConstr(outflow_p3.sum()  <= Volume_Links[2], "Total Mass Movement")
    m_p3.addConstrs((volume_p3[time] == volume_p3[time-1] - outflow_p3[time] + inflows[time] for time in range(1,20)), name="Mass balance")
    m_p3.addConstrs(outflow_p3[time] <= (volume_p3[time-1]*1.13047751*0.01 +  3.65949363)  for time in horizon[1:])

    # Objective 
    m_p3.setObjective(outflow_p3.sum(), obj)
    m_p3.optimize()

    data = {}
    data["P1"] = {"volume":to_numpy(volume_p1), "outflow":to_numpy(outflow_p1)}
    data["P2"] = {"volume":to_numpy(volume_p2), "outflow":to_numpy(outflow_p2)}
    data["P3"] = {"volume":to_numpy(volume_p3), "outflow":to_numpy(outflow_p3)}
    error = [Volume_Links[0] - data["P1"]["outflow"].sum(), Volume_Links[1] - data["P2"]["outflow"].sum(), Volume_Links[2] - data["P3"]["outflow"].sum()]
    return data, inflows, error


# Execute 
v_p = [500, 500, 500]
hor = 20 
inflows = np.zeros(30)

#inflows 
inflows_p12 = np.zeros((600))
#inflows_p12[:200] = np.sin(np.linspace(0,1.0,200)*np.pi) * 2.0

# Master 
v_link, v_pi = MasterProblem(hor, v_p, [inflows_p12[:20].sum(), inflows_p12[:20].sum(), 0.0 ])
# Inflow exchange 
temp = inflows[20:]
inflows = np.zeros(30)
inflows[:10] = temp
# Subploblem 
data, inflows, er = Subproblem(v_link, v_p, inflows, inflows_p12[:20])

v_p = [data["P1"]["volume"][-1], data["P2"]["volume"][-1], data["P3"]["volume"][-1]]

tor = 20 
DATA = [data]

inf = inflows
errors1 = [er[0]]
errors2 = [er[1]]
errors3 = [er[2]]

for i in range(0, 30):
    # Master 
    v_in = [inflows_p12[tor:tor+20].sum(), inflows_p12[tor:tor+20].sum(), 0.0 ]
    v_link, v_pi = MasterProblem(hor, v_p)
    # Inflow exchange 
    temp = inflows[20:]
    inflows = np.zeros(30)
    inflows[:10] = temp
    # Subploblem 
    inf = np.append(inf, inflows)
    data, inflows, er = Subproblem(v_link, v_p, inflows, inflows_p12[tor:tor+20])
    v_p = [data["P1"]["volume"][-1], data["P2"]["volume"][-1], data["P3"]["volume"][-1]]
    DATA.append(data)
    errors1.append(er[0])
    errors2.append(er[1])
    errors3.append(er[2])


# Stack the things
volume_1 = DATA[0]["P1"]["volume"] 
for i in range(1, len(DATA)):
    volume_1 = np.append(volume_1, DATA[i]["P1"]["volume"])
    
volume_2 = DATA[0]["P2"]["volume"] 
for i in range(1, len(DATA)):
    volume_2 = np.append(volume_2, DATA[i]["P2"]["volume"])
    
volume_3= DATA[0]["P3"]["volume"] 

for i in range(1, len(DATA)):
    volume_3 = np.append(volume_3, DATA[i]["P3"]["volume"] )

outflow_1 = DATA[0]["P1"]["outflow"] 
for i in range(1, len(DATA)):
    outflow_1 = np.append(outflow_1, DATA[i]["P1"]["outflow"] )
    
outflow_2 = DATA[0]["P2"]["outflow"] 
for i in range(1, len(DATA)):
    outflow_2 = np.append(outflow_2, DATA[i]["P2"]["outflow"] )
    
outflow_3= DATA[0]["P3"]["outflow"] 
for i in range(1, len(DATA)):
    outflow_3 = np.append(outflow_3, DATA[i]["P3"]["outflow"] )

plt.subplot(3,3,1)
plt.plot(inflows_p12, linewidth=2.0)
plt.title("P1")
plt.ylabel("Inflow")
plt.ylim([0, 20])


plt.subplot(3,3,2)
plt.plot(inflows_p12, linewidth=2.0)
plt.title("P2")
plt.ylim([0, 20])


plt.subplot(3,3,3)
plt.plot(outflow_1 + outflow_2 , linewidth=2.0)
plt.title("P3")
plt.ylim([0, 20])



plt.subplot(3,3,4)
plt.plot(volume_1, linewidth=2.0)
plt.ylabel("Volume")

plt.subplot(3,3,5)
plt.plot(volume_2, linewidth=2.0)

plt.subplot(3,3,6)
plt.plot(volume_3, linewidth=2.0)

plt.subplot(3,3,7)
plt.plot(outflow_1, linewidth=2.0)
plt.ylim([0, 20])
plt.ylabel("Outflows")

plt.subplot(3,3,8)
plt.plot(outflow_2, linewidth=2.0)
plt.ylim([0, 20])

plt.subplot(3,3,9)
plt.plot(outflow_3, linewidth=2.0)
plt.ylim([0, 20])

plt.figure(2)
plt.plot(np.cumsum(outflow_3))

plt.figure(3)
plt.plot(errors1)
plt.plot(errors2)
plt.plot(errors3)

plt.show()
