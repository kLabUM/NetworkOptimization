## Plotting Libraries 
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
import copy 
import numpy as np
from gurobipy import *
import copy


def master(data):
    for i in data["nodes"] :
        volumes = [] 
        outflows_t = []
        outflows_tt = []
        
    m = Model("MP")
    m.Params.OutputFlag = 0

    # Volume in basins 
    volume_basins = m.addVars(data["nodes"], lb=0, vtype=GRB.CONTINUOUS)
    # Volume in the links travelling 
    volume_lt = m.addVars(data["arcs"], lb=0, vtype=GRB.CONTINUOUS) # send to other with in this horizon 
    volume_ltt = m.addVars(data["arcs"], lb=0, vtype=GRB.CONTINUOUS) # send in the next horizon

    ## Node 0 + arc 01 
    m.addConstr(volume_basins[0] == data["vt0"][0] - volume_lt[0] - volume_ltt[0])
    m.addConstr(volume_lt[0] <= (data["horizon"] - data["horizon"]%data["travel_time"][0])*data["channel_flow"][0])
    m.addConstr(volume_ltt[0] <= (data["horizon"]%data["travel_time"][0])*data["channel_flow"][0])
    m.addConstr(volume_lt[0] + volume_ltt[0] <= data["horizon"] * data["channel_flow"][0])
    m.addConstr(volume_basins[0] <= data["max_volume"][0])


    ## Node 1 + arc 13 
    m.addConstr(volume_basins[1] == data["vt0"][1] - volume_lt[1] - volume_ltt[1] + volume_lt[0] + data["qin"][1])
    m.addConstr(volume_lt[1] <= (data["horizon"] - data["horizon"]%data["travel_time"][1])*data["channel_flow"][1])
    m.addConstr(volume_ltt[1] <= (data["horizon"]%data["travel_time"][1])*data["channel_flow"][1])
    m.addConstr(volume_lt[1] + volume_ltt[1] <= data["horizon"] * data["channel_flow"][1])
    m.addConstr(volume_basins[1] <= data["max_volume"][1])

    ## Node 2 + arc 23 
    m.addConstr(volume_basins[2] == data["vt0"][2] - volume_lt[2] - volume_ltt[2])
    m.addConstr(volume_lt[2] <= (data["horizon"] - data["horizon"]%data["travel_time"][2])*data["channel_flow"][2])
    m.addConstr(volume_ltt[2] <= (data["horizon"]%data["travel_time"][2])*data["channel_flow"][2])
    m.addConstr(volume_lt[2] + volume_ltt[2] <= data["horizon"] * data["channel_flow"][2])
    m.addConstr(volume_basins[2] <= data["max_volume"][2])

    ## Node 3 + arc 34 
    m.addConstr(volume_basins[3] == data["vt0"][3] - volume_lt[3] - volume_ltt[3] + volume_lt[1] + volume_lt[2] + data["qin"][3])
    m.addConstr(volume_lt[3] <= (data["horizon"] - data["horizon"]%data["travel_time"][3])*data["channel_flow"][3])
    m.addConstr(volume_ltt[3] <= (data["horizon"]%data["travel_time"][3])*data["channel_flow"][3])
    m.addConstr(volume_lt[3] + volume_ltt[3] <= data["horizon"] * data["channel_flow"][3])
    m.addConstr(volume_basins[3] <= data["max_volume"][3])

    ## Node 4 + arc 4out 
    m.addConstr(volume_basins[4] == data["vt0"][4] - volume_lt[4] - volume_ltt[4] + volume_lt[3] + data["qin"][4])
    m.addConstr(volume_lt[4] <= (data["horizon"] - data["horizon"]%data["travel_time"][4])*data["channel_flow"][4])
    m.addConstr(volume_ltt[4] <= data["horizon"]%data["travel_time"][4]*data["channel_flow"][4])
    m.addConstr(volume_lt[4] + volume_ltt[4] <= data["horizon"] * data["channel_flow"][4])
    m.addConstr(volume_basins[4] <= data["max_volume"][4])
    
    m.setObjective(data["weights"][0]*volume_basins[0] + data["weights"][1]*volume_basins[1] + data["weights"][2]*volume_basins[2] + data["weights"][3]*volume_basins[3] + data["weights"][4]*volume_basins[4], GRB.MINIMIZE)
    m.optimize()
    
    for i in data["nodes"]:
        volumes.append(volume_basins[i].x)
        outflows_t.append(volume_lt[i].x)
        outflows_tt.append(volume_ltt[i].x)

    return volumes, outflows_t, outflows_tt

def subproblem(data, volume_lt, volume_ltt, inflows_p1, inflows_p3, inflows_p4):
    # Pond0 -> Generate outflows
    ## Variables
    horizon_t = [i for i in range(0, data["horizon"])]

    m_p0 = Model("SP_0")
    m_p0.Params.OutputFlag = 0
    volume_p0  = m_p0.addVars(horizon_t, lb=0, ub=data["max_volume"][0], vtype=GRB.CONTINUOUS)
    outflow_p0 = m_p0.addVars(horizon_t, lb=0, ub=data["channel_flow"][0] , vtype=GRB.CONTINUOUS)

    # Constarints 

    # Upper limits for flows 
    m_p0.addConstr(outflow_p0.sum([str(i) for i in range(0,  90)]) <= volume_lt[0])
    m_p0.addConstr(outflow_p0.sum([str(i) for i in range(90, 100)]) <= volume_ltt[0])

    # Mass balance 
    # Loss of water is at current time
    m_p0.addConstrs((volume_p0[time] == volume_p0[time-1] - outflow_p0[time] for time in horizon_t[1:]))
    m_p0.addConstr(volume_p0[0] == data["vt0"][0] - outflow_p0[0]) # Inital condition 

    # Physical constraints
    m_p0.addConstr(outflow_p0[0] <= data["vt0"][0]*1.13047751*0.01 +  3.65949363)
    m_p0.addConstrs(outflow_p0[time] <= (volume_p0[time-1]*1.13047751*0.01 +  3.65949363) for time in horizon_t[1:])

    # Objective 
    m_p0.setObjective(outflow_p0.sum(), GRB.MAXIMIZE)
    m_p0.optimize()
    
    outflow_p0 = [outflow_p0[q].x for q in outflow_p0.keys()]
    volume_p0 = [volume_p0[q].x for q in volume_p0.keys()]
    
    # Pond1 -> Generate outflows
    ## Variables 
    m_p1       = Model("SP_1")
    m_p1.Params.OutputFlag = 0
    volume_p1  = m_p1.addVars(horizon_t, lb=0, ub=data["max_volume"][1], vtype=GRB.CONTINUOUS)
    outflow_p1 = m_p1.addVars(horizon_t, lb=0, ub=data["channel_flow"][1] , vtype=GRB.CONTINUOUS)

    # Constarints 
    # Upper limits for flows 
    m_p1.addConstr(outflow_p1.sum([str(i) for i in range(0,60)]) <= volume_lt[1])
    m_p1.addConstr(outflow_p1.sum([str(i) for i in range(60,100)]) <= volume_ltt[1])

    # Mass balance 
    # inflows_p1 carries water from previous time steps that were skipped in the last data["horizon"]. 
    assert(len(inflows_p1) == data["travel_time"][0]) # Test to ensure the length of inflows is good and if the data["horizon"] step is larger, we make it zeros for the rest.

    # This is where the magic happens.
    m_p1.addConstrs((volume_p1[time] == volume_p1[time-1] - outflow_p1[time] + inflows_p1[time] for time in horizon_t[1:data["travel_time"][0]]))
    m_p1.addConstrs((volume_p1[time] == volume_p1[time-1] - outflow_p1[time] + outflow_p0[time-data["travel_time"][0]] for time in horizon_t[data["travel_time"][0]:]))
    m_p1.addConstr(volume_p1[0] == data["vt0"][1] - outflow_p1[0] + inflows_p1[0]) # Inital condition 

    
    # Physical constraints 
    m_p1.addConstr(outflow_p1[0] <= data["vt0"][1]*1.13047751*0.01 +  3.65949363)
    m_p1.addConstrs(outflow_p1[time] <= (volume_p1[time-1]*1.13047751*0.01 +  3.65949363) for time in horizon_t[1:])

    # Objective 
    m_p1.setObjective(outflow_p1.sum(), GRB.MAXIMIZE)
    m_p1.optimize() 

    outflow_p1 = [outflow_p1[q].x for q in outflow_p1.keys()]
    volume_p1 = [volume_p1[q].x for q in volume_p1.keys()]

    # Pond2 -> Generate outflows
    ## Variables 
    m_p2       = Model("SP_2")
    m_p2.Params.OutputFlag = 0

    volume_p2  = m_p2.addVars(horizon_t, lb=0, ub=data["max_volume"][2], vtype=GRB.CONTINUOUS)
    outflow_p2 = m_p2.addVars(horizon_t, lb=0, ub=data["channel_flow"][2] , vtype=GRB.CONTINUOUS)

    # Constarints 
    # Upper limits for flows 
    m_p2.addConstr(outflow_p2.sum([str(i) for i in range(0,50)]) <= volume_lt[2])
    m_p2.addConstr(outflow_p2.sum([str(i) for i in range(50,100)]) <= volume_ltt[2])

    # Mass balance 
    # Loss of water is at current time
    m_p2.addConstrs((volume_p2[time] == volume_p2[time-1] - outflow_p2[time] for time in horizon_t[1:]))
    m_p2.addConstr(volume_p2[0] == data["vt0"][2] - outflow_p2[0]) # Inital condition 

    # Physical constraints
    m_p2.addConstr(outflow_p2[0] <= data["vt0"][2]*1.13047751*0.01 +  3.65949363)
    m_p2.addConstrs(outflow_p2[time] <= (volume_p2[time-1]*1.13047751*0.01 +  3.65949363) for time in horizon_t[1:])

    # Objective 
    m_p2.setObjective(outflow_p2.sum(), GRB.MAXIMIZE)
    m_p2.optimize()

    outflow_p2 = [outflow_p2[q].x for q in outflow_p2.keys()]
    volume_p2 = [volume_p2[q].x for q in volume_p2.keys()]

    # Pond3 -> Generate outflows
    ## Variables 
    m_p3       = Model("SP_3")
    m_p3.Params.OutputFlag = 0

    volume_p3  = m_p3.addVars(horizon_t, lb=0, ub=data["max_volume"][3], vtype=GRB.CONTINUOUS, name="ji")
    outflow_p3 = m_p3.addVars(horizon_t, lb=0, ub=data["channel_flow"][3] , vtype=GRB.CONTINUOUS, name="re")

    # Upper limits for flows 
    m_p3.addConstr(outflow_p3.sum([str(i) for i in range(0, 70)]) <= volume_lt[3]) # For the channel it is sending water to.
    m_p3.addConstr(outflow_p3.sum([str(i) for i in range(70,100)]) <= volume_ltt[3])

    # Mass balance 
    m_p3.addConstrs(volume_p3[time] == volume_p3[time-1] - outflow_p3[time] + inflows_p3[time] for time in horizon_t[1:40]) # assuming travel time 2 is > than tt1
    m_p3.addConstrs(volume_p3[time] == volume_p3[time-1] - outflow_p3[time] + inflows_p3[time] + outflow_p1[time-40] for time in horizon_t[40:50])
    m_p3.addConstrs(volume_p3[time] == volume_p3[time-1] - outflow_p3[time] + outflow_p1[time-40] + outflow_p2[time-50] for time in horizon_t[50:])
    m_p3.addConstr(volume_p3[0] == data["vt0"][3] - outflow_p3[0] + inflows_p3[0]) # Inital condition 

    # Physical constraints
    m_p3.addConstr(outflow_p3[0] <= data["vt0"][3]*1.13047751*0.01 +  3.65949363)
    m_p3.addConstrs(outflow_p3[time] <= (volume_p3[time-1]*1.13047751*0.01 +  3.65949363) for time in horizon_t[1:])

    # Objective 
    m_p3.setObjective(outflow_p3.sum(), GRB.MAXIMIZE)
    m_p3.optimize()
    
    outflow_p3 = [outflow_p3[q].x for q in outflow_p3.keys()]
    volume_p3 = [volume_p3[q].x for q in volume_p3.keys()]

    # Pond4 -> Generate outflows
    ## Variables 
    m_p4       = Model("SP_4")
    m_p4.Params.OutputFlag = 0

    volume_p4  = m_p4.addVars(horizon_t, lb=0, ub=data["max_volume"][4], vtype=GRB.CONTINUOUS)
    outflow_p4 = m_p4.addVars(horizon_t, lb=0, ub=data["channel_flow"][4] , vtype=GRB.CONTINUOUS)

    # Constarints 
    # Upper limits for flows 
    m_p4.addConstr(outflow_p4.sum([str(i) for i in range(0, 90)]) <= volume_lt[4])
    m_p4.addConstr(outflow_p4.sum([str(i) for i in range(90,100)]) <= volume_ltt[4])

    # Mass balance 
    # inflows_p4 carries water from previous time steps that were skipped in the last data["horizon"]. 

    # This is where the magic happens.
    m_p4.addConstrs((volume_p4[time] == volume_p4[time-1] - outflow_p4[time] + inflows_p4[time] for time in horizon_t[1:30]))
    m_p4.addConstrs((volume_p4[time] == volume_p4[time-1] - outflow_p4[time] + outflow_p3[time-data["travel_time"][3]] for time in horizon_t[30:]))
    m_p4.addConstr(volume_p4[0] == data["vt0"][4] - outflow_p4[0] + inflows_p4[0]) # Inital condition 

    
    # Physical constraints 
    m_p4.addConstr(outflow_p4[0] <= data["vt0"][4]*1.13047751*0.01 +  3.65949363)
    m_p4.addConstrs(outflow_p4[time] <= (volume_p4[time-1]*1.13047751*0.01 +  3.65949363) for time in horizon_t[1:])

    # Objective 
    m_p4.setObjective(outflow_p4.sum(), GRB.MAXIMIZE)
    m_p4.optimize()
    
    outflow_p4 = [outflow_p4[q].x for q in outflow_p4.keys()]
    volume_p4 = [volume_p4[q].x for q in volume_p4.keys()]

    # Update the inflows for the next iterations
    inflows_p1 = np.zeros(data["travel_time"][0])
    inflows_p1[:10] = outflow_p0[90:]

    inflows_p3 = copy.deepcopy(np.asarray(outflow_p2[50:]))
    inflows_p3[:40] = inflows_p3[:40] + np.asarray(outflow_p1[60:])

    inflows_p4 = np.zeros(data["travel_time"][3])
    inflows_p4 = outflow_p3[70:]

    # Residual water that would not enter in the data["horizon"]. 

    # Update the initial volume and inflows for master problem.
    data["qin"][1] = sum(outflow_p0[90:])
    data["qin"][3] = sum(outflow_p1[60:]) + sum(outflow_p2[50:])
    data["qin"][4] = sum(outflow_p3[70:])

    data["vt0"][0] = volume_p0[-1] 
    data["vt0"][1] = volume_p1[-1] 
    data["vt0"][2] = volume_p2[-1] 
    data["vt0"][3] = volume_p3[-1]
    data["vt0"][4] = volume_p4[-1]

    print(data["vt0"])
    # Collate data during the horizon 
    outflow_ponds = [outflow_p0, outflow_p1, outflow_p2, outflow_p3, outflow_p4]
    volume_ponds = [volume_p0, volume_p1, volume_p2, volume_p3, volume_p4]
    return data, outflow_ponds, volume_ponds, inflows_p1, inflows_p3, inflows_p4

# Initial Conditions 
INIT_VOLUME = [1000, 800, 1000, 1000, 0]
# Initial Conditions 
#INIT_VOLUME1 = [1000, 500, 500, 1000, 10]

# Physical Constraints
data = {}
data["nodes"] = [0, 1, 2, 3, 4]
data["arcs"] = [0, 1, 2, 3, 4]
data["max_volume"]   = [1000, 1000, 1000, 1200, 1000]
data["channel_flow"] =  [5, 10, 5, 10, 10]
data["travel_time"]  = [10, 40, 50, 30, 10]
data["weights"]  = [1.0, 1.0, 1.0, 0.5, 1.0]
data["vt0"] = INIT_VOLUME
data["qin"] = [0, 0, 0, 0, 0]
data["horizon"] = 100

volume = {}
outflow = {}
for i in data["nodes"]:
    volume[i] = []
for i in data["nodes"]:
    outflow[i] = []

inflows_p1 = np.zeros(data["travel_time"][0])
inflows_p3 = np.zeros(data["travel_time"][2])
inflows_p4 = np.zeros(data["travel_time"][3])

for i in range(0,6):
    volumes, outflows_t, outflows_tt = master(data)
    data, outflow_pond, volume_ponds, inflows_p1, inflows_p3, inflows_p4  = subproblem(data, outflows_t, outflows_tt, inflows_p1, inflows_p3, inflows_p4)
    for i in data["nodes"]:
        for v,o in zip(volume_ponds[i], outflow_pond[i]):
            volume[i].append(v)
            outflow[i].append(o)
print(sum(outflow[4]))

np.save("./volume_100.npy", volume)
np.save("./outflow_100.npy", outflow)

sns.set_style("whitegrid")
# Covnert to numbers

plt.figure(1)
p = 1
for i in data["nodes"]:
    plt.subplot(2,5,p)
    plt.plot(volume[i])
    #plt.plot(volume1[i])
    plt.title(str(i))
    plt.ylim([0, 1200])
    plt.xticks([])
    p = p+1

for i in data["nodes"]:
    plt.subplot(2,5,p)
    plt.ylim([0, 12])
    plt.plot(outflow[i])
    #plt.plot(outflow1[i])
    p+= 1

plt.show()

