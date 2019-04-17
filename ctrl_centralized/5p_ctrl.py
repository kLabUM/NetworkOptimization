from gurobipy import *
import numpy  as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
# Create a model 
m = Model("5Ponds")

# Physical Constraints
nodes = [0, 1, 2, 3, 4]
arcs = [0, 1, 2, 3, 4]
max_volume   = [1000, 1000, 1000, 1000, 1000]
channel_flow = [5, 10, 5, 10, 10]
travel_time  = [10, 40, 50, 30, 10]
weights      = [1.0, 1.0, 1.0, 0.5, 1.0]

analysis_horizon = 600 
temp_horizon = [i for i in range(0, analysis_horizon)]

# Initial Conditions 
INIT_VOLUME = [1000, 800, 1000, 1000, 0]

# Create variables for optimization 
volume_nodes = m.addVars(nodes,
        temp_horizon,
        lb=0.0,
        name="volume")

flow_arcs = m.addVars(nodes,
        temp_horizon,
        lb=0.0,
        name="flow")

# Create the constraints
# 1. Set the upper bound for the volume in nodes
m.addConstrs(volume_nodes[node, time] <= max_volume[node] for node in nodes for time in temp_horizon)

# 2. Set the upper bound for flow in the arcs
m.addConstrs(flow_arcs[arc, time] <= channel_flow[arc] for arc in arcs for time in temp_horizon)

# 3. Set the inital water level 
m.addConstrs(volume_nodes[node, 0] == INIT_VOLUME[node] for node in nodes)

# 4. Set the discharge constraint on the outflows
m.addConstrs(flow_arcs[arc, time] <= (volume_nodes[arc, time-1]*1.13047751*0.01  + 3.65949363) for arc in arcs for time in temp_horizon[1:])

# 5. Mass balance constraint for each node
## Node 0
m.addConstrs(volume_nodes[0, time] == volume_nodes[0, time-1] - flow_arcs[0, time] for time in temp_horizon[1:])

## Node 1
m.addConstrs(volume_nodes[1, time] == volume_nodes[1, time-1] - flow_arcs[1, time] for time in temp_horizon[1:travel_time[0]])
m.addConstrs(volume_nodes[1, time] == volume_nodes[1, time-1] - flow_arcs[1, time] + flow_arcs[0, time-travel_time[0]] for time in temp_horizon[travel_time[0]:])

## Node 2
m.addConstrs(volume_nodes[2, time] == volume_nodes[2, time-1] - flow_arcs[2, time] for time in temp_horizon[1:])

## Node 3
m.addConstrs(volume_nodes[3, time] == volume_nodes[3, time-1] - flow_arcs[3, time] for time in temp_horizon[1:travel_time[1]])
m.addConstrs(volume_nodes[3, time] == volume_nodes[3, time-1] - flow_arcs[3, time] + flow_arcs[1, time-travel_time[1]] for time in temp_horizon[travel_time[1]:travel_time[2]])
m.addConstrs(volume_nodes[3, time] == volume_nodes[3, time-1] - flow_arcs[3, time] + flow_arcs[1, time-travel_time[1]] + flow_arcs[2, time-travel_time[2]] for time in temp_horizon[travel_time[2]:])

## Node 4
m.addConstrs(volume_nodes[4, time] == volume_nodes[4, time-1] - flow_arcs[4, time] for time in temp_horizon[1:travel_time[3]])
m.addConstrs(volume_nodes[4, time] == volume_nodes[4, time-1] - flow_arcs[4, time] + flow_arcs[3, time-travel_time[3]] for time in temp_horizon[travel_time[3]:])

# Weighted objective function
# TODO : figure out how to loop the objective. 
m.setObjective(quicksum(volume_nodes[0, time] * weights[0] for time in temp_horizon) + quicksum(volume_nodes[1, time] * weights[1] for time in temp_horizon) + quicksum(volume_nodes[2, time] * weights[2] for time in temp_horizon)+ quicksum(volume_nodes[3, time] * weights[3] for time in temp_horizon)  + quicksum(volume_nodes[4, time] * weights[4] for time in temp_horizon) , GRB.MINIMIZE)

m.optimize()

sns.set_style("whitegrid")#, {'axes.spines.right': False, 'axes.spines.top': False, 'ytick.left': True, 'xtick.bottom': True})
# Covnert to numbers
volume = {}
for i in nodes:
    volume[i] = []
    for j in temp_horizon:
        volume[i].append(volume_nodes[(i, j)].x)

outflow = {}
for i in arcs:
    outflow[i] = []
    for j in temp_horizon:
        outflow[i].append(flow_arcs[(i, j)].x)

plt.figure(1)
p = 1
for i in nodes:
    plt.subplot(2,5,p)
    plt.plot(volume[i])
    plt.title(str(i))
    plt.ylim([0, 1200])
    plt.xticks([])
    p = p+1

for i in nodes:
    plt.subplot(2,5,p)
    plt.plot(outflow[i])
    plt.ylim([0, 12])
    p+= 1

plt.show()
