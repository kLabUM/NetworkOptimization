from gurobipy import *
import numpy  as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")

def to_numpy(gurobi_variable):
    temp = [gurobi_variable[i].x for i in gurobi_variable.keys()]
    return np.asarray(temp)

# Common Containers and variables 
horizon  = [time for time in range(0, 20)]  # Time steps

# Pond1 -> Generate outflows
## Variables 
m_p1       = Model("Pond1")
volume_p1  = m_p1.addVars(horizon, lb=0, ub=100, vtype=GRB.CONTINUOUS, name="Volume")
outflow_p1 = m_p1.addVars(horizon, lb=0, ub=10 , vtype=GRB.CONTINUOUS, name="Outflows")

# Constarints 
m_p1.addConstr(volume_p1[0] == 90, "Inital Volume")
m_p1.addConstr(outflow_p1.sum() <= 100, "Total Mass Movement")
m_p1.addConstrs((volume_p1[time] == volume_p1[time-1] - outflow_p1[time-1] for time in horizon[1:]), name="Mass balance")
m_p1.addConstrs(outflow_p1[time] <= (volume_p1[time-1]*6.14647148e-02 + 3.81)*0.447 for time in horizon[1:])

# Objective 
m_p1.setObjective(volume_p1.sum(), GRB.MINIMIZE)
m_p1.optimize()


# Pond2 -> Generate Outflows 
## Variables 
m_p2       = Model("Pond2")
volume_p2  = m_p2.addVars(horizon, lb=0, ub=100, vtype=GRB.CONTINUOUS, name="Volume")
outflow_p2 = m_p2.addVars(horizon, lb=0, ub=10 , vtype=GRB.CONTINUOUS, name="Outflows")

# Constarints 
m_p2.addConstr(volume_p2[0] == 90, "Inital Volume")
m_p2.addConstr(outflow_p2.sum() <= 100, "Total Mass Movement")
m_p2.addConstrs((volume_p2[time] == volume_p2[time-1] - outflow_p2[time-1] for time in horizon[1:]), name="Mass balance")
m_p2.addConstrs(outflow_p2[time] <= (volume_p2[time-1]*6.14647148e-02 + 3.81)*0.447 for time in horizon[1:])

# Objective 
m_p2.setObjective(volume_p2.sum(), GRB.MINIMIZE)
m_p2.optimize()


# Pond3 -> Take in the outflows and plan how you want to release them.
## Generate the outflows from the above into the downstream 
inflow_p1 = to_numpy(outflow_p1)
inflow_p2 = to_numpy(outflow_p2)
inflows = np.zeros(20)
inflows[10:] = inflows[10:] + inflow_p1[:10] + inflow_p2[:10] 

m_p3       = Model("Pond3")
volume_p3  = m_p3.addVars(horizon, lb=0, ub=100, vtype=GRB.CONTINUOUS, name="Volume")
outflow_p3 = m_p3.addVars(horizon, lb=0, ub=10 , vtype=GRB.CONTINUOUS, name="Outflows")
# Constarints 
m_p3.addConstr(volume_p3[0] == 90 + inflows[0], "Inital Volume")
m_p3.addConstr(outflow_p3.sum() <= 100, "Total Mass Movement")
m_p3.addConstrs((volume_p3[time] == volume_p3[time-1] - outflow_p3[time-1] + inflows[time] for time in horizon[1:]), name="Mass balance")
m_p3.addConstrs(outflow_p3[time] <= (volume_p3[time-1]*6.14647148e-02 + 3.81)*0.447 for time in horizon[1:])

# Objective 
m_p3.setObjective(volume_p3.sum(), GRB.MINIMIZE)
m_p3.optimize()

# Moment of truth 
plt.subplot(2,3,1)
plt.plot(to_numpy(volume_p1))
plt.ylabel("Volume")
plt.subplot(2,3,2)
plt.plot(to_numpy(volume_p2))
plt.subplot(2,3,3)
plt.plot(to_numpy(volume_p3))

plt.subplot(2,3,4)
plt.plot(to_numpy(outflow_p1))
plt.ylabel("Outflows")
plt.subplot(2,3,5)
plt.plot(to_numpy(outflow_p2))
plt.subplot(2,3,6)
plt.plot(to_numpy(outflow_p3))
plt.show()
