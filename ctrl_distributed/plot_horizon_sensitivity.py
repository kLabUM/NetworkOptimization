import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

data = {}
data["nodes"] = [0, 1, 2, 3, 4]
sns.set_style("whitegrid")
plt.figure(1)

for j in ["60", "100"]:
    volume = np.load("./volume_"+j+".npy").item()
    outflow = np.load("./outflow_"+j+".npy").item()

    p = 1
    for i in data["nodes"]:
        plt.subplot(2,5,p)
        plt.plot(volume[i], label=j)
        plt.title(str(i))
        plt.ylim([0, 1200])
        plt.xticks([])
        p = p+1

    for i in data["nodes"]:
        plt.subplot(2,5,p)
        plt.ylim([0, 12])
        plt.plot(outflow[i], label=j)
        p+= 1

plt.legend()
plt.show()
