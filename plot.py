import numpy as np
import matplotlib.pyplot as plt
 
data = np.loadtxt("results.csv",delimiter=",")

plt.hist(data)
plt.show()