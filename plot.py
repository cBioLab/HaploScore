import numpy as np
import matplotlib.pyplot as plt
 
data = np.loadtxt("result.csv",delimiter=",")

plt.hist(data)
plt.show()