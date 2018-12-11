import numpy as np
import matplotlib.pyplot as plt
import sys
 
args = sys.argv
data = np.loadtxt(args[1],delimiter=",")

plt.hist(data)
plt.show()
