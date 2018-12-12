import numpy as np
import matplotlib.pyplot as plt
 
for i in range(1,10):
  data = np.loadtxt("result"+str(i)+".csv",delimiter=",")

  plt.hist(data,50)
  plt.savefig("figure"+str(i)+".png");
