import numpy as np
import matplotlib.pyplot as plt
 
for i in range(1,10):
  data = np.loadtxt("result_cut"+str(i)+".csv",delimiter=",")

  plt.clf()
  plt.hist(data,50)
  plt.savefig("figure_cut"+str(i)+".png");

for i in range(1,10):
  data = np.loadtxt("result_shf"+str(i)+".csv",delimiter=",")

  plt.clf()
  plt.hist(data,50)
  plt.savefig("figure_shf"+str(i)+".png");
