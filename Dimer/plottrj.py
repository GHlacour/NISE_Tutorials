import numpy as np
import matplotlib.pyplot as plt

H=np.loadtxt('Energy.txt')

plt.plot(H[:,0],H[:,1])
plt.show()
