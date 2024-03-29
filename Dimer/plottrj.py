import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

H=np.loadtxt('Energy.txt')

plt.plot(H[:,0],H[:,1])
plt.plot(H[:,0],H[:,3])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Site Energy',fontsize=16)
plt.show()
