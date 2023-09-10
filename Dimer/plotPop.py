import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('Pop.dat')
Rate = np.loadtxt('RateMatrix.dat')

plt.plot(Data[:,0],0.5*np.exp(2*Rate[0,0]*Data[:,0]/1000)+0.5);
plt.plot(Data[:,0],Data[:,1])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Population',fontsize=16)
plt.show()

