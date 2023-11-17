import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('CorrelationMatrix.dat')

plt.plot(Data[:,0],Data[:,1])
plt.plot(Data[:,0],Data[:,2])
plt.plot(Data[:,0],Data[:,3])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Correlation function [cm $^{-2}$]',fontsize=16)
plt.show()


