import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True

DataIn = np.loadtxt('ExcitonStats.dat',skiprows=1)
#Data = DataIn[0:10000000,:]
Data=DataIn

plt.hist2d(Data[:,1],Data[:,2],bins=100,cmap='Greys',norm="log")
plt.xlim([-15,15])
plt.ylim([-15,15])
plt.xlabel(r'$\mu_x$',fontsize=16)
plt.ylabel(r'$\mu_y$',fontsize=16)
plt.axis('square')
plt.show()
