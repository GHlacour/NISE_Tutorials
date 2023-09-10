import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('PP.par.dat')
DataP = np.loadtxt('PP.per.dat')

plt.plot(Data[:,0],Data[:,1])
plt.plot(DataP[:,0],DataP[:,1]*3)
plt.xlabel('$\omega$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('Pump Probe [arb.u.]',fontsize=16)
plt.show()

