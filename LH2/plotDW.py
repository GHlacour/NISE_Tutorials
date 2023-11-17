import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('CG_2DES_doorway_re.dat')

plt.plot(Data[:,0],Data[:,1])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Response function [arb.u.]',fontsize=16)
plt.show()

Data = np.loadtxt('CG_2DES_windows_GB_re.dat')

plt.plot(Data[:,0],Data[:,1])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Response function [arb.u.]',fontsize=16)
plt.show()

Data = np.loadtxt('CG_2DES_windows_SE_re.dat')

plt.plot(Data[:,0],Data[:,1])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Response function [arb.u.]',fontsize=16)
plt.show()

Data = np.loadtxt('CG_2DES_windows_EA_re.dat')

plt.plot(Data[:,0],Data[:,1])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Response function [arb.u.]',fontsize=16)
plt.show()
