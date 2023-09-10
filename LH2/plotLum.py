import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('Luminescence.dat_org')
Data[:,1]=Data[:,1]/np.max(Data[:,1])
plt.plot(Data[:,0],Data[:,1])
#plt.xlabel('Time [fs]',fontsize=16)
#plt.ylabel('Response function [arb.u.]',fontsize=16)
#plt.show()

Data = np.loadtxt('Luminescence.dat')
Datax = np.loadtxt('Luminescence_x.dat')
Datay = np.loadtxt('Luminescence_y.dat')
Dataz = np.loadtxt('Luminescence_z.dat')
ymax=np.max(Data[:,1])
Data[:,1]=Data[:,1]/ymax
Datax[:,1]=Datax[:,1]/ymax
Datay[:,1]=Datay[:,1]/ymax
Dataz[:,1]=Dataz[:,1]/ymax
plt.plot(Data[:,0],Data[:,1])
plt.plot(Data[:,0],Datax[:,1])
plt.plot(Data[:,0],Datay[:,1])
plt.plot(Data[:,0],Dataz[:,1])
plt.xlabel('$\omega$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('Absorption [arb.u.]',fontsize=16)
plt.show()

