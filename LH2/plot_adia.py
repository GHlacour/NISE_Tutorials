import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

DataO = np.loadtxt('Absorption_org.dat')
DataN = np.loadtxt('Absorption_non.dat')
DataA = np.loadtxt('Absorption_adia.dat')

plt.plot(DataO[:,0],DataO[:,1],'k',label='ORG')
plt.plot(DataN[:,0],DataN[:,1],'r:',label='NON')
plt.plot(DataA[:,0],DataA[:,1],'b',label='ADIA')
plt.xlabel(r'$\omega$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('Absorption [arb.u.]',fontsize=16)
plt.legend()
plt.show()

