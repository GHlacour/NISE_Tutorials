import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('Absorption.dat')
DataLum = np.loadtxt('Luminescence.dat')
DataLD = np.loadtxt('LD.dat')
DataCD = np.loadtxt('CD.dat')
DataDOS = np.loadtxt('DOS.dat')

fig=plt.figure()
plt.xlim([11500,13000])
plt.plot(Data[:,0],Data[:,1],'k')
plt.plot(DataLum[:,0],DataLum[:,1]*6,'g')
plt.plot(DataLD[:,0],DataLD[:,1],'r')
plt.plot(DataCD[:,0],DataCD[:,1],'b')
plt.plot(DataDOS[:,0],DataDOS[:,1]*10,'k:')
plt.legend(('Absorption','Luminescence','LD','CD','DOS'))
plt.xlabel('$\omega$ [cm$^{-1}$]',fontsize=16)
loc,labels=plt.yticks()
plt.yticks(loc,[])
plt.ylabel('Intensity [arb.u.]',fontsize=16)
plt.axis('tight')
plt.xlim([11000,13000])
plt.show()
fig.savefig('LH2_spectra.eps')

