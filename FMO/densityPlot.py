import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('AbsoluteDensityMatrix.dat')
vmax=Data.max()
fig = plt.figure()
ax = fig.add_subplot(111)
caxes=ax.matshow(Data,cmap=plt.cm.bwr,vmin=-vmax,vmax=vmax)
ax.set_xlabel('Site',fontsize=16)
ax.set_ylabel('Site',fontsize=16)
ax.xaxis.set_ticks_position('bottom')
fig.colorbar(caxes)
plt.savefig('AbsoluteDensityMatrix.png',dpi=400)
plt.show()

Data = np.loadtxt('SpectralDensityMatrix.dat')
vmax=Data.max()
fig = plt.figure()
ax = fig.add_subplot(111)
caxes=ax.matshow(Data,cmap=plt.cm.bwr,vmin=-vmax,vmax=vmax)
ax.set_xlabel('Site',fontsize=16)
ax.set_ylabel('Site',fontsize=16)
ax.xaxis.set_ticks_position('bottom')
fig.colorbar(caxes)
plt.savefig('SpectralDensityMatrix.png',dpi=400)
plt.show()

