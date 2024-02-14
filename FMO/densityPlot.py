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
# Increment all x-axis tick labels by 1
#xticks = [tick + 1 for tick in ax.get_xticks()]
#ax.set_xticks(xticks)

# Increment all y-axis tick labels by 1
#yticks = [tick + 1 for tick in ax.get_yticks()]
#ax.set_yticks(yticks)

# Update the x-axis and y-axis labels
ax.set_xticklabels([int(label + 1) for label in ax.get_xticks()])
ax.set_yticklabels([int(label + 1) for label in ax.get_yticks()])
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

