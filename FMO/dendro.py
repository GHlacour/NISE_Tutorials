import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

plt.rcParams.update({
    'font.size': 14,          # Default font size
    'axes.titlesize': 16,     # Title font size
    'axes.labelsize': 14,     # Axes label font size
    'xtick.labelsize': 12,    # X-axis tick label font size
    'ytick.labelsize': 12,    # Y-axis tick label font size
    'legend.fontsize': 12,    # Legend font size
    'figure.titlesize': 18    # Figure title font size
})

Data = np.loadtxt('AbsoluteDensityMatrix.dat')
Data=Data/Data[0,0]
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

# Perform hierarchical clustering and obtain linkage matrix
print(Data)
pData=-np.log(Data)
pData=pData-np.diag(np.diag(pData))
print(pData)
DataS=squareform(pData)
print(DataS)
cut=-np.log(0.67)
print(cut)
linkage_matrix = linkage(DataS, method='complete')
linkage_matrix = linkage(DataS, method='single')
#linkage_matrix = linkage(Data, method='average')

# Plot the dendrogram
dendrogram(linkage_matrix,color_threshold=cut, labels=[f'Site {i+1}' for i in range(len(Data))])
#plt.plot([0,len(Data)],[cut,cut],'k--')
plt.axhline(c='k',linestyle='--',y=cut)
#plt.title('Hierarchical Clustering Dendrogram')
#plt.xlabel('Sites')
plt.ylabel('Segment Distance')
plt.savefig('FMO_Dendrogram.png',dpi=400)
plt.show()

