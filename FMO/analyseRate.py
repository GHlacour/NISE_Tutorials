import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Degeneracy=np.array([1,2,2,2])

Data_org=np.loadtxt('RateMatrix.dat')
eigenvalues, eigenvectors = np.linalg.eig(Data_org)

zero_eigenvalue_index = np.where(np.isclose(eigenvalues, 0, atol=1e-5))[0][0]
stationary_distribution_org = eigenvectors[:, zero_eigenvalue_index]
stationary_distribution_org /= np.sum(stationary_distribution_org)  # Normalize to ensure it's a valid probability distribution
# Find effective energy

Ener=-np.log(stationary_distribution_org/Degeneracy)*77/1.44
print(Ener)

Data = np.loadtxt('QC_RateMatrix.dat')
Energies=np.loadtxt('SegmentEnergies.dat',skiprows=1)
Degeneracy=np.array([1,2,2,2])
BZD=Degeneracy*np.exp(-(Energies[:,1])/77*1.44)
BZD=BZD/np.sum(BZD)
print(BZD)

# Eigenvalues and Eigenvectors
eigenvalues, eigenvectors = np.linalg.eig(Data)

zero_eigenvalue_index = np.where(np.isclose(eigenvalues, 0, atol=1e-5))[0][0]
stationary_distribution = eigenvectors[:, zero_eigenvalue_index]
stationary_distribution /= np.sum(stationary_distribution)  # Normalize to ensure it's a valid probability distribution
print(eigenvalues)
# Take care of the phase of the eigenvectors
for i in range(len(eigenvectors[0])):
    if np.isclose(eigenvectors[0, i], 1, atol=1e-5):
        stationary_distribution *= np.sign(eigenvectors[0, i])
        break

# Plotting the Equilibrium Distribution in Percent
states = np.arange(1, len(stationary_distribution) + 1)
percentage_distribution = stationary_distribution * 100
colors = ['red','orange','green','blue']

fig, axes = plt.subplots(4, 1, figsize=(6, 12),sharex=True)
axes[0].bar(states, percentage_distribution,color=colors)
axes[0].plot(states,BZD*100,'x',color='black') 
axes[0].set_ylabel('Eq. Dist. (%)',fontsize=16)
#plt.title('Equilibrium Distribution on FMO Segments')
normal=200/np.sum(np.abs(eigenvectors[:,0]))
axes[1].bar(states, normal*eigenvectors[:,0],color=colors)
axes[1].text(0.75, 0.95, str(round(-eigenvalues[0],2)) + ' ps$^{-1}$', transform=axes[1].transAxes, fontsize=16, verticalalignment='top')
normal=200/np.sum(np.abs(eigenvectors[:,1]))
axes[2].bar(states, normal*eigenvectors[:,1],color=colors)
axes[2].text(0.75, 0.95, str(round(-eigenvalues[1],2)) + ' ps$^{-1}$', transform=axes[2].transAxes, fontsize=16, verticalalignment='top')
normal=200/np.sum(np.abs(eigenvectors[:,3]))
axes[3].bar(states, normal*eigenvectors[:,3],color=colors)
axes[3].text(0.75, 0.95, str(round(-eigenvalues[3],2)) + ' ps$^{-1}$', transform=axes[3].transAxes, fontsize=16, verticalalignment='top')
axes[3].set_xlabel('Segment',fontsize=16)
axes[3].set_xticks([1,2,3,4])
axes[1].set_ylabel('Transfer (%)',fontsize=16)
axes[2].set_ylabel('Transfer (%)',fontsize=16)
axes[3].set_ylabel('Transfer (%)',fontsize=16)
plt.tight_layout()
plt.savefig('FMO_rate.png',dpi=400)
plt.show()
