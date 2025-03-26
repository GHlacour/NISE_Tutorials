import numpy as np
from tqdm import tqdm

print("")
print("Generating the binary trajectroy for the FMO system.")
print("The average Hamiltonian is take from:")
print("Brixner, T.; Stenger, J.; Vaswani, H. M.; Cho, M.; Blankenship, R. E.;")
print("Fleming, G. R. Two-Dimensional Spectroscopy of Electronic Couplings in")
print("Photosynthesis. Nature 2005, 434, 625−628.")
print("The frequency fluctuations are adapted from:")
print("Tempelaar, R.; Jansen, T.L.C.; Knoester, J. Vibrational Beatings")
print("Conceal Evidence of Electronic Coherence in the FMO Light-Harvesting")
print("Complex. J. Phys. Chem. B. 2014, 118, 12865-12872.")
print("") 
# Run this python script to generate the Hamiltonian trajectories
units=7
Horg=np.array([12420, -106, 8, -5, 6, -8, -4,  -106, 12560, 28, 6, 2, 13, 1, 8, 28, 12140, -62, -1, -9, 17, -5, 6, -62, 12315, -70, -19, -57, 6, 2, -1, -70, 12460, 40, -2, -8, 13, -9, -19, 40, 12500, 32, -4, 1, 17, -57, -2, 32, 12400])
H=Horg.reshape(units,units)

mu=np.array([0.74, 0.86, 0.20, 0.80, 0.74, 0.14, 0.50,  0.56, -0.50, -0.96, 0.53, -0.66, 0.88, 0.71,  0.37, 0.11, 0.21, 0.28, -0.16, -0.46, 0.50])

#  Hamiltonian
HH=np.zeros(int(units*(units+1)/2))
Hc=np.zeros(units*units)
#  Dipole moments
mu4bin=np.zeros((units*3))

# Find eigenbasis
E,c=np.linalg.eigh(H)
sigma=75 # Disorder at 77 K
tau=140 # Timescale for disorder from Tempelaar paper

dt=5
steps=1000000
alpha=np.exp(-dt/tau)
beta=np.sqrt(1-alpha**2)

file_H=open("Energy.bin","wb")
file_mu=open("Dipole.bin","wb")

# Copy Coupling
for ai in range(units):
    for aj in range(ai):
        ind=int(ai+units*aj-(aj*(aj+1)/2))
        HH[ind]=Horg[ai*units+aj]

# Copy Dipoles
#for ai in range(units):
#    for x in range(3):
#        mu4bin[x*units+ai]=mu[3*ai+x];
mu4bin=mu

# Create initial random numbers
diag=np.random.randn(units)
# Create diagonal elements
for st in tqdm(range(steps)):
    for ai in range(units):
        ind=int(ai+units*ai-(ai*(ai+1)/2))
        HH[ind]=diag[ai]*sigma+Horg[ai*units+ai]
        # Update random numbers according to J. Chem. Phys. 127:084507 (2007)
    diag=diag*alpha+np.random.randn(units)*beta
    # Save Hamiltonian, and dipole to binary files
    Hf=np.array(HH,'float32')
    step=np.array([0],'float32')
    step.tofile(file_H)
    Hf.tofile(file_H)

    step.tofile(file_mu)
    muf=np.array(mu4bin,'float32')
    muf.tofile(file_mu)

file_H.close
file_mu.close
