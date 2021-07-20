import numpy as np
import matplotlib.pyplot as plt
# This python script require vpython to be installed (see vpython.org)
# vpython is used for the structural visualization
from vpython import *
import sys
sys.path.append("../../Hamiltonians/Spectra")
sys.path.append("../../Hamiltonians/Structure")
import Spectra
import Structure
import CD
import visHam

# This file create the Hamiltonian trajectory for LH2 using the structure
# from the 1kzu pdb structure file included.
# To run the code use python 3.7 or newer and type:
# python GenHam.py

E0=12255 # cm-1 The B850 Chromophore average gap
E1=240 # cm-1 The extra energy for the B800 Chromophore gap
MU0=4.481 # Debye The transition dipole before the application of a scale factor
sigma=320 # cm-1 The width of the energy distribution for B850
sigma1=141 # cm-1 The width of the energy distribution for B800
tau=150 # fs The correlation time for the overdamped Brownian oscillators
fcs=1.25; # Factor to scale couplings (resulting in an effective dipole of 5.001 Debye)
dt=3 # fs Time between generated snapshots
units=27 # The number of chromophores
steps=200000 # Number of timesteps
N=units

# Set op constants for generating fluctuations
alpha=np.exp(-dt/tau)
beta=np.sqrt(1-alpha**2)
dsiatoicm=(0.393430307**2)*219474.631370515*(0.529177249**3) # Conv from deb^2/a^3 to cm-1

# Initialize arrays 
#  Positions of Mg, NB, and ND
x=np.zeros((units,3))
y=np.zeros((units,3))
z=np.zeros((units,3))
#  Dipole moment
mu=np.zeros((units,3))
#  Dye position
r=np.zeros((units,3))
#  Hamiltonian
H=np.zeros(int(units*(units+1)/2))
SH=np.zeros(units)
#  Dipole moments
mu4bin=np.zeros((units*3))
#  Positions
pos4bin=np.zeros((units*3))
pos41bin=np.zeros((units*3))
pos42bin=np.zeros((units*3))
# Helping arrays for Hamiltonian construction
B800=np.zeros((units))
ssigma=np.ones((units))*sigma

# Open files
file_input=open("1kzu.pdb","r")

# Symmetry constants to recover C3 symmetry operations
sa=-0.5
sb=0.866025

# Read structure from pdb file
index=0
while True:
    data=file_input.readline()
    if not data: break
    words=data.split()
    if len(words)>2:
       if words[2]=="MG":
          x[index][0]=float(words[6])
          y[index][0]=float(words[7])
          z[index][0]=float(words[8])
          x[index+9][0]=sa*x[index][0]-sb*y[index][0]
          y[index+9][0]=sb*x[index][0]+sa*y[index][0]
          z[index+9][0]=z[index][0]
          x[index+18][0]=sa*x[index][0]+sb*y[index][0]
          y[index+18][0]=-sb*x[index][0]+sa*y[index][0]
          z[index+18][0]=z[index][0]
       if words[2]=="NB":
          x[index][1]=float(words[6])
          y[index][1]=float(words[7])
          z[index][1]=float(words[8])
          x[index+9][1]=sa*x[index][1]-sb*y[index][1]
          y[index+9][1]=sb*x[index][1]+sa*y[index][1]
          z[index+9][1]=z[index][1]
          x[index+18][1]=sa*x[index][1]+sb*y[index][1]
          y[index+18][1]=-sb*x[index][1]+sa*y[index][1]
          z[index+18][1]=z[index][1]
       if words[2]=="ND":
          x[index][2]=float(words[6])
          y[index][2]=float(words[7])
          z[index][2]=float(words[8])
          x[index+9][2]=sa*x[index][2]-sb*y[index][2]
          y[index+9][2]=sb*x[index][2]+sa*y[index][2]
          z[index+9][2]=z[index][2]
          x[index+18][2]=sa*x[index][2]+sb*y[index][2]
          y[index+18][2]=-sb*x[index][2]+sa*y[index][2]
          z[index+18][2]=z[index][2]
          index=index+1

# Verify that the correct number (9) of unique dye atoms were read from file
print(index)

# Construct arrays with transition dipole moments and positions
for atom in range(27):
   # Write positions to human readable file
#   file_x.write(str(x[atom][0]) + " " + str(y[atom][0]) + " " + str(z[atom][0]) + "\n")
   mu[atom][0]=x[atom][2]-x[atom][1]
   mu[atom][1]=y[atom][2]-y[atom][1]
   mu[atom][2]=z[atom][2]-z[atom][1]  
   mum=np.sqrt(mu[atom][0]**2+mu[atom][1]**2+mu[atom][2]**2)
   mu[atom][0]=mu[atom][0]/mum
   mu[atom][1]=mu[atom][1]/mum
   mu[atom][2]=mu[atom][2]/mum
   # Write transition dipole moments to human readable file
#   file_dp.write(str(mu[atom][0]) + " " + str(mu[atom][1]) + " " + str(mu[atom][2]) + "\n")
   r[atom][0]=x[atom][0]
   r[atom][1]=y[atom][0]
   r[atom][2]=z[atom][0]
   # Store data in arrays for saving in binary files
   mu4bin[atom]=mu[atom][0]*MU0*np.sqrt(fcs)
   mu4bin[units+atom]=mu[atom][1]*MU0*np.sqrt(fcs)
   mu4bin[2*units+atom]=mu[atom][2]*MU0*np.sqrt(fcs)
   pos4bin[atom]=r[atom][0]
   pos4bin[units+atom]=r[atom][1]
   pos4bin[2*units+atom]=r[atom][2]
   pos41bin[atom]=x[atom][1]
   pos41bin[units+atom]=y[atom][1]
   pos41bin[2*units+atom]=z[atom][1]
   pos42bin[atom]=x[atom][2]
   pos42bin[units+atom]=y[atom][2]
   pos42bin[2*units+atom]=z[atom][2]

# Generate helping arrays for Hamiltonian construction.
# Accounting for difference between B850 and B800 chromophores
for atom in range(9):
      B800[3*atom+1]=E1
      ssigma[3*atom+1]=sigma1

# Create Hamiltonian first off-diagonal part
for ai in range(27):
   for aj in range(ai):
     if ai!=aj:
        rr=r[ai,:]-r[aj,:]
        rd=np.sqrt(rr[0]**2+rr[1]**2+rr[2]**2)
        # Equation for transition-dipole coupling
        J=sum(mu[ai,:]*mu[aj,:])/(rd**3)-3*sum(mu[ai,:]*rr[:])*sum(mu[aj,:]*rr[:])/(rd**5)
        # Convert to cm-1
        J=J*dsiatoicm*MU0*MU0*fcs
        # Do indexing for tridiagonal matrix
        ind=int(ai+units*aj-(aj*(aj+1)/2))
        # Print information to screen
        print(ai)
        print(aj)
        print(ind)
        H[ind]=J
        print(J)

HH=np.zeros((N,N))
# Construct square Hamiltonian to hunam readable file
for ai in range(27):
  for aj in range(27):
    if aj<ai:
      ind=int(ai+units*aj-(aj*(aj+1)/2))
    if aj>ai:
      ind=int(aj+units*ai-(ai*(ai+1)/2))
    if ai==aj:
       HH[ai,ai]=E0+B800[ai]
#      file_HH.write(str(E0+B800[ai]) + " ")
    if ai!=aj:
       HH[ai,aj]=H[ind]
       HH[aj,ai]=H[ind]
#      file_HH.write(str(H[ind]) + " ")
#  file_HH.write("\n")

# Close all files
file_input.close

# Plot structure
#Structure.visual(r,3*mu,N,1)
Spectra.absorption(HH,mu,N,10)
CD.CD(HH,r,mu,N,10)
visHam.visHam(HH,N)
