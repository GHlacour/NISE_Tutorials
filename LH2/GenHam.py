import numpy
from tqdm import tqdm

# This file create the Hamiltonian trajectory for LH2 using the structure
# from the 1kzu pdb structure file included.
# To run the code use python 3.7 or newer and type:
# python GenHam.py

# Note that in this file the B850 alpha and beta are split by 300 cm-1
# to improve the CD spectra. The sigma value for these chromopores is
# also reduced from 320 to 305 cm-1 compared to the study by Sardjan et al.

print("")
print("Generating the binary trajectroy for the FMO system.")
print("The average Hamiltonian adapted from:")
print("Sardjan, A. S. et al. J. Phys. Chem. B 2020, 124, 9420-9427.")
print("Note that here the B850 alpha and beta are split by 300 cm-1")
print("to improve the CD spectra. The sigma value for these chromopores is")
print("also reduced from 320 to 305 cm-1 compared to the study by Sardjan et al.")
print("")

E0=12255 # cm-1 The B850 Chromophore average gap
E1=240 # cm-1 The extra energy for the B800 Chromophore gap
E2=-300 # Gap between alpha and beta chromophores 
MU0=4.481 # Debye The transition dipole before the application of a scale factor
sigma=320 # cm-1 The width of the energy distribution for B850
sigma1=141 # cm-1 The width of the energy distribution for B800
sigma=305 # Adapted for CD
tau=150 # fs The correlation time for the overdamped Brownian oscillators
fcs=1.25 # Factor to scale couplings (resulting in an effective dipole of 5.001 Debye)

dt=3 # fs Time between generated snapshots
units=27 # The number of chromophores
steps=200000 # Number of timesteps

# Set op constants for generating fluctuations
alpha=numpy.exp(-dt/tau)
beta=numpy.sqrt(1-alpha**2)
dsiatoicm=(0.393430307**2)*219474.631370515*(0.529177249**3) # Conv from deb^2/a^3 to cm-1

# Initialize arrays 
#  Positions of Mg, NB, and ND
x=numpy.zeros((units,3))
y=numpy.zeros((units,3))
z=numpy.zeros((units,3))
#  Dipole moment
mu=numpy.zeros((units,3))
#  Dye position
r=numpy.zeros((units,3))
#  Hamiltonian
H=numpy.zeros(int(units*(units+1)/2))
SH=numpy.zeros(units)
#  Dipole moments
mu4bin=numpy.zeros((units*3))
#  Positions
pos4bin=numpy.zeros((units*3))
pos41bin=numpy.zeros((units*3))
pos42bin=numpy.zeros((units*3))
# Helping arrays for Hamiltonian construction
B800=numpy.zeros((units))
ssigma=numpy.ones((units))*sigma

# Open files
file_input=open("1kzu.pdb","r")
file_H=open("Energy.bin","wb")
file_SH=open("SiteEnergy.bin","wb")
file_x=open("Positions.txt","w")
file_mu=open("Dipole.bin","wb")
file_pos=open("Positions.bin","wb")
file_pos2=open("Positions2.bin","wb")
file_HH=open("Ham.txt","w")
file_dp=open("Dipole.txt","w")

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
#print(index)

# Generate helping arrays for Hamiltonian construction.
# Accounting for difference between B850 and B800 chromophores
for atom in range(9):
      B800[3*atom+1]=E1
      ssigma[3*atom+1]=sigma1
      B800[3*atom+2]=+E2/2
      B800[3*atom]=-E2/2

# Construct arrays with transition dipole moments and positions
for atom in range(27):
   # Write positions to human readable file
   file_x.write(str(x[atom][0]) + " " + str(y[atom][0]) + " " + str(z[atom][0]) + "\n")
   mu[atom][0]=x[atom][2]-x[atom][1]
   mu[atom][1]=y[atom][2]-y[atom][1]
   mu[atom][2]=z[atom][2]-z[atom][1]   
   mum=numpy.sqrt(mu[atom][0]**2+mu[atom][1]**2+mu[atom][2]**2)
   mu[atom][0]=mu[atom][0]/mum
   mu[atom][1]=mu[atom][1]/mum
   mu[atom][2]=mu[atom][2]/mum
   # Write transition dipole moments to human readable file
   file_dp.write(str(mu[atom][0]) + " " + str(mu[atom][1]) + " " + str(mu[atom][2]) + "\n")
   r[atom][0]=x[atom][0]
   r[atom][1]=y[atom][0]
   r[atom][2]=z[atom][0]
   # Store data in arrays for saving in binary files
   mu4bin[atom]=mu[atom][0]*MU0*numpy.sqrt(fcs)
   mu4bin[units+atom]=mu[atom][1]*MU0*numpy.sqrt(fcs)
   mu4bin[2*units+atom]=mu[atom][2]*MU0*numpy.sqrt(fcs)
   pos4bin[atom]=r[atom][0]
   pos4bin[units+atom]=r[atom][1]
   pos4bin[2*units+atom]=r[atom][2]
   pos41bin[atom]=x[atom][1]
   pos41bin[units+atom]=y[atom][1]
   pos41bin[2*units+atom]=z[atom][1]
   pos42bin[atom]=x[atom][2]
   pos42bin[units+atom]=y[atom][2]
   pos42bin[2*units+atom]=z[atom][2]

# Create Hamiltonian first off-diagonal part
for ai in range(27):
   for aj in range(ai):
     if ai!=aj:
        rr=r[ai,:]-r[aj,:]
        rd=numpy.sqrt(rr[0]**2+rr[1]**2+rr[2]**2)
        # Equation for transition-dipole coupling
        J=sum(mu[ai,:]*mu[aj,:])/(rd**3)-3*sum(mu[ai,:]*rr[:])*sum(mu[aj,:]*rr[:])/(rd**5)
        # Convert to cm-1
        J=J*dsiatoicm*MU0*MU0*fcs
        # Do indexing for tridiagonal matrix
        ind=int(ai+units*aj-(aj*(aj+1)/2))
        # Print information to screen
#        print(ai)
#        print(aj)
#        print(ind)
        H[ind]=J
#        print(J)

# Create initial random numbers
diag=numpy.random.randn(27)
# Create diagonal elements
for st in tqdm(range(steps)):
   for ai in range(27):
      ind=int(ai+units*ai-(ai*(ai+1)/2))
      # Find energy gap including shift for B800 chromophores
      H[ind]=diag[ai]*ssigma[ai]+E0+B800[ai]
      SH[ai]=H[ind]
   # Update random numbers according to J. Chem. Phys. 127:084507 (2007)   
   diag=diag*alpha+numpy.random.randn(27)*beta
   # Save Hamiltonian, dipoles, and positions to binary files
   # Full Hamiltonian
   Hf=numpy.array(H,'float32')
   step=numpy.array([0],'float32')
   step.tofile(file_H)  
   Hf.tofile(file_H)
   # Diagonal part of Hamiltonian
   HSf=numpy.array(SH,'float32')
   step=numpy.array([0],'float32')
   step.tofile(file_SH)
   HSf.tofile(file_SH)
   # Dipole moments
   step.tofile(file_mu)
   muf=numpy.array(mu4bin,'float32')
   muf.tofile(file_mu)
   # Positions of Mg
   step=numpy.array([100],'float32') # Box size in Angstrom
   step.tofile(file_pos)
   puf=numpy.array(pos4bin,'float32')
   puf.tofile(file_pos)
   # Positions of NB and ND
   step.tofile(file_pos2)
   puf2=numpy.array(pos41bin,'float32')
   puf2.tofile(file_pos2)
   step.tofile(file_pos2)
   puf2=numpy.array(pos42bin,'float32')
   puf2.tofile(file_pos2)

# Write square Hamiltonian to human readable file
for ai in range(27):
  for aj in range(27):
    if aj<ai:
      ind=int(ai+units*aj-(aj*(aj+1)/2))
    if aj>ai:
      ind=int(aj+units*ai-(ai*(ai+1)/2))
    if ai==aj:
      file_HH.write(str(E0+B800[ai]) + " ")
    if ai!=aj:
      file_HH.write(str(H[ind]) + " ")
  file_HH.write("\n")

# Close all files
file_input.close
file_H.close
file_SH.close
file_x.close
file_mu.close
file_pos.close
file_pos2.close
file_dp.close
