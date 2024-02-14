import numpy as np

# Create a trajectory with the following number of timesteps
steps=100000

# For NISE it is generally the best to use the z-axis as the unique axis
# for SFG the axis normal to the surface

# Polarizability matrix xx xy xz yy yz zz components for the two sites
alpha1=np.array([ 1.0,0.0,0.0,1.0,0.0,1.0])
alpha2=np.array([ 1.0,0.0,0.0,1.0,0.0,1.0])

# Dipole vectors x y z for the two sites
mu1=np.array([ 0.0,0.0,1.0])
mu2=np.array([ 0.0,0.0,1.0])

# Open the two files
file_alpha=open("Polarizability.bin","wb")
file_mu=open("Dipole.bin","wb")

# Loop over the timesteps and write the dipole and polarizability information
# for each step. (Here we assume the parameters to be fixed, but they could
# fluctuate in the for example due to rotation or wobbling-in-a-cone type
# motion
for st in range(steps):
    # Write timestep number to both files
    step=np.array([st],'float32')
    step.tofile(file_alpha)
    step.tofile(file_mu)
    # Write the transition polarizability for each component for both sites
    for x in range(6):
         alpha=np.array(alpha1[x],'float32')
         alpha.tofile(file_alpha)
         alpha=np.array(alpha2[x],'float32')
         alpha.tofile(file_alpha)
    # Write the transition dipole for each component for both sites
    for x in range(3):
         mu=np.array(mu1[x],'float32')
         mu.tofile(file_mu)
         mu=np.array(mu2[x],'float32')
         mu.tofile(file_mu)

# We are finished - close the two files
file_alpha.close
file_mu.close
