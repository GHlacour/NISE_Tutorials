import numpy as np

steps=100000

# Polarizability matrix xx xy xz yy yz zz components for the two sites
alpha1=np.array([ 1.0,0.0,0.0,1.0,0.0,1.0])
alpha2=np.array([ 1.0,0.0,0.0,1.0,0.0,1.0])

# Dipole vectors x y z for the two sites
mu1=np.array([ 0.0,0.0,1.0])
mu2=np.array([ 0.0,0.0,1.0])

file_alpha=open("Polarizability.bin","wb")
file_mu=open("Dipole.bin","wb")

for st in range(steps):
    # Write timestep number
    step=np.array([st],'float32')
    step.tofile(file_alpha)
    step.tofile(file_mu)
    # Write polarizability for each component for both sites
    for x in range(6):
         alpha=np.array(alpha1[x],'float32')
         alpha.tofile(file_alpha)
         alpha=np.array(alpha2[x],'float32')
         alpha.tofile(file_alpha)
    for x in range(3):
         mu=np.array(mu1[x],'float32')
         mu.tofile(file_mu)
         mu=np.array(mu2[x],'float32')
         mu.tofile(file_mu)

file_alpha.close
file_mu.close
