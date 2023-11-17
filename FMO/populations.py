import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

def propagate_population(initial_population, rate_matrix, time_points):
    populations = [initial_population]
    transition_matrix = np.exp(rate_matrix * time_points[-1])
    for i in range(len(time_points) - 1):
        dt = time_points[i + 1] - time_points[i]
        transition_matrix = expm(rate_matrix * dt)
        population = np.dot(transition_matrix, populations[-1])
        populations.append(population)
    return np.array(populations)

def analyze_rate_matrix(rate_matrix, tol=1e-6):
    eigenvalues, left_eigenvectors = np.linalg.eig(rate_matrix)
    left_eigenvectors_inv = np.linalg.inv(left_eigenvectors)
    equilibrium_vector = left_eigenvectors[:, np.abs(eigenvalues) < tol]
    equilibrium_distribution = equilibrium_vector / np.sum(equilibrium_vector)
    return eigenvalues, left_eigenvectors, equilibrium_distribution

# Example rate matrix and initial population
rate_matrix=np.loadtxt("QC_RateMatrix.dat")
#rate_matrix=np.loadtxt("RateMatrix.dat")

initial_population=np.array([0.0, 0.0, 0.0,1.0])

# Time points
time_points = np.linspace(0, 2, 100)

# Propagate population for the given time points
populations = propagate_population(initial_population, rate_matrix, time_points)

# Plot the populations
plt.figure()
plt.plot(time_points, populations[:, 0], label='Segment 1')
plt.plot(time_points, populations[:, 1], label='Segment 2')
plt.plot(time_points, populations[:, 2], label='Segment 3')
plt.plot(time_points, populations[:, 3], label='Segment 4')
plt.xlabel('Time [ps]')
plt.ylabel('Population')
plt.legend()
plt.show()

eigenvalues, eigenvectors,equilibrium_distribution = analyze_rate_matrix(rate_matrix)
print(np.abs(eigenvalues) < 1e-6)
print(eigenvectors[:, np.abs(eigenvalues) < 1e-6])
# Print the results
print("Eigenvalues:")
print(eigenvalues)
print()
print("Eigenvectors:")
print(eigenvectors)
print()
print("Equilibrium Distribution:")
print(equilibrium_distribution)

