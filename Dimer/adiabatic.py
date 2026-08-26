import numpy as np
import sys
import os

# Read Hamiltonain Trajectory in Upper Triangular Format
def read_hamiltonian_trajectory(filename, n, num_frames):
    """Reads the Hamiltonian trajectory from a binary file in upper triangular format."""
    hamiltonian_trajectory = []
    with open(filename, 'rb') as f:
        for i in range(num_frames):
            # Skip the frame index
            np.fromfile(f, dtype=np.float32, count=1)
            # Read the upper triangular part of the Hamiltonian matrix
            upper_triangular = np.fromfile(f, dtype=np.float32, count=(n * (n + 1)) // 2)
            # Convert upper triangular to full matrix hermitian square matrix
            hamiltonian = np.zeros((n, n))
            hamiltonian[np.triu_indices(n)] = upper_triangular
            # Mirror upper triangle to lower triangle
            hamiltonian = hamiltonian + np.triu(hamiltonian, k=1).T
            if i == -1:
                print("Read Hamiltonian")
                print(upper_triangular)
                print(hamiltonian)
            hamiltonian_trajectory.append(hamiltonian)
            
    return hamiltonian_trajectory

# Read Dipole Trajectory
def read_dipole_trajectory(filename, n, num_frames):
    """Reads the dipole trajectory from a binary file."""
    dipole_trajectory = []
    with open(filename, 'rb') as f:
        for _ in range(num_frames):
            # Skip the frame index
            np.fromfile(f, dtype=np.float32, count=1)
            # Read the dipole vector (x y z components for each site)
            dipole_vector = np.fromfile(f, dtype=np.float32, count=3*n)
            dipole_trajectory.append(dipole_vector)
    return dipole_trajectory

# Transform to Average Basis
def transform_to_average_basis(hamiltonian_trajectory, dipole_trajectory):
    """Transforms the Hamiltonian and dipole trajectories to the average exciton basis."""
    num_frames = len(hamiltonian_trajectory)
    n = hamiltonian_trajectory[0].shape[0]

    # Initialize arrays for results
    average_exciton_basis_dipole = np.zeros((num_frames, 3*n), dtype=np.float32)
    adiabatic_hamiltonian = np.zeros((num_frames, n, n), dtype=np.float32)
    non_adiabatic_hamiltonian = np.zeros((num_frames, n, n), dtype=np.float32)

    # Find the average Hamiltonian over all frames
    average_hamiltonian = np.mean(hamiltonian_trajectory, axis=0)
    # Find the corresponding eigenvectors of the average Hamiltonian
    _, average_eigenvectors = np.linalg.eigh(average_hamiltonian)
    print("Average Hamiltonian")
    print(average_hamiltonian)
    print(hamiltonian_trajectory[0])

    for i in range(num_frames):
        H = hamiltonian_trajectory[i]
        dipole = dipole_trajectory[i]

        # Transform the dipole vector to the exciton basis
        # We need to keep in mind that the dipole first has all
        # x components, then y components, then z components for each site.
        transformed_dipole = average_eigenvectors.T @ np.transpose(dipole.reshape(3, n))        
        
        # Transform the Hamiltonians
        transformed_hamiltonian = average_eigenvectors.T @ H @ average_eigenvectors

        # Store results keeping in mind that the dipole vector first
        # has all x components, then y components, then z components for each site.
        #average_exciton_basis_dipole[i] = transformed_dipole.flatten()
        transposed_transformed_dipole=transformed_dipole.T
        average_exciton_basis_dipole[i] = transposed_transformed_dipole.flatten()
        # The adiabatic Hamiltonian is the diagonalized form of the Hamiltonian in the exciton
        # but with the off diagonal elements set to zero, while the non-adiabatic Hamiltonian retains the full transformed Hamiltonian.
        adiabatic_hamiltonian[i] = transformed_hamiltonian * np.eye(n)
        non_adiabatic_hamiltonian[i] = transformed_hamiltonian
        if i == 0:
            print(average_eigenvectors)
            print(np.transpose(dipole.reshape(3, n)))
            print(transformed_dipole)

    return {
        'average_exciton_basis_dipole': average_exciton_basis_dipole,
        'adiabatic_hamiltonian': adiabatic_hamiltonian,
        'non_adiabatic_hamiltonian': non_adiabatic_hamiltonian
    }

# Write Binary Hamiltonian File with the original format
def write_binary_file(filename, n, data):
    """Writes the transformed data to a binary file in the specified format."""
    with open(filename, 'wb') as f:
        num_frames = data.shape[0]
        for i in range(num_frames):
            # Write the frame index
            np.array([i], dtype=np.float32).tofile(f)
            # Write the upper triangular part of the Hamiltonian matrix
            upper_triangular = data[i][np.triu_indices(n)]
            upper_triangular.tofile(f)
            # Print first frame to screen
            if i == 0:
                print(upper_triangular)

# Write Binary Dipole File with the original format
def write_binary_dipole_file(filename, n, data):
    """Writes the transformed dipole data to a binary file in the specified format."""
    with open(filename, 'wb') as f:
        num_frames = data.shape[0]
        for i in range(num_frames):
            # Write the frame index
            np.array([i], dtype=np.float32).tofile(f)
            # Write the dipole vector (x y z components for each site)
            data[i].tofile(f)
            if i == 0:
                print(data)


def main(hamiltonian_file, dipole_file, output_prefix, n, num_frames):
    """Main function to process and transform the trajectory files."""
    hamiltonian_trajectory = read_hamiltonian_trajectory(hamiltonian_file, n, num_frames)
    dipole_trajectory = read_dipole_trajectory(dipole_file, n, num_frames)

    results = transform_to_average_basis(hamiltonian_trajectory, dipole_trajectory)

    print("Dipole")
    write_binary_dipole_file(f"{output_prefix}_dipole_avg_basis.bin", n, results['average_exciton_basis_dipole'])
    print("Adiabatic")
    write_binary_file(f"{output_prefix}_adiabatic_hamiltonian.bin", n, results['adiabatic_hamiltonian'])
    print("Nonadiabatic")
    write_binary_file(f"{output_prefix}_nonadiabatic_hamiltonian.bin", n, results['non_adiabatic_hamiltonian'])

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python adiabatic.py <hamiltonian_file> <dipole_file> <output_prefix> <n> <num_frames>")
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))