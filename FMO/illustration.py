import numpy as np
import matplotlib.pyplot as plt

# Energy levels of the four states
#energies = [4, 3, 2, 1]
energies= [ 12140.389648, 12424.981445, 12617.533203, 12524.882812]

# Transition probabilities between states
transition_matrix = np.array([
    [0.0, 7.0, 0.38, 0.0],
    [0.038,0.0, 0.1, 2.5],
    [0.0, 0.007, 0.0, 0.22],
    [0.0, 4.0, 0.65, 0.0]
])
#transition_matrix=np.transpose(transition_matrix)

# Compute relaxation rates based on transition probabilities
relaxation_rates = np.max(transition_matrix, axis=1)

# Plot the energy relaxation
plt.figure(figsize=(8, 6))
plt.title('Energy Relaxation')
plt.xlabel('States')
plt.ylabel('Energy Level')

# Plot energy levels as horizontal lines
for i, energy in enumerate(energies):
    plt.hlines(energy, 0, 5*len(energies)-1, color='gray', linestyle='-', alpha=0.5)

# Plot state transitions as vertical arrows
for i in range(len(energies)):
    for j in range(len(energies)):
        if i != j:
            arrow_length = energies[i] - energies[j]
            arrow_direction = np.sign(arrow_length)
            arrow_length = abs(arrow_length) - 10
            arrow_head_length = min(10, arrow_length)
            plt.arrow(5*i+j, energies[j], 0, arrow_direction * arrow_length,
                      color='blue', alpha=0.7, width=0.02 + 0.1 * transition_matrix[i,j],
                      head_length=arrow_head_length)

#plt.xticks(range(len(energies)), [f'State {i+1}' for i in range(len(energies))])
plt.tight_layout()
plt.show()

