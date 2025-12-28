import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm, sinm, cosm, eig
import math
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


################ Decay matrix ################
# Create N+2xN+2 decay matrices
def decayMatrix(N, K_input, KR, KNR):
    K_decay = np.zeros([N + 2, N + 2])
    K_R = np.ones(N) * KR  # If all segments get the same value
    K_NR = np.ones(N) * KNR  # If all segments get the same value
    '''
    #K_R = np.arange(1, N + 1) * 10  # Radiative decay rates
    #K_NR = np.arange(1, N + 1) * 100  # Non-Radiative decay rates (10x bigger)
    '''

    # Copy rate matrix to decay matrix, fill in R/NR decay rates
    K_decay[0:N, 0:N] = K_input
    K_decay[N, 0:N] = K_NR  # Non-Radiative decay rates
    K_decay[N + 1, 0:N] = K_R  # Radiative decay rates (10x bigger)

    # Subtract R/NR decay rates from diagonal
    K_decay[0:N, 0:N] -= np.diag(K_NR)
    K_decay[0:N, 0:N] -= np.diag(K_R)
    return K_NR, K_R, K_decay

################ Double annihilation A0 matrix ################
def annihilationA0(f_states, g_states, KA0):
    K_A0 = np.zeros((g_states, f_states))
    K_A0[0, :] = np.ones(f_states) * KA0  # If all segments get the same value
    '''
    #K_A0[0, :] = np.arange(1, f_states + 1) * 0.1
    '''
    return K_A0

################ Single annihilation A1 matrix ################
def annihilationA1(N, c_w_r, f_states, KA1):
    K_A1 = np.zeros([N, f_states])
    for i in range(f_states):
        if c_w_r[i][0] == c_w_r[i][1]:
            K_A1[c_w_r[i][0]][i] += KA1
    #        for j in range(len(c_w_r)):
    #            K_A1[i,j] += np.sum(c_w_r[j].count(i))
    return K_A1

################ Non-Radiative matrix for each f state ################
# Doubly excited state e.g. 00 goes to 0NR with probability 2*(0->NR)
# Doubly excited state e.g. 01 goes to 0NR with probability (1->NR)
def NRmatrix(N, c_w_r, K_NR):
    K_NRf = np.zeros([N, f_states])
    for i in range(f_states):
        for j in range(2):
            K_NRf[c_w_r[i][j - 1], i] += K_NR[c_w_r[i][j]]
    return K_NRf

################ Radiative matrix for each f state ################
def Rmatrix(N, c_w_r, K_R, f_states):
    K_Rf = np.zeros([N, f_states])
    for i in range(f_states):
        for j in range(2):
            K_Rf[c_w_r[i][j - 1], i] += K_R[c_w_r[i][j]]
    return K_Rf

################ Doubly excited f states matrix ################
def doublyMatrix(c_w_r, K_input):
    mati = np.zeros([len(c_w_r), len(c_w_r)])
    for To in range(mati.shape[0]):
        for From in range(mati.shape[1]):
            if To == From:
                continue
            for i in range(2):
                for j in range(2):
                    if c_w_r[To][i] == c_w_r[From][j]:
                        ind_To, ind_From = c_w_r[To][i - 1], c_w_r[From][j - 1]
                        if not (c_w_r[To][i] == c_w_r[To][i - 1] and i == 1):
                            mati[To, From] += K_input[ind_To, ind_From]
    return mati

################ Puzzle out FD matrix ################
def FDMatrix(f_states, e_states, g_states, K_input, K_NR, K_R, K_A0, K_A1, K_NRf, K_Rf, mati):
    K_FD = np.zeros([f_states + e_states + g_states, f_states + e_states + g_states])
    # Add NR decay from e states+NR (to NR NR)
    K_FD[f_states + e_states, f_states:f_states + N] = K_NR
    # Add NR decay from e states+R
    K_FD[f_states + e_states + 1, f_states + N:f_states + 2 * N] = K_NR
    # Add R decay from e states+NR
    K_FD[f_states + e_states + 1, f_states:f_states + N] = K_R
    # Add R decay from e states+R
    K_FD[f_states + e_states + 2, f_states + N:f_states + 2 * N] = K_R
    # Add e states from transferMatrix
    K_FD[f_states:f_states + N, f_states:f_states + N] = K_input[0:N, 0:N]
    K_FD[f_states + N:f_states + 2 * N, f_states + N:f_states + 2 * N] = K_input[0:N, 0:N]
    # Add double annihilation
    K_FD[f_states + e_states:f_states + e_states + g_states, 0:f_states] += K_A0
    # Add single annihilation
    K_FD[f_states:f_states + e_states // 2, 0:f_states] += K_A1
    # Add non-radiative decay of f states
    K_FD[f_states:f_states + e_states // 2, 0:f_states] += K_NRf
    # Add radiative decay of f states
    K_FD[f_states + e_states // 2:f_states + e_states, 0:f_states] += K_Rf
    # Add f states matrix
    K_FD[0:f_states, 0:f_states] = mati
    # Subtract K_R(NR) from diagonal
    K_FD[f_states:f_states + N, f_states:f_states + N] -= np.diag(K_NR)
    K_FD[f_states + N:f_states + 2 * N, f_states + N:f_states + 2 * N] -= np.diag(K_NR)
    K_FD[f_states:f_states + N, f_states:f_states + N] -= np.diag(K_R)
    K_FD[f_states + N:f_states + 2 * N, f_states + N:f_states + 2 * N] -= np.diag(K_R)
    # Subtract f states from diagonal
    for i in range(f_states):  # Only on the diagonal, subtract all other elements per column
        K_FD[i, i] -= np.sum(K_FD[:, i])
    return K_FD

################ Matrix diagonalization ################
def diagonalize_matrix(matrix):
    # Ensure the input is a numpy array
    matrix = np.array(matrix)

    # Check if the matrix is square
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    # Calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = LA.eig(matrix)

    # Diagonal matrix of eigenvalues
    D = np.diag(eigenvalues)

    # Verify the diagonalization: matrix should be close to eigenvectors * D * inv(eigenvectors)
    A_reconstructed = eigenvectors @ D @ LA.inv(eigenvectors)

    return eigenvalues, eigenvectors, D, A_reconstructed

################ Long waiting times ################
def is_time_infinite(K_FD, eigenvectors, D, tol=1e-10, N_iter=40000, step=100):
    for t in range(10, N_iter, step):
        Product = expm(K_FD * t)
        sumRNR = sum(Product[f_states + e_states:f_states + e_states + g_states, 0])
        print("sumRNR is : ", sumRNR)
        if abs(sumRNR - 1) < tol:
            return t, Product
    print("The search did not converge!")
    print("Please, increase N_iter.")
    print("sumRNR is : ", sumRNR)
    exit(0)


# Given values
system = 'LH2'  # Name of the system
# Current directory (change if necessary by pwd)
path = '/scratch/p317440/Jun25NISE/Tutorials_NISE/F-2DES-CG_LH2/'
#path = '/Users/stephanie/Documents/GitHub/Current_NISE_OD/NISE_Tutorials/0QF/'  # Local
# Generate the filename with the system name or change name to filename
# filename = f'{path}{system}_transfer_matrix.dat'
filename = 'RateMatrix.dat'
# Read the matrix from the .dat file
K_input = np.loadtxt(filename)
N = len(K_input)

# Conversion Factors
ns2ps = 1e3
fs2ps = 1e-3

# Lifetimes from literature
tauF = 986  # 986 ps
tauR = 10*ns2ps  # 10 ns
tauA1 = 0.59  # 0.59 ps

# Compute Rates
kF = 1/tauF
kR = 1/tauR
kA1 = 1/tauA1
kNR = kF-kR  # kF = kR + kNR

# Input values for annihilation and decay
KA0 = 0  # Annihilation rate when no excitons are left behind
KA1 = kA1  # Annihilation rate leaves one exciton behind
KR = kR  # Radiative rate
KNR = kNR  # Non-radiative rate
tol = 1e-6
N_iter = 20000
step = 1000

# Find the segment configurations corresponding to each double excited state
# c_w_r[f_state_number][site_number(0/1)] returns the number of the excited segment in that site number
c_w_r = list(itertools.combinations_with_replacement(range(N), 2))

# All states
f_states = N * (N + 1) // 2
e_states = 2 * N
g_states = 3

# Construct K_FD from sub-matrices
K_NR, K_R, K_decay = decayMatrix(N, K_input, KR, KNR)
K_A0 = annihilationA0(f_states, g_states, KA0)
K_A1 = annihilationA1(N, c_w_r, f_states, KA1)
K_NRf = NRmatrix(N, c_w_r, K_NR)
K_Rf = Rmatrix(N, c_w_r, K_R, f_states)
mati = doublyMatrix(c_w_r, K_input)
K_FD = FDMatrix(f_states, e_states, g_states, K_input, K_NR, K_R, K_A0, K_A1, K_NRf, K_Rf, mati)

########################################################################################################################
eigenvalues, eigenvectors, D, A_reconstructed = diagonalize_matrix(K_FD)
######### Now compute Q1 and Q2 #########

try:
    it, Prod = is_time_infinite(K_FD, eigenvectors, D, tol, N_iter, step)
    roundProd = np.round(Prod, 5)
    Q2 = np.zeros(f_states)
    Q1 = np.zeros(N)
    for i in range(f_states):
        Q2[i] += np.dot(Prod[-3:, i], np.array([0, 1, 0])) + 2 * np.dot(Prod[-3:, i], np.array([0, 0, 1]))
    for i in range(N):
        Q1[i] += np.dot(Prod[-3:, range(f_states, f_states + N)[i]], np.array([0, 1, 0]))
    np.savetxt('Q1s.dat', Q1, delimiter=' ')
    np.savetxt('Q2s.dat', Q2, delimiter=' ')
    f = open('QFactors_Log.txt', 'a')
    f.write(
        f'For system {system} with N={int(N)} segments and rates R={round(KR, 5)}, NR={round(KNR, 5)}, A0={round(KA0, 5)}, and A1={round(KA1, 5)},\n'
        f'Q1:{Q1}\n'f'Q2:{Q2}\n')
    f.close()

except ValueError as err:
    print('ValueError:', err)

plotYN = 0
if plotYN != 0:
    #################### PLOTTING ##################
    factor = 1  # 1000
    # Visualize Matrix as plot #
    fig1, ax1 = plt.subplots(figsize=(2.5, 2.5))
    fig2, ax2 = plt.subplots()
    kmin = np.min(np.min(K_input))
    im1 = ax1.matshow(K_input, cmap='RdBu', vmin=kmin, vmax=-kmin)
    fig1.colorbar(im1, ax=ax1)
    im2 = ax2.matshow(K_FD, cmap='RdBu', vmin=kmin, vmax=-kmin)
    fig2.suptitle(
        'N:' + str(int(N)) + ' with $k_R$:' + str(round(KR, 5)) + ', $k_{NR}$:' + str(round(KNR, 5)) + ', $k_{A0}$:' + str(
            round(KA0, 5)) + ', $k_{A1}:$' + str(round(KA1, 5)))
    fig2.colorbar(im2, ax=ax2)
    print('kmin = ', kmin)
    # Loop over data dimensions and create text annotations.
    '''
    #for i in range(len(K_FD)):
    #    for j in range(len(K_FD)):
    #        text = ax2.text(j, i, int(np.round(K_FD[i, j])), ha="center", va="center", color="k")
    '''
    for i in range(len(K_input)):
        for j in range(len(K_input)):
            text = ax1.text(j, i, int(np.round(factor * K_input[i, j])), ha="center", va="center", color="k")
    filename2 = f'{path}{system}_N{int(N)}-R_{round(KR, 4)}-NR_{round(KNR, 4)}-A0_{round(KA0, 4)}-A1_{round(KA1, 4)}.pdf'
    fig2.savefig(filename2, dpi=400)
    filename1 = f'{path}{system}_N{int(N)}.pdf'
    fig1.savefig(filename1, dpi=400)
    plt.show()

print('Done')

np.set_printoptions(formatter={'float_kind': '{:6.1f}'.format})
print('K_FD')
print(K_FD)
print('t:', it)
print('A1:', KA1, ', A0:', KA0)
print('kNR:', KNR, ', kR:', KR)
print('Q1:', Q1)
print('Q2:', Q2)