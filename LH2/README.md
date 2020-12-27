# LH2_Hamiltonian
The GenHam.py file can be used to generate a Hamiltonain trajectory for the excitonic states of LH2 using the 1kzu.pdb file to define the structure. Run the python script with python 3.7 using the command:\
python GenHam.py\
\
Input file:\
1kzu.pdb\
\
Binary output files:\
Energy.bin - The Hamiltonian Trajectory\
Dipole.bin - The Transition dipole Trajectory\
Positions.bin - The Position Trajectory\
\
Human readable output files:\
Ham.txt - The average Hamiltonian\
Dipole.txt - The transition dipoles\
Positions.txt - Dye center positions (at Magnesium)\
\
The binary output files can be used as input for the NISE spectral calculation programme, which can be found at: https://github.com/GHlacour/NISE_2017

Please, refer to the paper "Observation of ultrafast coherence transfer and degenerate states with polarization controlled two-dimensional electronic spectroscopy" by Andy S. Sardjan, Floris P. Westerman, Jennifer P. Ogilvie, and Thomas L.C. Jansen (https://pubs.acs.org/doi/10.1021/acs.jpcb.0c08126 J. Phys. Chem. B Tanimura Feestschrift (2020)) if you use this code.
