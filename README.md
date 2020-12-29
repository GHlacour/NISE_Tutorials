This repository contain examples for NISE spectral simulations with the NISE\_2017 code.

The folder LH2 contain the files needed to construct a Hamiltonian tranjectory for the LH2 light harvesting system of purple bacteria. The parameters in this file originate from the paper: https://doi.org/10.1021/acs.jpcb.0c08126
The Hamiltonian trajectory is created from the 1kzu.pdb protein data bank file. 
1. Use python (for example from Anaconda3) to generate the binary files using the command: python GenHamp.py
2. Submit a linear absorption calculation using the script script1D. You will have to change the location of the NISE\_2017 binary file.
3. Examine the linear response in TD\_Absorption.dat and the absorption spectrum in Absorption.dat for example by plotting with the plot.py script.
4. Analyse the Hamiltonian by submitting the scriptAnalyse script. Examine the output files.
5. Submit a Luminescence calculation using the script scriptLum and plot the output.
6. Submit the 2DUVvis calculation using the submit2D\_MPI script.
7. Find the 2D spectra from the calculated response functions with: 2DFFT input2D
8. Plot the 2D spectra with the plot2D.py script. Make sure that the plotted spectral region correspond to the region of interest.
9. Calculate the circular dichroim by submitting the scriptCD script.




