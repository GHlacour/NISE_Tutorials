This repository contain examples for NISE spectral simulations with the NISE\_2017 code.

### Dimer
The folder Dimer contains the file dimer.py which when executed with *python dimer.py* open a GUI for creating a Hamiltonian trajectory for a dimer system.
1. Use python dimer.py to run the GUI. The default settings creates a homo dimer trajectory when the save button is pressed. 
2. Look at the created Hamiltonian text files: Energy.txt and Dipole.txt. Generate the binary files needed for NISE\_2017 by running *translate inpTra*.
3. Submit the linear absorption calculation using the script script1D. You will have to change the location of the NISE\_2017 binary file.
4. Examine the linear response in TD\_Absorption.dat and the absorption spectrum in Absorption.dat for example by plotting with the plot.py script.
5. Analyse the Hamiltonian by submitting the scriptAnalyse script. Examine the output files. Verify that the information on the average frecuencies and coupling match the input you provided in the GUI.
6. Submit the 2DUVvis calculation using the submit2D\_MPI script.
7. Find the 2D spectra from the calculated response functions with: 2DFFT input2D
8. Plot the 2D spectra with the plot2D.py script. Make sure that the plotted spectral region correspond to the region of interest.

### LH2
The folder LH2 contains the files needed to construct a Hamiltonian tranjectory for the LH2 light harvesting system of purple bacteria. The parameters in this file originate from the paper: https://doi.org/10.1021/acs.jpcb.0c08126
The Hamiltonian trajectory is created from the 1kzu.pdb protein data bank file. 
1. Use python (for example from Anaconda3) to generate the binary files using the command: *python GenHam.py*
2. Submit a linear absorption calculation using the script script1D. You will have to change the location of the NISE\_2017 binary file.
3. Examine the linear response in TD\_Absorption.dat and the absorption spectrum in Absorption.dat for example by plotting with the plot.py script.
4. Analyse the Hamiltonian by submitting the scriptAnalyse script. Examine the output files.
5. Submit a Luminescence calculation using the script scriptLum and plot the output.
6. Submit the 2DUVvis calculation using the submit2D\_MPI script.
7. Find the 2D spectra from the calculated response functions with: 2DFFT input2D
8. Plot the 2D spectra with the plot2D.py script. Make sure that the plotted spectral region correspond to the region of interest.
9. Calculate the circular dichroim by submitting the scriptCD script.

Examples are also given for calculating the couplings on the fly with the transition-dipole (TDC) and extended-dipole (EDC) coupling schemes. This reduce the size of the file needed to store the Hamiltonian considerably as only the diagonal part is stored. The couplings are generated from the dipoles and positions of magnesium atoms in the TDC case and the position of the NB and ND in the EDC case.

### LHCII
The folder LHCII contains the files needed for the LHCII light-harvesting system. The codes are mainly written in MATLAB. A python version is also available, but with more limited functionality.
1. Run *genNISEinput.m* using MATLAB. The input files for NISE will be generated. Parameters can be adjusted in the "Files and Parameters" section.
2. Run NISE calculation for the desired techniques (Abs, Lum, 2D, etc.) using the corresponding input files (input1D, inputLum, input2D, etc.).
3. The code *plotNISE.m* can be used to plot basic spectra.

Alternatively, the script *runNISE.sh* can be used to automatically submit steps 1 and 2 to the cluster (if neccessary, adjust the NISE installation path).

The codes can also be used to calculate for other light-harvesting systems, provided that the corresponding site energies and pdb structure are supplied in the "Energy" and "pdb" folders.

### LH2-CG2DES
Here are the tutorials for the CG_2DES calculation with the LH2 system as an example. Here, provide one snapshot Hamiltonian and one snapshot dipole file, which you need as an input file to generate the input binary file for the NISE calculation. 
1. First, you need to make sure that your system has installed Python and a required package, numpy. 
2. Second, you just run the two files with Python: gen_energy_LH2.py and gen_dipole.py, and it generates Energy.txt and Dipole.txt files.
3. Then you need to convert the txt file to the binary file in order to run the NISE calculation with just run the comment: ~/NISE/NISE_2017/bin/translate inpTra(Here you need to change the location to where your NISE program installed)
4. Now, you can run the MC-FRET calculation to generate the rate matrix by submitting the submit_mcfret.sh file, which will generate the thermal corrected rate matrix without the thermal corrected matrix. You can choose which you want. You need to note that the rate matrix name for the CG-2DES calculation is "RateMatrix.dat." you need to make sure the name is correct for the next calculation.
5. Finally, you can just submit the job and run the calculation for CG-2DES calculations by submitting the submit_CG2DES.sh file.
