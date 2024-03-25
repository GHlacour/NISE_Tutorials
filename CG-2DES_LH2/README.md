Here are the tutorials for the CG_2DES calculation with the LH2 system as an example.
Here, provide one snapshot Hamiltonian and one snapshot dipole file, which you need as an input file to generate the input binary file for the NISE calculation.
First, you need to make sure that your system has installed Python and a required package, numpy.
Second, you just run the two files with Python:  gen_energy_LH2.py and gen_dipole.py, and it generates Energy.txt and Dipole.txt files.
Then you need to convert the txt file to the binary file in order to run the NISE calculation with just run the comment:  ~/NISE/NISE_2017/bin/translate  inpTra(Here you need to change the location to where your NISE program installed)
Now, you can run the MC-FRET calculation to generate the rate matrix by submitting the submit_mcfret.sh file, which will generate the thermal corrected rate matrix without the thermal corrected matrix. You can choose which you want. You need to note that the rate matrix name for the CG-2DES calculation is "RateMatrix.dat." you need to make sure the name is correct for the next calculation.
Finally, you can just submit the job and run the calculation for CG-2DES calculations by submitting the submit_CG2DES.sh file.
