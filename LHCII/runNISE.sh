#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=runNISE
#SBATCH --mem-per-cpu=2000

module load foss
module load MATLAB
matlab -nodisplay < genNISEinput.m
~/git-lacourjansenlab/NISE_2017/bin/NISE input1D
~/git-lacourjansenlab/NISE_2017/bin/NISE inputMCFRET
