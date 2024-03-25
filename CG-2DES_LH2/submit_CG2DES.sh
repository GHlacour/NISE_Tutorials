#!/bin/bash
#SBATCH --time=1-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=CG_2DES
#SBATCH --mem=4000

module load  foss
~/NISE/NISE_2017/bin/NISE    inputCG2DES  #here you need to change the folder for NISE
