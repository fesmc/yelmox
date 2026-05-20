#!/bin/bash 
#SBATCH -p short
#SBATCH -J yelmo-pyr 
#SBATCH -o pyr.out
#SBATCH -e pyr.err
#SBATCH --mem=1800

# Run the job
./yelmox.x yelmo_Pyrenees.nml

# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ydyn="cb_z0=-400 cb_z1=500 beta_u0=100"
