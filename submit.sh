#!/bin/bash 
#SBATCH -p long
#SBATCH -J yelmo 
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#SBATCH --mem=18000
#SBATCH -t 15-02:30:15

##SBATCH -p long_nolimits
##SBATCH -q qnolimits
##SBATCH -J ismip6_AIS 
##SBATCH -o yelmo.out
##SBATCH -e yelmo.err
##SBATCH -t 14-23:30:15

# Run the job
./yelmox_ismip6.x yelmo_ismip6_Antarctica_spinup.nml
#./yelmox_ismip6.x yelmo_ismip6_Antarctica.nml



# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ydyn="cb_z0=-400 cb_z1=500 beta_u0=100"
