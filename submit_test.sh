#!/bin/bash 
#SBATCH -p long
#SBATCH -J yelmo 
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#SBATCH --mem=1800
#SBATCH -t 7-02:30:15

./runylmox -r -e yelmox -o output/deglaciation_ant_quad -n par/yelmo_Antarctica_deglaciation_quad-nl-slope.nml
