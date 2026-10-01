## How to run

### Spinup
./runme -rs -q shared -e yelmox -w 2-00:00:00 -m 10G -n par/yelmo_Antarctica_spinup_32KM_good.nml -o "/scratch/b/b383705/outputs/spinup"

### Transient run
./runme -rs -q shared -e yelmox -w 2-00:00:00 -m 10G -n par/yelmo_Antarctica_lgp_32KM_good.nml -o "/scratch/b/b383705/outputs/lgp" -p yelmo.restart="/scratch/b/b383705/outputs/spinup/restart-15.000-kyr/yelmo_restart.nc" isos.restart="/scratch/b/b383705/outputs/spinup/restart-15.000-kyr/isos_restart.nc" marine_shelf.restart="/scratch/b/b383705/outputs/spinup/restart-15.000-kyr/marine_shelf.nc"

## Forcing index

The file alpha_combined_125kyr_interp.dat should be stored in input/ to be accessed by the model. See XXXXX for how the index was created. 
