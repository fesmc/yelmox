#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/geothermal.o: $(libdir)/geothermal.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/hyster.o: $(libdir)/hyster.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/latinhypercube.o: $(libdir)/latinhypercube.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/marine_shelf.o: $(libdir)/marine_shelf.f90 $(objdir)/pico.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/ismip6.o: $(libdir)/ismip6.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/esm_forcing.o: $(libdir)/esm_forcing.f90 $(objdir)/marine_shelf.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# $(objdir)/isostasy.o: $(libdir)/isostasy.f90 $(objdir)/nml.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

#$(objdir)/ncio.o: $(libdir)/ncio.f90
#	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

# $(objdir)/nml.o: $(libdir)/nml.f90
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/gaussian_filter.o: $(libdir)/gaussian_filter.f90
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sediments.o: $(libdir)/sediments.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snapclim.o: $(libdir)/snapclim.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snapesm.o: $(libdir)/snapesm.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/climate_out.o: $(libdir)/climate_out.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# REMBO for the rembo climate backend: the adapter over rembo1 with rembo=1,
# else a stub with the same interface that stops if climate = "rembo".
$(objdir)/climate_rembo.o: $(libdir)/climate_rembo.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_REMBO) -c -o $@ $<

$(objdir)/climate_rembo_stub.o: $(libdir)/climate_rembo_stub.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# The climate backend of a domain ([comps] climate = snapclim | snapesm | esm
# | rembo), chosen at runtime; the domain reads dom%clim, filled by yelmox_climate.
$(objdir)/yelmox_climate.o: $(libdir)/yelmox_climate.f90 $(objdir)/climate_out.o \
						$(objdir)/snapclim.o $(objdir)/snapesm.o $(objdir)/esm_forcing.o \
						$(objdir)/marine_shelf.o $(objdir)/kryos_forcing.o $(climate_rembo_obj)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) -c -o $@ $<

# Kryos: the domain (kryos_domain + config + init + remap) and the modules built
# on it -- region-specific physics, per-step coupling, cold start + restarts,
# output, and the driver-owned transient forcing.
$(objdir)/kryos.o: $(libdir)/kryos.f90 $(objdir)/marine_shelf.o \
						$(objdir)/climate_out.o $(objdir)/yelmox_climate.o \
						$(objdir)/smbpal.o $(objdir)/smb_simple.o \
						$(objdir)/surface_chion.o \
						$(objdir)/sediments.o $(objdir)/geothermal.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) -c -o $@ $<

$(objdir)/kryos_regions.o: $(libdir)/kryos_regions.f90 $(objdir)/kryos.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) -c -o $@ $<

$(objdir)/kryos_coupling.o: $(libdir)/kryos_coupling.f90 $(objdir)/kryos.o \
						$(objdir)/kryos_regions.o $(objdir)/kryos_forcing.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) -c -o $@ $<

$(objdir)/kryos_forcing.o: $(libdir)/kryos_forcing.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) -c -o $@ $<

$(objdir)/kryos_startup.o: $(libdir)/kryos_startup.f90 $(objdir)/kryos.o \
						$(objdir)/kryos_regions.o $(objdir)/kryos_coupling.o \
						$(objdir)/kryos_forcing.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) -c -o $@ $<

# CMIP/ISMIP-formatted output ([output] write_cmip)
$(objdir)/cmip_output.o: $(libdir)/cmip_output.f90 $(objdir)/marine_shelf.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) -c -o $@ $<

$(objdir)/kryos_output.o: $(libdir)/kryos_output.f90 $(objdir)/kryos.o $(objdir)/cmip_output.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) -c -o $@ $<

# Bipolar ocean coupling: bridge over kryos_domain + the obm box model. Lives
# alongside the bipolar driver in yelmox_bipolar/ -- it is only pertinent to
# that program -- and is linked via obm_libs (bipolar targets only).
$(objdir)/obm_coupling.o: yelmox_bipolar/obm_coupling.f90 $(objdir)/kryos.o \
						$(objdir)/obm_defs.o $(objdir)/ice2ocean.o $(objdir)/ocean2ice.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) -c -o $@ $<

# $(objdir)/stommel.o: $(libdir)/stommel.f90 $(objdir)/yelmo_defs.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/timeout.o: $(libdir)/timeout.f90 $(objdir)/nml.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/timer.o: $(libdir)/timer.f90
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/timestepping.o: $(libdir)/timestepping.f90 $(objdir)/nml.o $(objdir)/ncio.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/varslice.o: $(libdir)/varslice.f90 $(objdir)/nml.o $(objdir)/ncio.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/xarray.o: $(libdir)/xarray.f90 $(objdir)/nml.o $(objdir)/ncio.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# insol library 
$(objdir)/interp1D.o: $(libdir)/insol/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(libdir)/insol/insolation.f90 $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<


# smbpal library
$(objdir)/smbpal_precision.o: $(libdir)/smbpal/smbpal_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_itm.o: $(libdir)/smbpal/smb_itm.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/smb_pdd.o: $(libdir)/smbpal/smb_pdd.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(libdir)/smbpal/interp_time.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smbpal.o: $(libdir)/smbpal/smbpal.f90 $(objdir)/smbpal_precision.o $(objdir)/insolation.o  \
					$(objdir)/interp1D.o  $(objdir)/interp_time.o \
					$(objdir)/smb_pdd.o $(objdir)/smb_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# chion surface model wrapper (surface_method = "chion"): annual cycle of daily
# chion steps on grid_surface, host-supplied insolation.
$(objdir)/surface_chion.o: $(libdir)/surface_chion.f90 $(objdir)/insolation.o \
						$(CHIONROOT)/libchion/include/libchion.a
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_CHION) -c -o $@ $<

# smb_simple library
$(objdir)/smb_simple.o: $(libdir)/smb_simple.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# pico library
$(objdir)/pico_geometry.o: $(libdir)/pico/pico_geometry.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/pico_physics.o: $(libdir)/pico/pico_physics.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/pico.o: $(libdir)/pico/pico.f90 $(objdir)/pico_geometry.o $(objdir)/pico_physics.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# oceanic models for bipolar mode
$(objdir)/obm_defs.o: $(libdir)/obm/obm_defs.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/ice2ocean.o: $(libdir)/obm/ice2ocean.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/ocean2ice.o: $(libdir)/obm/ocean2ice.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/atm2ocean.o: $(libdir)/obm/atm2ocean.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/stommel.o: $(libdir)/obm/stommel.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/nautilus.o: $(libdir)/obm/nautilus.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/obm.o: $(libdir)/obm/obm.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# General yelmox helper modules for different applications
$(objdir)/yelmox_hysteresis_help.o: yelmox_hysteresis_help.f90 $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_YELMO) -c -o $@ $^

#############################################################
##							
## List of library files
##
#############################################################

yelmox_libs = 			$(objdir)/geothermal.o \
					    $(objdir)/hyster.o \
					    $(objdir)/interp1D.o \
					    $(objdir)/insolation.o \
					    $(objdir)/ismip6.o \
					    $(objdir)/esm_forcing.o \
					    $(objdir)/marine_shelf.o \
                        $(objdir)/pico_geometry.o \
                        $(objdir)/pico_physics.o \
                        $(objdir)/pico.o \
					    $(objdir)/sediments.o \
			 		    $(objdir)/smbpal_precision.o \
						$(objdir)/interp_time.o \
					    $(objdir)/smb_itm.o \
					    $(objdir)/smb_pdd.o \
					    $(objdir)/smbpal.o \
					    $(objdir)/smb_simple.o \
					    $(objdir)/surface_chion.o \
					    $(objdir)/climate_out.o \
					    $(climate_rembo_obj) \
					    $(objdir)/yelmox_climate.o \
					    $(objdir)/snapclim.o \
					    $(objdir)/snapesm.o \
						$(objdir)/kryos.o \
						$(objdir)/kryos_regions.o \
						$(objdir)/kryos_coupling.o \
						$(objdir)/kryos_forcing.o \
						$(objdir)/kryos_startup.o \
						$(objdir)/cmip_output.o \
						$(objdir)/kryos_output.o

# Ocean box model stack + its kryos_domain coupling bridge: bipolar-only, linked
# on top of yelmox_libs by the yelmox_bipolar targets.
obm_libs = 				$(objdir)/obm_defs.o\
						$(objdir)/ice2ocean.o\
						$(objdir)/ocean2ice.o\
						$(objdir)/atm2ocean.o\
						$(objdir)/stommel.o\
						$(objdir)/nautilus.o\
						$(objdir)/obm.o\
						$(objdir)/obm_coupling.o

yelmox_help = 			$(objdir)/yelmox_hysteresis_help.o

# yelmox objects embed fesm-utils, yelmo, isostasy (and rembo) types: rebuild
# them when one of these archives changes, or objects keep a stale type layout
# (fesmc/FastHydrology#10). libyelmo.a covers FastHydrology, elsa and tracer.
$(yelmox_libs) $(obm_libs): $(FESMUTILSLIBDIR)/libfesmutils.a \
                            $(YELMOROOT)/libyelmo/include/libyelmo.a \
                            $(ISOSTASYROOT)/libisostasy/include/libisostasy.a

$(objdir)/climate_rembo.o: $(REMBOROOT)/librembo/include/librembo.a
