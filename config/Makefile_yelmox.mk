#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/geothermal.o: $(libdir)/geothermal.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/marine_shelf.o: $(libdir)/marine_shelf.f90 $(objdir)/pico.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/esm_forcing.o: $(libdir)/esm_forcing.f90 $(objdir)/marine_shelf.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/sediments.o: $(libdir)/sediments.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# Namelist helpers of the forcing modules (var_defs lookup, defaults file).
$(objdir)/yelmox_nml.o: $(libdir)/yelmox_nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# The terms of the forcing (coefficient x pattern) and the backend capabilities.
$(objdir)/climate_terms.o: $(libdir)/climate_terms.f90 $(objdir)/climate_out.o \
						$(objdir)/kryos_forcing.o $(objdir)/yelmox_nml.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snapshots.o: $(libdir)/snapshots.f90 $(objdir)/kryos_forcing.o $(objdir)/yelmox_nml.o $(objdir)/climate_terms.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/climate_out.o: $(libdir)/climate_out.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# Climate transforms on climate_out after the backend and the terms (taper, bias), fixed order.
$(objdir)/climate_mods.o: $(libdir)/climate_mods.f90 $(objdir)/climate_out.o \
						$(objdir)/climate_terms.o $(objdir)/yelmox_nml.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# REMBO for the rembo climate backend: the adapter over rembo1 with rembo=1,
# else a stub with the same interface that stops if climate = "rembo".
$(objdir)/climate_rembo.o: $(libdir)/climate_rembo.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_REMBO) -c -o $@ $<

$(objdir)/climate_rembo_stub.o: $(libdir)/climate_rembo_stub.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# The climate backend of a domain ([comps] climate = snapshots | esm | rembo),
# chosen at runtime; the domain reads dom%clim, filled by yelmox_climate.
$(objdir)/yelmox_climate.o: $(libdir)/yelmox_climate.f90 $(objdir)/climate_out.o $(objdir)/climate_terms.o $(objdir)/climate_mods.o \
						$(objdir)/snapshots.o $(objdir)/esm_forcing.o \
						$(objdir)/marine_shelf.o $(objdir)/kryos_forcing.o $(climate_rembo_obj)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) -c -o $@ $<

# Prescribed ice thickness in Yelmo's fixed/relaxation cells ([ice_bnd]), from an ice history.
$(objdir)/ice_bnd.o: $(libdir)/ice_bnd.f90 $(objdir)/kryos_forcing.o $(objdir)/yelmox_nml.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) -c -o $@ $<

# Kryos: the domain (kryos_domain + config + init + remap) and the modules built
# on it -- region-specific physics, per-step coupling, cold start + restarts,
# output, and the driver-owned coefficient registry.
$(objdir)/kryos.o: $(libdir)/kryos.f90 $(objdir)/marine_shelf.o \
						$(objdir)/climate_out.o $(objdir)/yelmox_climate.o $(objdir)/kryos_forcing.o \
						$(objdir)/ice_bnd.o $(objdir)/smbpal.o $(objdir)/smb_simple.o \
						$(objdir)/surface_chion.o \
						$(objdir)/sediments.o $(objdir)/geothermal.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

$(objdir)/kryos_regions.o: $(libdir)/kryos_regions.f90 $(objdir)/kryos.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

$(objdir)/kryos_coupling.o: $(libdir)/kryos_coupling.f90 $(objdir)/kryos.o \
						$(objdir)/kryos_regions.o $(objdir)/kryos_forcing.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

$(objdir)/kryos_forcing.o: $(libdir)/kryos_forcing.f90 $(objdir)/yelmox_nml.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/kryos_startup.o: $(libdir)/kryos_startup.f90 $(objdir)/kryos.o \
						$(objdir)/kryos_regions.o $(objdir)/kryos_coupling.o \
						$(objdir)/kryos_forcing.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

# CMIP/ISMIP-formatted output ([output] write_cmip)
$(objdir)/cmip_output.o: $(libdir)/cmip_output.f90 $(objdir)/marine_shelf.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) -c -o $@ $<

$(objdir)/kryos_output.o: $(libdir)/kryos_output.f90 $(objdir)/kryos.o $(objdir)/cmip_output.o \
						$(objdir)/kryos_forcing.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

# Bipolar ocean coupling: bridge over kryos_domain + the obm box model. Lives
# alongside the bipolar driver in yelmox_bipolar/ -- it is only pertinent to
# that program -- and is linked via obm_libs (bipolar targets only).
$(objdir)/obm_coupling.o: yelmox_bipolar/obm_coupling.f90 $(objdir)/kryos.o \
						$(objdir)/kryos_coupling.o $(objdir)/kryos_forcing.o $(objdir)/climate_terms.o $(objdir)/climate_mods.o \
						$(objdir)/yelmox_nml.o \
						$(objdir)/obm_defs.o $(objdir)/ice2ocean.o $(objdir)/ocean2ice.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_YELMO) $(INC_ISOSTASY) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

# The orbit of the insolation ([orbit]): elements from the insol library.
$(objdir)/orbit_forcing.o: $(libdir)/orbit_forcing.f90 $(objdir)/yelmox_nml.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_INSOL) -c -o $@ $<

# smbpal library
$(objdir)/smbpal_precision.o: $(libdir)/smbpal/smbpal_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_itm.o: $(libdir)/smbpal/smb_itm.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/smb_pdd.o: $(libdir)/smbpal/smb_pdd.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(libdir)/smbpal/interp_time.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/smbpal.o: $(libdir)/smbpal/smbpal.f90 $(objdir)/smbpal_precision.o $(objdir)/orbit_forcing.o \
					$(objdir)/interp_time.o $(objdir)/smb_pdd.o $(objdir)/smb_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_INSOL) -c -o $@ $<

# chion surface model wrapper (surface_method = "chion"): annual cycle of daily
# chion steps on grid_surface, host-supplied insolation.
$(objdir)/surface_chion.o: $(libdir)/surface_chion.f90 $(objdir)/orbit_forcing.o \
						$(CHIONROOT)/libchion/include/libchion.a
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_CHION) $(INC_INSOL) -c -o $@ $<

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
$(objdir)/nautilus.o: $(libdir)/obm/nautilus.f90 $(objdir)/obm_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<
$(objdir)/obm.o: $(libdir)/obm/obm.f90 $(objdir)/obm_defs.o $(objdir)/nautilus.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

#############################################################
##							
## List of library files
##
#############################################################

yelmox_libs = 			$(objdir)/geothermal.o \
					    $(objdir)/orbit_forcing.o \
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
					    $(objdir)/yelmox_nml.o \
					    $(objdir)/climate_terms.o $(objdir)/climate_mods.o \
					    $(climate_rembo_obj) \
					    $(objdir)/yelmox_climate.o \
					    $(objdir)/snapshots.o \
						$(objdir)/ice_bnd.o \
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
						$(objdir)/nautilus.o\
						$(objdir)/obm.o\
						$(objdir)/obm_coupling.o

# yelmox objects embed fesm-utils, yelmo, isostasy (and rembo) types: rebuild
# them when one of these archives changes, or objects keep a stale type layout
# (fesmc/FastHydrology#10). libyelmo.a covers FastHydrology, elsa and tracer.
$(yelmox_libs) $(obm_libs): $(FESMUTILSLIBDIR)/libfesmutils.a \
                            $(YELMOROOT)/libyelmo/include/libyelmo.a \
                            $(ISOSTASYROOT)/libisostasy/include/libisostasy.a

# The insol types (insol_class) are embedded in the surface models' state.
$(objdir)/orbit_forcing.o $(objdir)/smbpal.o $(objdir)/surface_chion.o: $(INSOLROOT)/libinsol/include/libinsol.a

$(objdir)/climate_rembo.o: $(REMBOROOT)/librembo/include/librembo.a
