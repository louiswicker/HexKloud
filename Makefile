#
FC      = gfortran
NCAR_FC = ncargf90
#
FFLAGS  = -O2 -fopenmp

OUTPUTINC = ${NETCDF_INC}
LINKOPTS  = ${NETCDF_LIB}

NCDF_CODE  = ncdf.inc.f90

NCARG_CODE = plotting.inc.f90

.SUFFIXES: .f90 .F90

OBJS = module_mp_nssl_2mom.o weno.o rk3_grid.o ncdf_utils.o\
       rk3_init.o rk3_rhss.o rk3_rhsu1.o \
       rk3_rhsu2.o rk3_rhsu3.o  \
       rk3_rhsw.o rk3_smlstep.o rk3_coefs.o kessler.o 
OBJS_NCDF = rk3_main.o $(OBJS)
OBJS_NCAR = rk3_main_ncarg.o rk3_plot.o $(OBJS)

%.o : %.f90 ; $(FC) $(FFLAGS) -c $*.f90

%.o : %.F90 ; $(FC) $(FFLAGS) -c $*.F90

rk3_ncdf: $(OBJS_NCDF) $(NCDF_CODE)
	$(FC) $(FFLAGS) -I$(OUTPUTINC) -o rk3_ncdf $(OBJS_NCDF) $(LINKOPTS)

rk3_ncar:  $(OBJS_NCAR) $(OUTPUT_CODE)
	$(NCAR_FC) $(FFLAGS) -I$(OUTPUTINC) -o rk3_ncar $(OBJS_NCAR) $(LINKOPTS)

rk3_main.o: rk3_main.f90 $(OUTPUT_CODE) module_mp_nssl_2mom.o rk3_grid.o \
             initialize.inc.f90 boundaries.inc.f90 ncdf.inc.f90
	$(FC) -c $(FFLAGS) rk3_main.F90 

rk3_main_ncarg.o: rk3_main.f90 $(OUTPUT_CODE) initialize.inc.f90 boundaries.inc.f90 plotting.inc.f90
	$(NCAR_FC) -c $(FFLAGS) -DUSENCARG rk3_main.F90 -o rk3_main_ncarg.o

rk3_grid.o: rk3_grid.f90 
	$(FC) -c $(FFLAGS) rk3_grid.f90 

ncdf_utils.o: ncdf_utils.f90
	$(FC) $(FFLAGS) -I$(OUTPUTINC) -c $*.F90

weno.o: weno.f90 
	$(FC) -c $(FFLAGS) weno.f90 

clean:
	rm -f *.o *.mod rk3_ncar rk3_ncdf
