#
FC      = gfortran
NCAR_FC = ncargf90
#
FFLAGS  = -O2

OUTPUTINC = ${NETCDF_INC}
LINKOPTS  = ${NETCDF_LIB}

NCDF_CODE  = ncdf.inc.f90

NCARG_CODE = plotting.inc.f90

.SUFFIXES: .f90 .F90

OBJS = module_mp_nssl_2mom.o weno.o ncdf_utils.o\
       rk3_main.o rk3_init.o rk3_rhss.o rk3_rhsu1.o \
       rk3_rhsu2.o rk3_rhsu3.o \
       rk3_rhsw.o rk3_smlstep.o rk3_coefs.o kessler.o 

.f90.o:
	$(FC) $(FFLAGS) -I$(OUTPUTINC) -c $*.f90

.F90.o:
	$(FC) $(FFLAGS) -I$(OUTPUTINC) -c $*.F90

rk3_ncdf:  $(OBJS) $(NCDF_CODE) initialize.inc.f90 boundaries.inc.f90
	$(FC) $(FFLAGS) -I$(OUTPUTINC) -o rk3_ncdf $(OBJS) $(LINKOPTS)

rk3_ncar:  $(OBJS) $(OUTPUT_CODE) initialize.inc.f90 boundaries.inc.f90
	$(NCAR_FC) $(FFLAGS) -I$(OUTPUTINC) -o rk3_ncar $(OBJS) $(LINKOPTS)

rk3_main.o: rk3_main.f90 $(OUTPUT_CODE) initialize.inc.f90 boundaries.inc.f90
	$(FC) -c $(FFLAGS) rk3_main.f90 

weno.o: weno.f90 
	$(FC) -c $(FFLAGS) weno.f90 

clean:
	rm -f *.o *.mod rk3_ncar rk3_ncdf
