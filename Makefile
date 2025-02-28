VPATH = ../ncarg_arm/bin

.f.o:
	gfortran -c  $*.f
rk3_moist: rk3_main.o rk3_init.o rk3_plot.o rk3_rhss.o rk3_rhsu1.o \
           rk3_rhsu2.o rk3_rhsu3.o\
           rk3_rhsw.o rk3_smlstep.o rk3_coefs.o kessler.o \
           plotting.inc initialize.inc boundaries.inc
	ncarg_arm/bin/ncargf90 -o rk3_moist  rk3_main.o \
           rk3_init.o rk3_plot.o rk3_rhss.o rk3_rhsu1.o \
           rk3_rhsu2.o rk3_rhsu3.o \
           rk3_rhsw.o rk3_smlstep.o rk3_coefs.o kessler.o

rk3_main.o: rk3_main.f plotting.inc initialize.inc boundaries.inc
	ncarg_arm/bin/ncargf90 -c rk3_main.f 

rk3_plot.o: rk3_plot.f
	ncarg_arm/bin/ncargf90 -c rk3_plot.f 

clean:
	rm -f *.o rk3_moist
