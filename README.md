# HexKloud

Joe Klemp's toy cloud model modified by NSSL scientists


To compile you need the following environment variables set:

NETCDF_INC  --> points to your netCDF4 include directory
NETCDF_LIB  --> points to your netCDF4 library directory

!=========================================
Compiler

->requirements: GNU fortran compiler system (tested with GNU 12)

->cmd:  'make'  --> outputs NCDF version of code

!=========================================
NCARG version (makes ncdf and ncar graphics)

->requirements:  valid ncar graphics installation and environment variables

->cmd:  'make rk3_ncar'  --> output NCDF+NCARG version of code


!=========================================
# NAMELIST INFO (NSSL added quite a few namelist variables to make it more runtime configurable)

!=========================================
! "main" namelist
!=========================================
!
!  HORIZONTAL/VERTICAL ADVECTION CONTROLS
!-----------------------------------------
!  h_mom_adv, v_mom_adv, 
!  h_sca_adv, v_sca_adv 
!
!   Valid horizontal 2,3,4,5,6,33,55 (third-order,5th-order WENO)
!   Vertical advection is 2nd or 3rd order
!
!  MP_PHYSICS: microphysics: 
!
!    0=dry (vapor only)
!    1=kessler
!   18= NSSL 2-moment
!-----------------------------------------
!  INITAL CONDITIONS
!
!    delt : bubble temperature
!-----------------------------------------
!
!  OUTPUT controls
!-----------------------------------------
!
!  doplot: True/False --> flag for ncargraphics output
!
!  iwty:   type of ncar output 
!     1=gmeta
!    20=postscript for ncargraphics output
!
!  writenc = flag for netcdf output
!  ncuopt  = option for U output in netcdf (2=all hex centers; 3=averaged to x points)
!  ncupert = Whether netcdf output has full U (0) or perturbation U (1)

&main
  mp_physics = 18,
  nssl_2moment_on = 1,
!  nssl_3moment = 1,
  h_mom_adv = 5
  v_mom_adv = 3
  h_sca_adv = 55
  v_sca_adv = 3
  delt = 3.0
  iwty = 20
  debug = .false.
  runname = 'n2mweno',
  writenc = .true.,
  doplot = .true.,
  ncuopt = 3,
  ncupert = 0,
/
!
!=========================================
! "gridtime" namelist
!=========================================
!
!   dt   : time step
!   tstp : total time (seconds) (multiple of tip and dt)
!   tip  : output interval (seconds) (must be divisible by dt)
!
!-----------------------------------------
!
!   nx,ny,nz : grid dimensions

!     xl : east-west domain size
!                                d = 2.*xl/(sqrt(3.)*float(nx-1))
!     xl = ((sqrt(3.)/2.)*(nx-1)*d
!     yl = d*(ny-1)
!-----------------------------------------
!
&gridtime
  dt   = 10.0
  tstp = 7200,
  tip  =  600.,
  nz   = 41,
  nx   = 91,
  ny   = 79,
  xl   = 84000.,
/
!
!=========================================
! "nssl_mp_params" namelist
!=========================================
!
&nssl_mp_params
!  infall = 0,
/
