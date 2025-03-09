
!===========================================================
!
!   WRITE_NC4_FILE (simple netCDF4 write utility)
!
!===========================================================


    SUBROUTINE WRITE_NC4_FILE(filename, nx, ny, nz, nv, x, y, z, var, label, time)

    USE netcdf

    implicit none
 
    character(len=*),    intent(in) :: filename
    integer,             intent(in) :: nx, ny, nz, nv
    real,                intent(in) :: time

    real, dimension(nx), intent(in)            :: x
    real, dimension(ny), intent(in)            :: y
    real, dimension(nz), intent(in)            :: z
    real, dimension(nx,ny,nz,nv), intent(in)   :: var
!   real, dimension(nx,ny,nz), intent(in)      :: z

    character(len=*), dimension(nv), intent(in) :: label
    
! Local declarations

    integer :: n
    integer :: shuffle = 1, deflate= 1, deflate_level = 2

    integer :: ncid, status
    integer :: nxDimID, nyDimID, nzDimID, ntDimID
    integer ::  xVarID,  yVarID,  zVarID, tVarID

    integer :: VarID(nv)
    
    real, dimension(nx,ny,nz,1) :: tmp

    write(6,*) 'SUBROUTINE WRITE_NC4_FILE: ', filename
    
!----------------- Create and open netCDF -----------------

!    status = nf90_create(trim(filename),NF90_64BIT_OFFSET,ncid)
    status = nf90_create(trim(filename),NF90_NETCDF4,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Define dimensions and variables -------------

    status = nf90_def_dim(ncid,"xh",nx,nxDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
    
    status = nf90_def_dim(ncid,"yh",ny,nyDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"zh",nz,nzDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"time",1,ntDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"xh",nf90_float,(/nxDimID/), xVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"yh",nf90_float,(/nyDimID/), yVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
    
!   status = nf90_def_var(ncid,"zh",nf90_float,(/nxDimID,nyDimID,nzDimID/), zVarID)
    status = nf90_def_var(ncid,"zh",nf90_float,(/nzDimID/), zVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"time",nf90_float,(/ntDimID/), tVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
      
    do n = 1,nv
       ! write(0,*) 'write variable n = ',n
       status = nf90_def_var(ncid, label(n), nf90_float, (/nxDimID,nyDimID,nzDimID,ntDimID/),VarID(n))
       if(status /= nf90_NoErr) write(*,*) label(n), nf90_strerror(status)
       status = NF90_DEF_VAR_DEFLATE(ncid, VarID(n), shuffle, deflate, deflate_level)
       if(status /= nf90_NoErr) write(*,*) label(n), nf90_strerror(status)
      
    enddo

    status = nf90_enddef(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Write coordinates and variables -------------

    status = nf90_put_var(ncid, xVarID, x)
    if(status /= nf90_NoErr) write(*,*) 'X-COORD: ', nf90_strerror(status)

    status = nf90_put_var(ncid, yVarID, y)
    if(status /= nf90_NoErr) write(*,*) 'Y-COORD: ', nf90_strerror(status)

    status = nf90_put_var(ncid, zVarID, z)
    if(status /= nf90_NoErr) write(*,*) 'Z-COORD: ', nf90_strerror(status)

    status = nf90_put_var(ncid, tVarID, time)
    if(status /= nf90_NoErr) write(*,*) 'Z-COORD: ', nf90_strerror(status)

    do n = 1,nv
       tmp(:,:,:,1) = var(:,:,:,n)
        
       status = nf90_put_var(ncid, VarID(n), tmp)
       if(status /= nf90_NoErr) write(*,*) label(n), nf90_strerror(status)
      
    enddo

    status = nf90_close(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    RETURN
    END SUBROUTINE WRITE_NC4_FILE

!===========================================================
!
!   READ_NC4_FILE (simple netCDF4 read utility)
!
!===========================================================

    SUBROUTINE READ_NC4_FILE( filename, nt, nx, ny, nz, xc, yc, zc, time, den ) 

    use netcdf

    implicit none

    integer, intent(in) :: nt, nx, ny, nz
    character(len=*), intent(in) :: filename

    real, dimension(nx), intent(out) :: xc
    real, dimension(ny), intent(out) :: yc
    real, dimension(nt), intent(out) :: time

    real, dimension(nt,nx,ny,nz), intent(out) :: den, zc

    integer :: k
    integer :: varid, ncid, status

!------------------Open netCDF-----------------------------

    status = nf90_open(trim(filename),nf90_nowrite,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!----------Get 1D variables needed from netcdf----------------

    status = nf90_inq_varid(ncid,"time",varid)
    status = nf90_get_var(ncid,varid,time,start=(/1/),count=(/nt/))

    status = nf90_inq_varid(ncid,"zh",varid)
    status = nf90_get_var(ncid,varid,zc,start=(/1/),count=(/nx,ny,nz,nt/))

    status = nf90_inq_varid(ncid,"yh",varid)
    status = nf90_get_var(ncid,varid,yc,start=(/1/),count=(/ny/))

    status = nf90_inq_varid(ncid,"xh",varid)
    status = nf90_get_var(ncid,varid,xc,start=(/1/),count=(/nx/))

!----------Get 3D variables needed from netcdf----------------

    status = nf90_inq_varid(ncid,"den",varid)

    IF( status /= nf90_NoErr) THEN
        write(*,*) ' ----> Retrieve_Beta/READNC2: No 3D density in file, stopping'
        stop 999
    ELSE
      status = nf90_get_var(ncid,varid,den,start=(/1,1,1,1/),count=(/nx,ny,nz,nt/))
    ENDIF

!------------------Close netCDF----------------------------

    status = nf90_close(ncid)

    RETURN
    END SUBROUTINE READ_NC4_FILE

!===========================================================
!
!   READ_NC4_DIMS (simple netCDF4 read utility to return dims)
!
!===========================================================

    SUBROUTINE READ_NC4_DIMS( filename, nt, nx, ny, nz )

    use netcdf

    implicit none
 
    integer, intent(out) :: nt, nx, ny, nz
    character(len=*), intent(in) :: filename
 
    integer :: status, ncid, dimid
 
!------------------Open netCDF-----------------------------
 
    status = nf90_open(trim(filename),nf90_nowrite, ncid)
 
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
 
! Read NT (unlimited dimension)
 
    status = nf90_inquire(ncid, unlimiteddimid = dimid)
    if (status /= nf90_noerr) write(*,*) 'No unlimited dimension, ', nf90_strerror(status)

    status = nf90_inq_dimid( ncid, "nt", dimid )
    if (status /= nf90_noerr) write(*,*) "Cannot find NT, ", nf90_strerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len = nt)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
! Read NZ
 
    status = nf90_inq_dimid( ncid, "nz", dimid )
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len = nz)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
! Read NY
 
    status = nf90_inq_dimid( ncid, "ny", dimid )
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len = ny)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
! READ NX

    status = nf90_inq_dimid( ncid, "nx", dimid )
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len = nx)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)

    status = nf90_close(ncid)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
    WRITE(6,*) nx, ny, nz, nt

    RETURN
    END SUBROUTINE READ_NC4_DIMS

!===========================================================
!
!   READ_NC4_ATTR (simple netCDF4 utility to read global attrs)
!
!===========================================================

    SUBROUTINE READ_NC4_ATT( filename, wbc, ebc, sbc, nbc, model )

    use netcdf

    implicit none
 
    character(len=*), intent(in) :: filename

    integer, intent(out) :: wbc, ebc, sbc, nbc
    character(len=*), intent(out) :: model
 
    integer :: status, ncid, dimid
 
!------------------Open netCDF-----------------------------
 
    status = nf90_open(trim(filename),nf90_nowrite,ncid)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
 
! Read WBC 
 
    status = nf90_get_att(ncid, NF90_GLOBAL, "wbc", wbc)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
! Read EBC
 
    status = nf90_get_att(ncid, NF90_GLOBAL, "ebc", ebc)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
! Read SBC
 
    status = nf90_get_att(ncid, NF90_GLOBAL, "sbc", sbc)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
! Read NBC
 
    status = nf90_get_att(ncid, NF90_GLOBAL, "nbc", nbc)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)

! Read MODEL type
 
    status = nf90_get_att(ncid, NF90_GLOBAL, "model", model)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)

! Close file
 
    status = nf90_close(ncid)
    if (status /= nf90_noerr) write(*,*) nf90_strerror(status)
 
!   WRITE(6,*) wbc, ebc, sbc, nbc

    RETURN
    END SUBROUTINE READ_NC4_ATT
