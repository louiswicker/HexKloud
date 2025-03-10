MODULE RK3_GRID
     
     implicit none
     
!       common /grid/ xh, xu1, xu2, xu3, yh, yu1, yu2, yu3

      real, allocatable, dimension(:,:) :: xh,xu1,xu2,xu3, yh,yu1,yu2,yu3
      integer ::  nx= 91, ny= 79, nz= 41
     
END MODULE RK3_GRID
