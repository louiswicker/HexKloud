!                                                                     72
     iwmax = nxc
     jwmax = nyc
     kwmax = 1
     wmax (nit+1) = 0.
     waxis(nit+1) = nit 
     do j=1,ny
        do i=1,nx
           do k=1,nz1
              IF( w(k,i,j) .gt. wmax(nit+1) ) THEN
                     wmax(nit+1) = w(k,i,j)
                     iwmax = i
                     jwmax = j
                     kwmax = k
               END IF
            end do
         end do
      end do

      IF( (kkk.ge.ip) .or. (nit.eq.0) ) THEN      

        kkk = 0

     if ( mp_physics == 1 ) then ! set dbz for Kessler micro
       do j = 1,ny
         do i = 1,nx
           do k = 1,nz1
             tmp = 2.46e4*(1000.*rqx(k,i,j,lr))**1.27
             if ( tmp > 0.0 ) then
                dbz(k,i,j) = max(0.0, 10.0 * log10(tmp))
             else
                dbz(k,i,j) = 0.0
             endif
           enddo
          enddo
         enddo
     endif

!************  Maximum vertical velocity
!
        write(6,*) 'i,j,k,wmax: ', IWMAX,JWMAX,KWMAX,WMAX(nit+1)
        rtot=0.
        do j=1,ny-jper
           do i=1,nx-iper
              do k=1,nz1
                     rtot = rtot+rr(k,i,j)
              end do 
           end do
        end do

        rtot  = 100.*(rtot-rrtot)/ritot
        write(6,*) '       percent mass change =',rtot

!===============================================================================
! netcdf writeout


     if ( writenc ) then
       write(timestr,'(i5.5)') int(time)
       ncdf_file = trim(runname)//'.'//timestr//'.nc'

! for KESSLER or NSSL-2M runs, these variables will always be outputted

      varlabel(1)  = 'U  '
      varlabel(2)  = 'V  '
      varlabel(3)  = 'W  '
      varlabel(4)  = 'THP'

      varlabel(5)  = 'QvP'
      varlabel(6)  = 'Qc '
      varlabel(7)  = 'Qr'

      IF ( mp_physics == 18 ) THEN  ! NSSL - 2M

       n = 7
       if ( lhl > 1 ) then
         varlabel(n+1)  = 'Qi '
         varlabel(n+2)  = 'Qs '
         varlabel(n+3)  = 'Qh '
         varlabel(n+4)  = 'Qhl'
         n = n+4
       endif

       if ( nssl_2moment_on == 1 ) then
         varlabel(n+1)  = 'ccw'
         varlabel(n+2)  = 'crw'
         varlabel(n+3)  = 'cci'
         varlabel(n+4)  = 'csw'
         varlabel(n+5)  = 'chw'
         varlabel(n+6)  = 'chl'
         varlabel(n+7)  = 'ccn'
         varlabel(n+8)  = 'vh '
         varlabel(n+9)  = 'vhl'
         n = n+9
       endif

       if ( nssl_3moment == 1 ) then
         varlabel(n+1)  = 'zrw'
         varlabel(n+2)  = 'zhw'
         varlabel(n+3)  = 'zhl'
         n = n+3
       endif

       varlabel(n+1) = 'REF'
       n = n+1

       ncdf_nvar = n

       allocate(ncdf_var(nxpl,nypl,nz1,ncdf_nvar) )

      ELSE  ! Kessler microphysics

         ncdf_nvar = 8
         varlabel(8)  = 'REF'

         allocate(ncdf_var(nxpl,nypl,nz1,ncdf_nvar) )

      ENDIF

! Create filename for netCDF file
 
!     ncdf_file = 'hexkloud.XXXX.nc'
!     write(ncdf_file(10:13),100) int(time)
! 100 format(i4.4)

!===============================================================================

     write(6,*) 'Now writing ', ncdf_file
     write(6,*) ''

!===============================================================================

     DO k = 1,nz1      ! outer k-loop for 3D arrays

! U & V reconstruction

        IF ( ncuopt == 1 ) THEN ! version 1 for U

        do j=1,nypl

!                  jp1 = min(j+1,ny)
!                  if(jper*j.eq.ny )  jp1 = 2
!                  jm1 = max(j-1,1)
!                  if(jper*j.eq.1  )  jm1 = ny1

           jj = j+jpi-1

           if( jper*j.eq.1  )  jjm1 = ny1
           if(jper*j.eq.ny )  jp1 = 2

           do i=1,nxpl

              ii = i+ipi-1

              if(iper*i.eq.1) iim1 = nx1

              IF( mod(ii,2) == 0) THEN

                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2

                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,iim1,jj)    &
                                                   +u3(k,ii,jj)+u3(k,iim1,jjm1)  &
                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))

!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,ii-1,jj)  &
!                                                  +u3(k,ii,jj)+u3(k,ii-1,jj)  &
!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))

              ELSE

                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,iim1,jj)+u1(k,ii,jjp1)  &
                                                   +u3(k,ii,  jj)+u3(k,iim1,jj)  &
                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))

!!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii-1,jj+1)+u1(k,ii,jj+1)  &
!!                                                  +u3(k,ii-1,jj)+u3(k,ii,jj+1)  &
!!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))


              ENDIF
           end do
        end do

        
        ELSEIF ( ncuopt == 2 ) THEN ! version 2 for U
        ! copied from x-y ncarge U plot
        do j=1,nypl

                  jp1 = min(j+1,ny)
                  if(jper*j.eq.ny )  jp1 = 2
                  jm1 = max(j-1,1)
                  if(jper*j.eq.1  )  jm1 = ny1

!           jj = j+jpi-1

!           if( jper*j.eq.1  )  jjm1 = ny1

           do i=1,nxpl

                     if(mod(i,2).eq.0)  then
                        jpj = j
                        jpm = jm1
                     else
                        jpj = jp1
                        jpm = j
                     end if
				     im1=i-1
				     if(iper*i.eq.1) im1 = nx1
				     ncdf_var(i,j,k,1) = ampl*.5/sqrt(3.)            &
                                         *(u1(k,i,j)+u1(k,im1,jpj)  &
                                          +u3(k,i,j)+u3(k,im1,jpm)  &
                                       -2.*(u1z(k)+u1m+u3z(k)+u3m))

!              ii = i+ipi-1
!
!              if(iper*i.eq.1) iim1 = nx1
!
!              IF( mod(ii,2) == 0) THEN
!
!                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,ii-1,jj)    &
!                                                   +u3(k,ii,jj)+u3(k,ii-1,jj-1)  &
!                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))
!
!!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,ii-1,jj)  &
!!                                                  +u3(k,ii,jj)+u3(k,ii-1,jj)  &
!!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))
!
!              ELSE
!
!                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii-1,jj)+u1(k,ii,jj+1)  &
!                                                   +u3(k,ii,  jj)+u3(k,ii-1,jj)  &
!                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))
!
!!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii-1,jj+1)+u1(k,ii,jj+1)  &
!!                                                  +u3(k,ii-1,jj)+u3(k,ii,jj+1)  &
!!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))
!
!
!              ENDIF

           end do
        end do

        ELSEIF ( ncuopt == 3 ) THEN ! version 3 for U

        do j=1,nypl

!                  jp1 = min(j+1,ny)
!                  if(jper*j.eq.ny )  jp1 = 2
!                  jm1 = max(j-1,1)
!                  if(jper*j.eq.1  )  jm1 = ny1

           jj = j+jpi-1

           if( jper*j.eq.1  )  jjm1 = ny1
           if(jper*j.eq.ny )  jp1 = 2

           do i=1,nxpl

              ii = i+ipi-1

              if(iper*i.eq.1) iim1 = nx1

              IF( mod(ii,2) == 0) THEN

                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2

                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,iim1,jj)    &
                                                   +u3(k,ii,jp1)+u3(k,iim1,jj)  &
                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))

!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,ii-1,jj)  &
!                                                  +u3(k,ii,jj)+u3(k,ii-1,jj)  &
!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))

              ELSE

                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,iim1,jj)+u1(k,ii,jj)  &
                                                   +u3(k,ii,  jj)+u3(k,iim1,jj)  &
                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))

!!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii-1,jj+1)+u1(k,ii,jj+1)  &
!!                                                  +u3(k,ii-1,jj)+u3(k,ii,jj+1)  &
!!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))


              ENDIF
           end do
        end do

        

           ENDIF ! ncuopt

! W & V (u2 is the merdional velocity due to orientiation of hexagons)

        do j=1,nypl
           jj = j+jpi-1
           do i=1,nxpl
              ii = i+ipi-1
              if(mod(ii,2).eq.0.) then
                 jp1 =min(jj+1,ny)
                 if(jper*jj.eq.ny)  jp1 = 2

                 ncdf_var(i,j,k,2) = 0.5*(u2(k,ii,jj)+u2(k,ii,jp1))
                 ncdf_var(i,j,k,3) = .25*(w(k,ii,jj )+w(k+1,ii,jj )  &
                                         +w(k,ii,jp1)+w(k+1,ii,jp1))
              else
                 ncdf_var(i,j,k,2) = u2(k,ii,jj)
                 ncdf_var(i,j,k,3) = .5*(w(k,ii,jj)+w(k+1,ii,jj))
              end if
           end do
        end do

        do j=1,nypl
           jj = j+jpi-1

           do i=1,nxpl
              ii = i+ipi-1
              IF( mod(ii,2) == 0 ) THEN
                 jp1 =min(jj+1,ny)
                 if(jper*jj.eq.ny)  jp1 = 2
                 ncdf_var(i,j,k,4) = t0*.5*(t(k,ii,jj ) / (1.+1.61*qx(k,ii,jj ,lv))  &
                                           +t(k,ii,jp1) / (1.+1.61*qx(k,ii,jp1,lv))) - t0*tz(k)
              ELSE
                 ncdf_var(i,j,k,4) = t0*t(k,ii,jj)/(1.+1.61*qx(k,ii,jj,lv)) - t0*tz(k)
              END IF
            end do
         end do

         do j=1,nypl
           jj = j+jpi-1

           do i=1,nxpl
             ii = i+ipi-1

             IF( mod(ii,2) == 0 )  THEN

               jp1 =min(jj+1,ny)
               if(jper*jj.eq.ny)  jp1 = 2
               ncdf_var(i,j,k,5) = 1000.*0.5*(qx(k,ii,jj,lv) + qx(k,ii,jp1,lv) - 2.0*qvzv(k))
               do n = 2,nmoist
                  ncdf_var(i,j,k,4+n) = 1000.0*0.5*(qx(k,ii,jj,n)+qx(k,ii,jp1,n))
               enddo
               if ( nscalar > 0 ) then
                 do n = 1,nscalar
                   ncdf_var(i,j,k,4+nmoist+n) = 0.5*(sx(k,ii,jj,n)+sx(k,ii,jp1,n))
                 enddo
               endif
               ncdf_var(i,j,k,ncdf_nvar)= 0.5*(dbz(k,ii,jj)+dbz(k,ii,jp1))

             ELSE

                ncdf_var(i,j,k,5) = 1000.0 * (qx(k,ii,jj,lv) - qvzv(k))
                do n = 2,nmoist
                  ncdf_var(i,j,k,4+n) = 1000.0*qx(k,ii,jj,n)
                enddo
                IF ( nscalar > 0 ) THEN
                  do n = 1,nscalar
                     ncdf_var(i,j,k,4+nmoist+n) = sx(k,ii,jj,n)
                  enddo
                ENDIF
                ncdf_var(i,j,k,ncdf_nvar)= dbz(k,ii,jj)

             ENDIF
           ENDDO
         ENDDO

     ENDDO   ! outer k-loop

     CALL WRITE_NC4_FILE(ncdf_file, nxpl, nypl, nz1, ncdf_nvar, x, y, zu, ncdf_var, varlabel, time)

     write(6,*) 'Finished writing ', ncdf_file
     write(6,*) ''

     deallocate(ncdf_var)
     
     endif ! writenc

     END IF
