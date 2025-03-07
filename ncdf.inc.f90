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

! for KESSLER or NSSL-2M runs, these variables will always be outputted

     varlabel(1)  = 'U  '
     varlabel(2)  = 'V  '
     varlabel(3)  = 'W  '
     varlabel(4)  = 'THP'

     varlabel(5)  = 'QvP'
     varlabel(6)  = 'Qc '
     varlabel(7)  = 'Qr'

     IF ( mp_physics == 18 ) THEN  ! NSSL - 2M

         ncdf_nvar = 11

         allocate(ncdf_var(nxpl,nypl,nz1,ncdf_nvar) )

         varlabel(8)  = 'Qi '
         varlabel(9)  = 'Qs '
         varlabel(10) = 'Qh'
         varlabel(11) = 'REF'

     ELSE  ! Kessler microphysics

         ncdf_nvar = 7

         allocate(ncdf_var(nxpl,nypl,nz1,ncdf_nvar) )

     ENDIF

! Create filename for netCDF file
 
     ncdf_file = 'hexkloud.XXXX.nc'
     write(ncdf_file(10:13),100) int(time)
100 format(i4.4)

!===============================================================================

     write(6,*) 'Now writing ', ncdf_file
     write(6,*) ''

!===============================================================================

     DO k = 1,nz1      ! outer k-loop for 3D arrays

! U & V reconstruction

        do j=1,nypl

           jj = j+jpi-1

           if( jper*j.eq.1  )  jjm1 = ny1

           do i=1,nxpl

              ii = i+ipi-1

              if(iper*i.eq.1) iim1 = nx1

              IF( mod(ii,2) == 0) THEN

                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,ii-1,jj)    &
                                                   +u3(k,ii,jj)+u3(k,ii-1,jj-1)  &
                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))

!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii,jj)+u1(k,ii-1,jj)  &
!                                                  +u3(k,ii,jj)+u3(k,ii-1,jj)  &
!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))

              ELSE

                ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii-1,jj)+u1(k,ii,jj+1)  &
                                                   +u3(k,ii,  jj)+u3(k,ii-1,jj)  &
                                               -2.*(u1z(k)+u1m+u3z(k)+u3m))

!               ncdf_var(i,j,k,1) = 0.5/sqrt(3.) * (u1(k,ii-1,jj+1)+u1(k,ii,jj+1)  &
!                                                  +u3(k,ii-1,jj)+u3(k,ii,jj+1)  &
!                                              -2.*(u1z(k)+u1m+u3z(k)+u3m))


              ENDIF

           end do
        end do

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
               ncdf_var(i,j,k,5) = 0.5*(qx(k,ii,jj,lv) + qx(k,ii,jp1,lv) - 2.0*qvzv(k))
               ncdf_var(i,j,k,6) = 1000.0*0.5*(qx(k,ii,jj,lc)+qx(k,ii,jp1,lc))
               ncdf_var(i,j,k,7) = 1000.0*0.5*(qx(k,ii,jj,lr)+qx(k,ii,jp1,lr))

             ELSE

                ncdf_var(i,j,k,5) = 1000.0 * (qx(k,ii,jj,lv) - qvzv(k))
                ncdf_var(i,j,k,6) = 1000.0 * qx(k,ii,jj,lc)
                ncdf_var(i,j,k,7) = 1000.0 * qx(k,ii,jj,lr)

             ENDIF
           ENDDO
         ENDDO

         IF( mp_physics == 18 ) THEN

           DO j=1,nypl
              jj = j+jpi-1
              DO i=1,nxpl
                ii = i+ipi-1
                IF(mod(ii,2) == 0) THEN

                   jp1 =min(jj+1,ny)
                   if(jper*jj.eq.ny)  jp1 = 2

                   ncdf_var(i,j,k,8) = 1000.0*0.5*(qx(k,ii,jj,li)+qx(k,ii,jp1,li))
                   ncdf_var(i,j,k,9) = 1000.0*0.5*(qx(k,ii,jj,ls)+qx(k,ii,jp1,ls))
                   ncdf_var(i,j,k,10)= 1000.0*0.5*(qx(k,ii,jj,lh)+qx(k,ii,jp1,lh))
                   ncdf_var(i,j,k,11)= 0.5*(dbz(k,ii,jj)+dbz(k,ii,jp1))

                ELSE

                   ncdf_var(i,j,k,8) = 1000.0 * qx(k,ii,jj,li)
                   ncdf_var(i,j,k,9) = 1000.0 * qx(k,ii,jj,ls)
                   ncdf_var(i,j,k,10)= 1000.0 * qx(k,ii,jj,lh)
                   ncdf_var(i,j,k,11)= dbz(k,ii,jj)

                ENDIF

              ENDDO
           ENDDO

         ENDIF

     ENDDO   ! outer k-loop

     CALL WRITE_NC4_FILE(ncdf_file, nxpl, nypl, nz1, ncdf_nvar, x, y, zu, ncdf_var, varlabel)

     write(6,*) 'Finished writing ', ncdf_file
     write(6,*) ''

     deallocate(ncdf_var)

     END IF
