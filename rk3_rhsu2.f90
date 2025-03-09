!
!                                                                     72
      subroutine rhs_u2( u2,u21,ru2,fu2,ww,rho,ru1,ru3,u1z,u2z,u3z,  &
     &                u1m,u2m,u3m,ds,dtsa,dtsd,dtsf,dts,c1f,  &
     &                c2f,rdz,xnus,xnusz,nz1,nx,ny,iper,jper,  &
     &                   flux1,flux2,flux3,fluxz,h_adv,v_adv)
      implicit none

      integer nz1,nx,ny,iper,jper
      real    u2 (nz1,nx,0:ny),u21(nz1,nx,0:ny),ru1(nz1,0:nx,ny)  &
     &       ,fu2(nz1,nx,0:ny),ru2(nz1,nx,0:ny),ru3(nz1,0:nx,ny)  &
     &       ,flux1(nz1,0:nx,ny),flux2(nz1,nx,0:ny),flux3(nz1,0:nx,ny)  &
     &       ,fluxz(0:nz1,nx,ny),ww (nz1+1,nx,ny),rho(nz1,  nx,ny)  &
     &       ,ds(nz1),u1z(nz1),u2z(nz1),u3z(nz1),u1m,u2m,u3m  &
     &       ,dtsa,dtsd,dtsf,dts,rdz,xnus,xnusz,c1f,c2f

      integer nx1,i,ip1,im1,ip2,im2,ip3  &
     &       ,ny1,j,jp1,jm1,jp2,jm2,jp3,jpj,jpm,jpp,jmm  &
     &       ,nz2,k,nz
      integer h_adv, v_adv

      nx1 = nx-1
      ny1 = ny-1
      nz2 = nz1-1
      nz  = nz1 + 1
!
!     right hand side of u2 equation
!
    do j=1,ny
       jp1 = j+1
       if(jper*j.eq.ny  )  jp1 = 2
        jp2 = jp1+1
        if(jper*j.eq.ny1 )  jp2 = 2
        jp3 = jp2+1
        if(jper*j.eq.ny-2)  jp3 = 2
        jm1 = j-1
        if(jper*j.eq.1  )  jm1 = ny1
        jm2 = jm1-1
        if(jper*j.eq.2  )  jm2 = ny1
        do i=2-iper,nx1+iper
           if(mod(i,2).eq.0)  then
            jpp = jp1
            jpj = j
            jpm = jm1
            jmm = jm2
           else
            jpp = jp2
            jpj = jp1
            jpm = j
            jmm = jm1
           end if
           ip1 = i+1
           if(iper*i.eq.nx  ) ip1 = 2
           ip2 = ip1+1
           if(iper*i.eq.nx1 ) ip2 = 2
           ip3 = ip2+1
           if(iper*i.eq.nx-2) ip3 = 2
           im1 = i-1
           if(iper*i.eq.1  )  im1 = nx1
           im2 = im1-1
           if(iper*i.eq.2  )  im2 = nx1
!
!  vertical flux calculation                                          72
!
           fluxz(0  ,i,j) = 0.
           fluxz(nz1,i,j) = 0.
           fluxz(1  ,i,j) = 0.25*(ww(2  ,i,j)+ww(2  ,i,jp1))  &
     &                              *(u2(2  ,i,j)+u2(1  ,i,j  ))
           fluxz(nz2,i,j) = 0.25*(ww(nz1,i,j)+ww(nz1,i,jp1))  &
     &                              *(u2(nz1,i,j)+u2(nz2,i,j  ))

           if ( v_adv == 2 .or. v_adv == 4 ) then

            do k=2,nz1-2  
               fluxz(k,i,j) = .25*(ww(k+1,i,j)+ww(k+1,i,jp1))  &
     &                                 *(u2(k  ,i,j)+u2(k+1,i,j  ))
            end do
           else
            do k=2,nz1-2  
                     fluxz(k,i,j) = (1./12.)*  &
     &                             (.5*(ww(k+1,i,j)+ww(k+1,i,jp1)) *(  &
     &                              7.*(u2(k+1,i,j)+u2(k  ,i,j  ))  &
     &                                -(u2(k+2,i,j)+u2(k-1,i,j  )))  &
     &                         +abs(.5*(ww(k+1,i,j)+ww(k+1,i,jp1)))*(  &
     &                             -3.*(u2(k+1,i,j)-u2(k  ,i,j  ))  &
     &                                +(u2(k+2,i,j)-u2(k-1,i,j  ))))

            end do
           end if
!
!  horizontal flux calculations
!
            if( h_adv == 2 ) then
               do k=1,nz1
                  flux1(k,i,j) = .25*(ru1(k,i,j)+ru1(k,i  ,jp1))  &
     &                                      *(u2 (k,i,j)+u2 (k,ip1,jpm))
                  flux2(k,i,j) = .25*(ru2(k,i,j)+ru2(k,i  ,jp1))  &
     &                                      *(u2 (k,i,j)+u2 (k,i  ,jp1))
                  flux3(k,i,j) = .25*(ru3(k,i,j)+ru3(k,i  ,jp1))  &
     &                                      *(u2 (k,i,j)+u2 (k,ip1,jpj))
               end do
            else if( h_adv == 3 .or. h_adv == 4 )  then
               do k=1,nz1
                  flux1(k,i,j) =  1./12.*  &
     &                               .5*(ru1(k,i  ,j  )+ru1(k,i  ,jp1))  &
     &                               *(7.*(u2(k,ip1,jpm)+u2(k,i  ,j  ))  &
     &                                   -(u2(k,ip2,jm1)+u2(k,im1,jpj)))
                  flux2(k,i,j) =  1./12.*  &
     &                               .5*(ru2(k,i  ,j  )+ru2(k,i  ,jp1))  &
     &                               *(7.*(u2(k,i  ,jp1)+u2(k,i  ,j  ))  &
     &                                   -(u2(k,i  ,jp2)+u2(k,i  ,jm1)))
                  flux3(k,i,j) =  1./12.*  &
     &                               .5*(ru3(k,i  ,j  )+ru3(k,i  ,jp1))  &
     &                               *(7.*(u2(k,ip1,jpj)+u2(k,i  ,j  ))  &
     &                                   -(u2(k,ip2,jp1)+u2(k,im1,jpm)))
               end do
               if( h_adv == 3 ) then
                  do k=1,nz1
                     flux1(k,i,j) =  flux1(k,i,j)+1./12.*  &
     &                           abs(.5*(ru1(k,i,j)+ru1(k,i  ,jp1)))  &
     &                            *(-3.*(u2 (k,ip1,jpm)-u2 (k,i  ,j  ))  &
     &                                 +(u2 (k,ip2,jm1)-u2 (k,im1,jpj)))
                     flux2(k,i,j) =  flux2(k,i,j)+1./12.*  &
     &                           abs(.5*(ru2(k,i,j)+ru2(k,i  ,jp1)))  &
     &                            *(-3.*(u2 (k,i  ,jp1)-u2 (k,i  ,j  ))  &
     &                                 +(u2 (k,i  ,jp2)-u2 (k,i  ,jm1)))
                     flux3(k,i,j) =  flux3(k,i,j)+1./12.*  &
     &                           abs(.5*(ru3(k,i,j)+ru3(k,i  ,jp1)))  &
     &                            *(-3.*(u2 (k,ip1,jpj)-u2 (k,i  ,j  ))  &
     &                                 +(u2 (k,ip2,jp1)-u2 (k,im1,jpm))) 
                  end do
               end if
            else if( h_adv == 5 .or. h_adv == 6 ) then
               do k=1,nz1
                        flux1(k,i,j) = (1./60.)*  &
     &                              .5*(ru1(k,i  ,j  )+ru1(k,i  ,jp1))  &
     &                           *(37.*(u2 (k,ip1,jpm)+u2 (k,i  ,j  ))  &
     &                             -8.*(u2 (k,ip2,jm1)+u2 (k,im1,jpj))   &
     &                                +(u2 (k,ip3,jmm)+u2 (k,im2,jp1)))
                        flux2(k,i,j) = (1./60.)*  &
     &                              .5*(ru2(k,i  ,j  )+ru2(k,i  ,jp1))  &
     &                           *(37.*(u2 (k,i  ,jp1)+u2 (k,i  ,j  ))  &
     &                             -8.*(u2 (k,i  ,jp2)+u2 (k,i  ,jm1))   &
     &                                +(u2 (k,i  ,jp3)+u2 (k,i  ,jm2)))
                        flux3(k,i,j) = (1./60.)*  &
     &                              .5*(ru3(k,i  ,j  )+ru3(k,i  ,jp1))  &
     &                           *(37.*(u2 (k,ip1,jpj)+u2 (k,i  ,j  ))  &
     &                             -8.*(u2 (k,ip2,jp1)+u2 (k,im1,jpm))   &
     &                                +(u2 (k,ip3,jpp)+u2 (k,im2,jm1)))
               end do
               if( h_adv == 5 )  then
                  do k=1,nz1
                   flux1(k,i,j) =flux1(k,i,j)-(1./60.)*  &
     &                          abs(.5*(ru1(k,i,j)+ru1(k,i  ,jp1)))  &
     &                               *((u2 (k,ip3,jmm)-u2 (k,im2,jp1))  &
     &                             -5.*(u2 (k,ip2,jm1)-u2 (k,im1,jpj))   &
     &                            +10.*(u2 (k,ip1,jpm)-u2 (k,i  ,j  )))
                           flux2(k,i,j) = flux2(k,i,j)-(1./60.)*  &
     &                          abs(.5*(ru2(k,i,j)+ru2(k,i  ,jp1)))  &
     &                               *((u2 (k,i  ,jp3)-u2 (k,i  ,jm2))  &
     &                             -5.*(u2 (k,i  ,jp2)-u2 (k,i  ,jm1))   &
     &                            +10.*(u2 (k,i  ,jp1)-u2 (k,i  ,j  )))
                           flux3(k,i,j) = flux3(k,i,j)-(1./60.)*  &
     &                          abs(.5*(ru3(k,i,j)+ru3(k,i  ,jp1)))  &
     &                               *((u2 (k,ip3,jpp)-u2 (k,im2,jm1))  &
     &                             -5.*(u2 (k,ip2,jp1)-u2 (k,im1,jpm))   &
     &                            +10.*(u2 (k,ip1,jpj)-u2 (k,i  ,j  ))) 
                  end do
               end if
            end if
         end do
      end do
!
!   Flux divergence calculation
!
      do j=1,ny
         jp1 = j+1
         if(jper*j.eq.ny )  jp1 = 2
         jm1 = j-1
         if(jper*j.eq.1  )  jm1 = ny1
         do i=2-iper,nx1+iper
            im1 = i-1
            if(iper*i.eq.1  )  im1 = nx1
            if(mod(i,2).eq.0)  then
               jpj = j
               jpm = jm1
            else
               jpj = jp1
               jpm = j
            end if
            do k=1,nz1
               fu2(k,i,j) =- dtsa*(flux1(k,i,j)-flux1(k,im1,jpj)  &
     &                                 + flux2(k,i,j)-flux2(k,i  ,jm1)  &
     &                                 + flux3(k,i,j)-flux3(k,im1,jpm))  &
     &                        - dts*rdz*(fluxz(k,i,j)-fluxz(k-1,i,  j))
            end do
         end do
      end do
!
!  mixing and Coriolis terms
!
         do j=1,ny
            jp1 = min(j+1,ny)
            if(jper*j.eq.ny )  jp1 = 2
            jp2 = min(jp1+1,ny)
            if(jper*j.eq.ny1)  jp2 = 2
            jm1 = max(j-1,1)
            if(jper*j.eq.1  )  jm1 = ny1
            jm2 = jm1-1
            if(jper*j.eq.2  )  jm2 = ny1
            do i=1,nx1+iper
               if(mod(i,2).eq.0)  then
                  jpp = jp1
                  jpj = j
                  jpm = jm1
                  jmm = jm2
               else
                  jpp = jp2
                  jpj = jp1
                  jpm = j
                  jmm = jm1
               end if
               ip1 = i+1
               if(iper*i.eq.nx) ip1 = 2
               im1 = i-1
               if(iper*i.eq.1)  im1 = nx1
               do k=1,nz1
                  fu2(k,i,j) = fu2(k,i,j)   &
     &                       - ds(k)*.5*(rho(k,i,j)+rho(k,i,jp1))  &
     &                                 *(u21(k,i,j)-u2z(k)-u2m)   &
     &                       + dtsf*(c1f*(ru3(k,i  ,j  )+ru3(k,im1,jpj))  &
     &                              +c2f*(ru3(k,i  ,jp1)+ru3(k,im1,jpm))  &
     &                              +c1f*(ru1(k,i  ,jp1)+ru1(k,im1,jpj))  &
     &                              +c2f*(ru1(k,i  ,j  )+ru1(k,im1,jpp))  &
     &                                 -(rho(k,i  ,j  )+rho(k,i  ,jp1))  &
     &                                 *(u3z(k)+u3m+u1z(k)+u1m))
               end do
               do k=1,nz1
                  fu2(k,i,j) = fu2(k,i,j)  &
     &                            + xnus*.5*(rho(k,i,j)+rho(k,i  ,jp1))  &
     &                   *(u21(k,ip1,jpm)-2.*u21(k,i,j)+u21(k,im1,jpj)  &
     &                    +u21(k,ip1,jpj)-2.*u21(k,i,j)+u21(k,im1,jpm)  &
     &                    +u21(k,i  ,jp1)-2.*u21(k,i,j)+u21(k,i  ,jm1))
               end do
               do k=2,nz2
                  fu2(k,i,j) = fu2(k,i,j)  &
     &                         + xnusz*.5*(rho(k,i,j)+rho(k  ,i,jp1))  &
     &                   *(u21(k+1,i,j)-2.*u21(k,i,j)+u21(k-1,i,j  )  &
     &                    -u2z(k+1)    +2.*u2z(k)    -u2z(k-1))
               end do
            end do
         end do
!
      return
      end
