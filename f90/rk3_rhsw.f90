!
!-----------------------------------------------------------------------
!                                                                     72
      subroutine rhs_w( w,w1,fw,ww,p,pb,rt,rtb,rho,ru1,ru2,ru3,  &
     &                  rcv,rb,rqx,nmoist,rqvb,dtsa,g,ds,dts,   &
     &                  rdz,f,xnus,xnusz,nz1,nx,ny,iper,         &
     &                  jper,flux1,flux2,flux3,fluxz,order)

      implicit none

      integer nz1,nx,ny,iper,jper,nmoist
      real    w    (nz1+1,nx,ny),fw   (nz1,nx,  ny),rho  (nz1,nx,ny)  &
     &       ,w1   (nz1+1,nx,ny),p    (nz1,nx,  ny),pb   (nz1,nx,ny)  &
     &       ,ww   (nz1+1,nx,ny),rt   (nz1,nx,  ny),rtb  (nz1,nx,ny)  &
     &       ,ru1  (nz1,0:nx,ny),ru2  (nz1,nx,0:ny),rb   (nz1,nx,ny)  &
     &       ,ru3  (nz1,0:nx,ny),rqx  (nz1,  nx,ny,nmoist)            &
     &       ,ds   (nz1)        ,rqvb (nz1)  &
     &       ,flux1(nz1,0:nx,ny),flux2(nz1,nx,0:ny),flux3(nz1,0:nx,ny)  &
     &       ,fluxz(0:nz1,nx,ny),rcv,dts,dtsa,g,rdz,f,xnus,xnusz,cmoist

      integer nx1,i,ip1,im1,ip2,im2,ip3  &
     &       ,ny1,j,jp1,jm1,jp2,jm2,jp3,jpj,jpm,jpp,jmm  &
     &       ,nz2,k,nz
      character*6 order

      nx1 = nx-1
      ny1 = ny-1
      nz  = nz1+1
!
!        right hand side of w equation
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
!  vertical flux calculation
!
            fluxz(1  ,i,j) = 0.25*(ww(1  ,i,j)+ww(2 ,i,j))  &
     &                                 *(w (1  ,i,j)+w (2 ,i,j))
            fluxz(nz1,i,j) = 0.25*(ww(nz1,i,j)+ww(nz,i,j))  &
     &                                 *(w (nz1,i,j)+w (nz,i,j))

            if (order.eq.'second'.or.order.eq.'fourth') then

               do k=2,nz1-1  
                  fluxz(k,i,j) = .25*(ww(k,i,j)+ww(k+1,i,j))  &
     &                                    *(w (k,i,j)+w (k+1,i,j))
               end do
            else
               do k=2,nz1-1  
                        fluxz(k,i,j) = (1./12.)*  &
     &                              (.5*(ww(k,i,j)+ww(k+1,i,j)) *(  &
     &                7.*(w(k+1,i,j)+w(k,i,j))-(w(k+2,i,j)+w(k-1,i,j)))  &
     &                          +abs(.5*(ww(k,i,j)+ww(k+1,i,j)))*(  &
     &               -3.*(w(k+1,i,j)-w(k,i,j))+(w(k+2,i,j)-w(k-1,i,j))))

!cc               fluxz(k,i,j) = .25*(ww(k,i,j)+ww(k+1,i,j))
!cc     &                                    *(w (k,i,j)+w (k+1,i,j))
               end do
            end if
!
!  horizontal flux calculations
!
            if(order.eq.'second')  then
               do k=2,nz1
                  flux1(k,i,j) = .25*(ru1(k,i,j)+ru1(k-1,i,j  ))  &
     &                                      *(w  (k,i,j)+w  (k,ip1,jpm))
                  flux2(k,i,j) = .25*(ru2(k,i,j)+ru2(k-1,i,j  ))  &
     &                                      *(w  (k,i,j)+w  (k,i  ,jp1))
                  flux3(k,i,j) = .25*(ru3(k,i,j)+ru3(k-1,i,j  ))  &
     &                                      *(w  (k,i,j)+w  (k,ip1,jpj))
               end do
            else if(order.eq.'third '.or.order.eq.'fourth')  then
               do k=1,nz1
                  flux1(k,i,j) =  1./12.*  &
     &                                .5*(ru1(k,i  ,j  )+ru1(k-1,i,j))  &
     &                                *(7.*(w(k,ip1,jpm)+w(k,i  ,j  ))  &
     &                                    -(w(k,ip2,jm1)+w(k,im1,jpj)))
                  flux2(k,i,j) =  1./12.*  &
     &                                .5*(ru2(k,i  ,j  )+ru2(k-1,i,j))  &
     &                                *(7.*(w(k,i  ,jp1)+w(k,i  ,j  ))  &
     &                                    -(w(k,i  ,jp2)+w(k,i  ,jm1)))
                  flux3(k,i,j) =  1./12.*  &
     &                               .5*(ru3(k,i  ,j  )+ru3(k-1,i,j))  &
     &                                *(7.*(w(k,ip1,jpj)+w(k,i  ,j  ))  &
     &                                    -(w(k,ip2,jp1)+w(k,im1,jpm)))
               end do
               if(order.eq.'third ')  then
                  do k=1,nz1
                     flux1(k,i,j) =  flux1(k,i,j)+1./12.*  &
     &                                abs(.5*(ru1(k,i,j)+ru1(k-1,i,j)))  &
     &                               *(-3.*(w(k,ip1,jpm)-w(k,i  ,j  ))  &
     &                                    +(w(k,ip2,jm1)-w(k,im1,jpj)))
                     flux2(k,i,j) =  flux2(k,i,j)+1./12.*  &
     &                                abs(.5*(ru2(k,i,j)+ru2(k-1,i,j)))  &
     &                               *(-3.*(w(k,i  ,jp1)-w(k,i  ,j  ))  &
     &                                    +(w(k,i  ,jp2)-w(k,i  ,jm1)))
                     flux3(k,i,j) =  flux3(k,i,j)+1./12.*  &
     &                                abs(.5*(ru3(k,i,j)+ru3(k-1,i,j)))  &
     &                               *(-3.*(w(k,ip1,jpj)-w(k,i  ,j  ))  &
     &                                    +(w(k,ip2,jp1)-w(k,im1,jpm))) 
                  end do
               end if
            else if(order.eq.'fifth '.or.order.eq.'sixth ')  then
               do k=1,nz1
                  flux1(k,i,j) = (1./60.)*.5*(ru1(k,i,j)+ru1(k-1,i,j))  &
     &                              *(37.*(w(k,ip1,jpm)+w(k,i  ,j  ))  &
     &                                -8.*(w(k,ip2,jm1)+w(k,im1,jpj))   &
     &                                   +(w(k,ip3,jmm)+w(k,im2,jp1)))
                  flux2(k,i,j) = (1./60.)*.5*(ru2(k,i,j)+ru2(k-1,i,j))  &
     &                              *(37.*(w(k,i  ,jp1)+w(k,i  ,j  ))  &
     &                                -8.*(w(k,i  ,jp2)+w(k,i  ,jm1))   &
     &                                   +(w(k,i  ,jp3)+w(k,i  ,jm2)))
                  flux3(k,i,j) = (1./60.)*.5*(ru3(k,i,j)+ru3(k-1,i,j))  &
     &                              *(37.*(w(k,ip1,jpj)+w(k,i,j))  &
     &                                -8.*(w(k,ip2,jp1)+w(k,im1,jpm))   &
     &                                   +(w(k,ip3,jpp)+w(k,im2,jm1)))
               end do
               if(order.eq.'fifth ')  then
                  do k=1,nz1
                   flux1(k,i,j) =flux1(k,i,j)-(1./60.)*  &
     &                               abs(.5*(ru1(k,i,j)+ru1(k-1,i,j)))  &
     &                                  *((w(k,ip3,jmm)-w(k,im2,jp1))  &
     &                                -5.*(w(k,ip2,jm1)-w(k,im1,jpj))   &
     &                               +10.*(w(k,ip1,jpm)-w(k,i  ,j  )))
                           flux2(k,i,j) =flux2(k,i,j)-(1./60.)*  &
     &                               abs(.5*(ru2(k,i,j)+ru2(k-1,i,j)))  &
     &                                  *((w(k,i  ,jp3)-w(k,i  ,jm2))  &
     &                                -5.*(w(k,i  ,jp2)-w(k,i  ,jm1))   &
     &                               +10.*(w(k,i  ,jp1)-w(k,i  ,j  )))
                           flux3(k,i,j) =flux3(k,i,j)-(1./60.)*  &
     &                               abs(.5*(ru3(k,i,j)+ru3(k-1,i,j)))  &
     &                                  *((w(k,ip3,jpp)-w(k,im2,jm1))  &
     &                                -5.*(w(k,ip2,jp1)-w(k,im1,jpm))   &
     &                               +10.*(w(k,ip1,jpj)-w(k,i  ,j  )))
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
            do k=2,nz1
               fw(k,i,j) = - dtsa*(flux1(k,i,j)-flux1(k,im1,jpj)  &
     &                                  +flux2(k,i,j)-flux2(k,i  ,jm1)  &
     &                                  +flux3(k,i,j)-flux3(k,im1,jpm))  &
     &                        - dts*rdz*(fluxz(k,i,j)-fluxz(k-1,i,  j))
            end do
         end do
      end do
!
!     buoyancy term
!
!  123 continue
  
      do j=1,ny
         do i=1,nx
            do k=2,nz1
               cmoist = 0.5*(rho(k,i,j)+rho(k-1,i,j))*  &
     &            (rb (k  ,i,j)+rb (k-1,i,j)+rqvb(k)+rqvb(k-1))/  &
!     &            (rho(k  ,i,j)+rqv(k  ,i,j)+rqc(k  ,i,j)+rqr(k  ,i,j)+  &
     &            (rho(k  ,i,j)+ Sum(rqx(k  ,i,j,1:nmoist)) +  &
!     &             rho(k-1,i,j)+rqv(k-1,i,j)+rqc(k-1,i,j)+rqr(k-1,i,j))
     &             rho(k-1,i,j)+Sum(rqx(k-1,i,j,1:nmoist)))
               fw(k,i,j ) = fw(k,i,j) + 0.5*dts*g*cmoist*  &
     &                      ( ((p(k  ,i,j)-pb(k  ,i,j))/pb (k  ,i,j)  &
     &                                -rcv*rt(k  ,i,j) /rtb(k  ,i,j))  &
     &                      + ((p(k-1,i,j)-pb(k-1,i,j))/pb (k-1,i,j)  &
     &                                -rcv*rt(k-1,i,j)/ rtb(k-1,i,j))  &
     &           -(Sum(rqx(k  ,i,j,1:nmoist))-rqvb(k  ))  &
     &            /(rb(k  ,i,j)+rqvb(k))                            &
     &           -(Sum(rqx(k-1,i,j,1:nmoist))-rqvb(k-1))  &
     &            /(rb(k-1,i,j)+rqvb(k-1))                       )  
            end do
         end do
      end do
         
!        return
!
!  horizontal and vertical mixing terms
!

      do j=1,ny
         jp1 = min(j+1,ny)
         if(jper*j.eq.ny )  jp1 = 2
         jm1 = max(j-1,1)
         if(jper*j.eq.1  )  jm1 = ny1
         do i=1,nx
            if(mod(i,2).eq.0)  then
               jpj = j
               jpm = jm1
            else
               jpj = jp1
               jpm = j
            end if
            ip1 = i+1
            if(iper*i.eq.nx) ip1 = 2
            im1 = i-1
            if(iper*i.eq.1)  im1 = nx1
            do k=2,nz1
               fw(k,i,j) = fw(k,i,j) + xnus  &
     &                  *.5*(rho(k,i  ,j  )+rho(k-1,i,j))  &
     &                      *(w1(k,ip1,jpj)+w1(k,im1,jpm)  &
     &                       +w1(k,ip1,jpm)+w1(k,im1,jpj)  &
     &                       +w1(k,i  ,jp1)+w1(k,i  ,jm1)-6.*w1(k,i,j))
            end do
         end do
      end do
!
!     vertical mixing terms
!
      do j=1,ny
         do i=1,nx
            do k=2,nz1
               fw(k,i,j) = fw(k,i,j)+.5*(rho(k,i,j)+rho(k-1,i,j))  &
     &                   *(xnusz*(w1(k+1,i,j)-2.*w1(k,i,j)+w1(k-1,i,j))  &
     &                    - ds(k)*w1  (k,i,j))
            end do
         end do
      end do
!
      return
      end
!
