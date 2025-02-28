c
c                                                                     72
      subroutine calc_scoef( dtseps, c2, hh, rdz, t, p, tb, 
     *                       rho, rb, rqv, rqc, rqr, rqvb,
     *                       g, rcv, cofwz, coftz, cofwt,
     *                       cofwr, cofwrr, cofrz,
     *                       a, b, c, alpha, gamma, nx,ny,nz1      )
  
      implicit none

      integer nx,ny,nz1
      real dtseps, c2, hh(nx,ny), rdz, g, rcv
      real t    (nz1,nx,ny), p    (nz1,nx,ny), tb    (nz1  ,nx,ny)
     &    ,cofwz(nz1,nx,ny), cofwt(nz1,nx,ny), cofwrr(nz1  ,nx,ny)     
     &    ,a    (nz1,nx,ny), b    (nz1,nx,ny), c     (nz1  ,nx,ny)
     &    ,alpha(nz1,nx,ny), gamma(nz1,nx,ny), coftz (nz1+1,nx,ny)
     &    ,rho  (nz1,nx,ny), rb   (nz1,nx,ny), rqv   (nz1  ,nx,ny)
     &    ,rqc  (nz1,nx,ny), rqr  (nz1,nx,ny)
      real rqvb(nz1), cofwr(nx,ny), cofrz, cmoist
      integer i,j,k,nz,nz2
c
c        coefficients for tri-diagonal matrix
c
      nz = nz1+1
      nz2 = nz1-1
c
      do j=1,ny
         do i=1,nx
            do k=2,nz1
               cmoist =
     &                (rho(k,i,j)+rho(k-1,i,j))/
     &             (rho(k  ,i,j)+rqv(k  ,i,j)+rqc(k  ,i,j)+rqr(k  ,i,j)
     &             +rho(k-1,i,j)+rqv(k-1,i,j)+rqc(k-1,i,j)+rqr(k-1,i,j))

               cofwz (k,i,j) = cmoist*
     &                   dtseps*c2*hh(i,j)**2*rdz*(p(k,i,j)+p(k-1,i,j))
               coftz (k,i,j) = dtseps*rdz*(t(k,i,j)+t(k-1,i,j))
               cofwrr(k,i,j) = cmoist*cofwr(i,j)
            end do
            coftz   (1,i,j) = dtseps*rdz*(3.*t(1  ,i,j)-t(2  ,i,j))
            coftz  (nz,i,j) = dtseps*rdz*(3.*t(nz1,i,j)-t(nz2,i,j))
            do k=1,nz1

               cmoist =
     &                 (rho(k,i,j))/
     &                 (rho(k,i,j)+rqv(k,i,j)+rqc(k,i,j)+rqr(k,i,j))
               cofwt(k,i,j) =   cmoist*(rb(k,i,j)+rqvb(k))*
     &                           dtseps*rcv*hh(i,j)*g/tb(k,i,j)

            end do
            do k=2,nz1
                a(k,i,j) = - cofwz (k  ,i,j)* coftz(k-1,i,j)
     &                     + cofwrr(k  ,i,j)* cofrz
     &                     - cofwt (k-1,i,j)* coftz(k-1,i,j)
                b(k,i,j) =1.+cofwz (k,i,j)*(coftz(k,i,j)+coftz(k  ,i,j))
     &                     - coftz (k,i,j)*(cofwt(k,i,j)-cofwt(k-1,i,j))
                c(k,i,j) = - cofwz (k  ,i,j)* coftz(k+1,i,j)
     &                     - cofwrr(k  ,i,j)* cofrz
     &                     + cofwt (k  ,i,j)* coftz(k+1,i,j)
            end do
            do k=2,nz1
               alpha(k,i,j) = 1./(b(k,i,j)-a(k,i,j)*gamma(k-1,i,j))
               gamma(k,i,j) = c(k,i,j)*alpha(k,i,j)
            end do
         end do
      end do

      return
      end
c
