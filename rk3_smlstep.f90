!
!-----------------------------------------------------------------------
!                                                                     72
!
      subroutine smlstep( ru1,ru11,fu1,ru2,ru21,fu2,ru3,ru31,fu3,   &
     &                    rw,rw1,fw,t,ts,rt1,ft,rs,rr1,fr,  &
     &                    p,ww,div,du1,du2,du3,hh,gu1,gu2,gu3,  &
     &                    a,alpha,gamma,dhh1,dhh2,dhh3,  &
     &                    cofrz,coftz,cofwz,cofwr,cofwrr,cofwt,  &
     &                    rdz,dts,dtsa,dtsd,c2,smdivx,smdivz,resm,  &
     &                    nx,ny,nz1,ns_input,iper,jper )

      implicit none

      integer nx,ny,nz1,ns,iper,jper,ns_input
      real    ru1  (nz1,0:nx,ny),ru11(nz1,0:nx,ny),fu1 (nz1,0:nx,ny)  &
     &       ,du1  (nz1,0:nx,ny),gu1 (nz1,0:nx,ny),t   (nz1  ,nx,ny)  &
     &       ,ru2  (nz1,nx,0:ny),ru21(nz1,nx,0:ny),fu2 (nz1,nx,0:ny)  &
     &       ,du2  (nz1,nx,0:ny),gu2 (nz1,nx,0:ny),p   (nz1  ,nx,ny)  &
     &       ,ru3  (nz1,0:nx,ny),ru31(nz1,0:nx,ny),fu3 (nz1,0:nx,ny)  &
     &       ,du3  (nz1,0:nx,ny),gu3 (nz1,0:nx,ny),ww  (nz1+1,nx,ny)  &
     &       ,rw   (nz1+1,nx,ny),rw1 (nz1+1,nx,ny),fw  (nz1  ,nx,ny)  &
     &       ,ts   (nz1  ,nx,ny),rt1 (nz1  ,nx,ny),ft  (nz1  ,nx,ny)  &
     &       ,rs   (nz1  ,nx,ny),rr1 (nz1  ,nx,ny),fr  (nz1  ,nx,ny)  &
     &       ,div  (nz1  ,nx,ny),alpha(nz1,nx,ny),gamma(nz1  ,nx,ny)  &
     &       ,a    (nz1  ,nx,ny),cofwz(nz1,nx,ny),coftz(nz1+1,nx,ny)  &
     &       ,cofwt(nz1  ,nx,ny),dhh1     (nx,ny),dhh2       (nx,ny)  &
     &       ,cofwrr(nz1 ,nx,ny),hh       (nx,ny),dhh3       (nx,ny)  &
     &       ,cofwr      (nx,ny),cofrz  &
     &       ,rdz,dts,dtsa,dtsd,c2,smdivx,smdivz,resm

      integer i,j,k,im1,jm1,nx1,ny1,iu1,jv1,ip1,jp1,km1,kp1,jpj,jpm
      integer m,nz2
      logical test0,test1
      parameter (test0 = .false.)
      parameter (test1 = .false.)

!
!  set a few constants
!
      nz2 = nz1-1
      nx1 = nx-1
      ny1 = ny-1
      ns = ns_input

      if(test0) then
!
!  test section
!
         do j=1,ny
            do i=1,nx
               do k=1,nz1
                  ft(k,i,j) = 0.
                  fr(k,i,j) = 0.
                  fw(k,i,j) = 0.
               end do
            end do
         end do
         do j=1,ny
            do i=0,nx
               do k=1,nz1
                  fu1(k,i,j) = 0.
                  fu3(k,i,j) = 0.
               end do
            end do
         end do
         do j=0,ny
            do i=1,nx
               do k=1,nz1
                  fu2(k,i,j) = 0.
               end do
            end do
         end do
         do j=1,ny
            do i=1,nx
               do k=1,nz1+1
                  ww(k,i,j) = 0.
               end do
            end do
         end do

      end if
!
!     remove large time step values from ru1, ru2, ru3, rw1, and omega 
!
      do j=1,ny
         do i=iper,nx
            do k=1,nz1
               ru11(k,i,j) = ru11(k,i,j)-ru1(k,i,j)
               ru31(k,i,j) = ru31(k,i,j)-ru3(k,i,j)
            end do
         end do
      end do
      do j=jper,ny
         do i=1,nx
            do k=1,nz1
               ru21(k,i,j) = ru21(k,i,j)-ru2(k,i,j)
            end do
         end do
      end do
      do j=1,ny
         jv1 = j-1
         if(jper*j.eq.1)  jv1 = ny1
         jm1 = max(j-1,1)
         if(jper*j.eq.1)  jm1=ny1
         jp1 = min(j+1,ny)
         if(jper*j.eq.ny)  jp1=2
         do i=1,nx
            if(mod(i,2).eq.0)  then
               jpj=j
               jpm=jm1
            else
               jpj=jp1
               jpm=j
            end if
            im1=i-1
            if(iper*i.eq.1)  im1=nx1
            rw1(1,i,j) = 0.
            do k=2,nz1
               rw1(k,i,j) = .25*(gu1(k  ,i  ,j  )*ru11(k  ,i  ,j  )  &
     &                          +gu1(k  ,im1,jpj)*ru11(k  ,im1,jpj)  &
     &                          +gu1(k-1,i  ,j  )*ru11(k-1,i  ,j )  &
     &                          +gu1(k-1,im1,jpj)*ru11(k-1,im1,jpj))  &
     &                    + .25*(gu3(k  ,i  ,j  )*ru31(k  ,i  ,j )  &
     &                          +gu3(k  ,im1,jpm)*ru31(k  ,im1,jpm)  &
     &                          +gu3(k-1,i  ,j  )*ru31(k-1,i  ,j )  &
     &                          +gu3(k-1,im1,jpm)*ru31(k-1,im1,jpm))  &
     &                    + .25*(gu2(k  ,i  ,j  )*ru21(k  ,i  ,j  )  &
     &                          +gu2(k  ,i  ,jv1)*ru21(k  ,i  ,jv1)  &
     &                          +gu2(k-1,i  ,j  )*ru21(k-1,i  ,j  )  &
     &                          +gu2(k-1,i  ,jv1)*ru21(k-1,  i,jv1))  &
     &                    + hh(i,j)*(rw1(k,i,j  ) - rw(k  ,i  ,j  ))
               fw(k,i,j) = hh(i,j)*fw(k,i,j)
            end do
         end do
      end do
!
!     take the steps
!
      do m=1,ns
!
!        calculation of divergence
!
         do j=1,ny
            jm1 = max(j-1,1)
            if(jper*j.eq.1 ) jm1 = ny1
            jv1 = j-1
            if(jper*j.eq.1 ) jv1 = ny1
            jp1 = min(j+1,ny)
            if(jper*j.eq.ny) jp1 = 2
            do i=1,nx
               if(mod(i,2).eq.0)  then
                  jpj = j
                  jpm = jm1
               else
                  jpj = jp1
                  jpm = j
               end if
               im1 = max(i-1,1)
               if(iper*i.eq.1 ) im1 = nx1
               iu1 = i-1
               if(iper*i.eq.1 ) iu1 = nx1
               ip1 = min(i+1,nx)
               if(iper*i.eq.nx) ip1 = 2
               do k=1,nz1
                  km1=max(k-1,1)
                  kp1=min(k+1,nz1)
                  div(k,i,j) =  dtsa/dts*.5*  &
     &                         ((ru11(k,i  ,j  )+ru1(k,i  ,j  ))  &
     &                            *(t(k,i  ,j  )+t  (k,ip1,jpm))  &
     &                         -(ru11(k,iu1,jpj)+ru1(k,iu1,jpj))  &
     &                            *(t(k,i  ,j  )+t  (k,im1,jpj))  &
     &                         +(ru21(k,i  ,j  )+ru2(k,i  ,j  ))  &
     &                            *(t(k,i  ,j  )+t  (k,i  ,jp1))  &
     &                         -(ru21(k,i  ,jv1)+ru2(k,i  ,jv1))  &
     &                            *(t(k,i  ,j  )+t  (k,i  ,jm1))  &
     &                         +(ru31(k,i  ,j  )+ru3(k,i  ,j  ))  &
     &                            *(t(k,i  ,j  )+t  (k,ip1,jpj))  &
     &                         -(ru31(k,iu1,jpm)+ru3(k,iu1,jpm))  &
     &                            *(t(k,i  ,j  )+t  (k,im1,jpm)))  &
     &                       + .5*rdz*  &
     &                         ((rw1 (k+1,i,j  )+ww (k+1,i,j  ))  &
     &                            *(t(kp1,i,j  )+t  (k  ,i,j  ))  &
     &                         -(rw1 (k  ,i,j  )+ww (k  ,i,j  ))  &
     &                            *(t(k  ,i,j  )+t  (km1,i,j  )))
!ccd                  div(k,i,j)=  div(k,i,j)/t(k,i,j)

               end do
            end do
         end do
!
!        time step for u1 and u3 equations
!
         if(iper.eq.0)  then
            do j=1,ny
               do k=1,nz1
                  ru11(k,0 ,j) = ru11(k,0 ,j) + fu1(k,0 ,j)
                  du1 (k,0 ,j) = fu1 (k,0 ,j)
                  ru31(k,0 ,j) = ru31(k,0 ,j) + fu3(k,0 ,j)
                  du3 (k,0 ,j) = fu3 (k,0 ,j)
               end do
            end do
            do j=1,ny
               do k=1,nz1
                  ru11(k,nx,j) = ru11(k,nx,j) + fu1(k,nx,j)
                  du1 (k,nx,j) = fu1 (k,nx,j)
                  ru31(k,nx,j) = ru31(k,nx,j) + fu3(k,nx,j)
                  du3 (k,nx,j) = fu3 (k,nx,j)
               end do
            end do
         end if
         do j=1,ny
            jm1 = max(j-1,1)
            if(jper*j.eq.1 ) jm1 = ny1
            jp1 = min(j+1,ny)
            if(jper*j.eq.ny) jp1 = 2
            do i=1,nx1+iper              
               if(mod(i,2).eq.0)  then
                  jpj = j
                  jpm = jm1
               else
                  jpj = jp1
                  jpm = j
               end if
               ip1=i+1
               if(i.eq.nx)  ip1=2
               do k=1,nz1
                  du1 (k,i,j) = ru11(k,i,j)
                  ru11(k,i,j) = ru11(k,i,j) + fu1(k,i,j) - dtsd*c2  &
     &                                 *.5*(p  (k,ip1,jpm)+p  (k,i,j))  &
     &                                    *(rt1(k,ip1,jpm)-rt1(k,i,j)  &
     &                          +dhh1(i,j)*(rt1(k,ip1,jpm)+rt1(k,i,j)))  &
     &                        +dtsd*smdivx*(div(k,ip1,jpm)-div(k,i,j))/  &
     &                                (0.5*(t  (k,ip1,jpm)+t  (k,i,j)))
                  du3 (k,i,j) = ru31(k,i,j)
                  ru31(k,i,j) = ru31(k,i,j) + fu3(k,i,j) - dtsd*c2  &
     &                                 *.5*(p  (k,ip1,jpj)+p  (k,i,j))  &
     &                                    *(rt1(k,ip1,jpj)-rt1(k,i,j)  &
     &                          +dhh3(i,j)*(rt1(k,ip1,jpj)+rt1(k,i,j)))  &
     &                        +dtsd*smdivx*(div(k,ip1,jpj)-div(k,i,j))/  &
     &                                (0.5*(t  (k,ip1,jpj)+t  (k,i,j)))
               end do
               k=1
                  ru11(k,i,j) = ru11(k,i,j) - .125*dts*rdz*c2*  &
     &                       (p(k,ip1,jpm)+p(k,i,j))*gu1(k  ,i  ,  j)*  &
     &                       ((-3.*rt1(k,ip1,jpm)+4.*rt1(k+1,ip1,jpm)  &
     &                                              -rt1(k+2,ip1,jpm))  &
     &                       +(-3.*rt1(k,i  ,j  )+4.*rt1(k+1,i  ,j  )  &
     &                                              -rt1(k+2,i  ,j  )))
                  ru31(k,i,j) = ru31(k,i,j) - .125*dts*rdz*c2*  &
     &                       (p(k,ip1,jpj)+p(k,i,j))*gu3(k  ,i  ,j  )*  &
     &                       ((-3.*rt1(k,ip1,jpj)+4.*rt1(k+1,ip1,jpj)  &
     &                                              -rt1(k+2,ip1,jpj))  &
     &                       +(-3.*rt1(k,i  ,j  )+4.*rt1(k+1,i  ,j  )  &
     &                                              -rt1(k+2,i  ,j  )))
               k=nz1
                  ru11(k,i,j) = ru11(k,i,j) - .125*dts*rdz*c2*  &
     &                       (p(k,ip1,jpm)+p(k,i,j))*gu1(k  ,i  ,j  )*  &
     &                       (( 3.*rt1(k,ip1,jpm)-4.*rt1(k-1,ip1,jpm)  &
     &                                              +rt1(k-2,ip1,jpm))  &
     &                       + (3.*rt1(k,i  ,j  )-4.*rt1(k-1,i  ,j  )  &
     &                                              +rt1(k-2,i  ,j  )))
                  ru31(k,i,j) = ru31(k,i,j) - .125*dts*rdz*c2*  &
     &                       (p(k,ip1,jpj)+p(k,i,j))*gu3(k  ,i  ,j  )*  &
     &                       (( 3.*rt1(k,ip1,jpj)-4.*rt1(k-1,ip1,jpj)  &
     &                                              +rt1(k-2,ip1,jpj))  &
     &                       + (3.*rt1(k,i  ,j  )-4.*rt1(k-1,i  ,j  )  &
     &                                              +rt1(k-2,i  ,j  )))
               do k=2,nz2
                  ru11(k,i,j) = ru11(k,i,j) - .125*dts*rdz*c2*  &
     &                       (p(k,ip1,jpm)+p(k,i,j))*gu1(k  ,i  ,j  )*  &
     &                            ((rt1(k+1,ip1,jpm)-rt1(k-1,ip1,jpm))  &
     &                            +(rt1(k+1,i  ,j  )-rt1(k-1,i  ,j  )))
                  ru31(k,i,j) = ru31(k,i,j) - .125*dts*rdz*c2*  &
     &                       (p(k,ip1,jpj)+p(k,i,j))*gu3(k  ,i  ,j  )*  &
     &                            ((rt1(k+1,ip1,jpj)-rt1(k-1,ip1,jpj))  &
     &                            +(rt1(k+1,i  ,j  )-rt1(k-1,i  ,j  )))
               end do
               do k=1,nz1
                  du1 (k,i,j) = ru11(k,i,j)-du1(k,i,j)
                  du3 (k,i,j) = ru31(k,i,j)-du3(k,i,j)
               end do
            end do
         end do
!
!        time step for u2 equation
!
         if(jper.eq.0)  then
            do i=1,nx
               do k=1,nz1
                  ru21(k,i,0 ) = ru21(k,i,0 ) + fu2(k,i,0 )
                  du2 (k,i,0 ) = fu2 (k,i,0 )
               end do
            end do
            do i=1,nx
               do k=1,nz1
                  ru21(k,i,ny) = ru21(k,i,ny) + fu2(k,i,ny)
                  du2 (k,i,ny) = fu2 (k,i,ny)
               end do
            end do
         end if
         do j=1,ny1+jper              
            jp1=j+1
            if(j.eq.ny)  jp1=2
            do i=1,nx
               do k=1,nz1
                  du2 (k,i,j) = ru21(k,i,j)
                  ru21(k,i,j) = ru21(k,i,j) + fu2(k,i,j) - dtsd*c2*  &
     &                                    .5*(p(k,i,jp1)+p  (k,i,j))  &
     &                                    *(rt1(k,i,jp1)-rt1(k,i,j)  &
     &                          +dhh2(i,j)*(rt1(k,i,jp1)+rt1(k,i,j)))  &
     &                        +dtsd*smdivx*(div(k,i,jp1)-div(k,i,j))/  &
     &                                (0.5*(t  (k,i,jp1)+t  (k,i,j)))
               end do
               k=1
                  ru21(k,i,j) = ru21(k,i,j) - .125*dts*rdz*c2*  &
     &                          (p(k,i,jp1)+p(k,i,j))*gu2(k  ,i,j  )*  &
     &                          ((-3.*rt1(k,i,jp1)+4.*rt1(k+1,i,jp1)  &
     &                                               -rt1(k+2,i,jp1))  &
     &                          +(-3.*rt1(k,i,j  )+4.*rt1(k+1,i,j  )  &
     &                                               -rt1(k+2,i,j  )))
               k=nz1
                  ru21(k,i,j) = ru21(k,i,j) - .125*dts*rdz*c2*  &
     &                          (p(k,i,jp1)+p(k,i,j))*gu2(k  ,i,j  )*  &
     &                           ((3.*rt1(k,i,jp1)-4.*rt1(k-1,i,jp1)  &
     &                                               +rt1(k-2,i,jp1))  &
     &                           +(3.*rt1(k,i,j  )-4.*rt1(k-1,i,j  )  &
     &                                               +rt1(k-2,i,j  )))
               do k=2,nz2
                  ru21(k,i,j) = ru21(k,i,j) - .125*dts*rdz*c2*  &
     &                          (p(k,i,j)+p(k,i,jp1))*gu2(k  ,i,j  )*  &
     &                               ((rt1(k+1,i,jp1)-rt1(k-1,i,jp1))  &
     &                               +(rt1(k+1,i,j  )-rt1(k-1,i,j  )))
               end do
               do k=1,nz1
                  du2 (k,i,j) = ru21(k,i,j)-du2(k,i,j)
               end do
            end do
         end do

         if(.not.test1) then

            do j=1,ny              
               jv1 = j-1
               if(jper*j.eq.1)  jv1 = ny1
               jm1 = max(j-1,1)
               if(jper*j.eq.1)  jm1 = ny1
               jp1 = min(j+1,ny)
               if(jper*j.eq.ny)  jp1 = 2
               do i=1,nx              
                  if(mod(i,2).eq.0)  then
                     jpj = j
                     jpm = jm1
                  else
                     jpj = jp1
                     jpm = j
                  end if
                  iu1 = i-1
                  if(iper*i.eq.1)  iu1 = nx1
                  im1 = max(i-1,1)
                  if(iper*i.eq.1)  im1 = nx1
                  ip1 = min(i+1,nx)
                  if(iper*i.eq.nx) ip1 = 2
!
!                 explicit portion of rr1 time step
!
                  do k=1,nz1
                     rs(k,i,j) = rr1(k,i,j) + fr(k,i,j)  &
     &                          - dtsa*(ru11(k  ,i,j)-ru11(k,iu1,jpj)  &
     &                                 +ru31(k  ,i,j)-ru31(k,iu1,jpm)  &
     &                                 +ru21(k  ,i,j)-ru21(k,i  ,jv1))  &
     &                     - cofrz*resm*(rw1(k+1,i,j)-rw1 (k,i  ,j  ))
                  end do
!
!                 explicit portion of rt1 time step
!
                  do k=1,nz1
                     ts(k,i,j) = rt1(k,i,j) + ft(k,i,j) - dtsa*.5*  &
     &                     (ru11(k,i  ,j  )*(t(k,i  ,j  )+t(k,ip1,jpm))  &
     &                     -ru11(k,iu1,jpj)*(t(k,i  ,j  )+t(k,im1,jpj))  &
     &                     +ru21(k,i  ,j  )*(t(k,i  ,j  )+t(k,i  ,jp1))  &
     &                     -ru21(k,i  ,jv1)*(t(k,i  ,j  )+t(k,i  ,jm1))  &
     &                     +ru31(k,i  ,j  )*(t(k,i  ,j  )+t(k,ip1,jpj))  &
     &                     -ru31(k,im1,jpm)*(t(k,i  ,j  )+t(k,im1,jpm)))  &
     &                           - resm*(coftz(k+1,i,j)*rw1(k+1,i,j  )  &
     &                                 - coftz(k,i  ,j)*rw1(k  ,i,j  ))
                  end do
!
!                 explicit portion of rw1 time step
!
                  do k=2,nz1
                     rw1(k,i,j) = rw1(k,i,j) + fw(k,i,j)  &
     &                       - cofwz (k,i,j)*((ts (k,i,j)-ts (k-1,i,j))  &
     &                                  +resm*(rt1(k,i,j)-rt1(k-1,i,j)))  &
     &                       - cofwrr(k,i,j)*((rs (k,i,j)+rs (k-1,i,j))  &
     &                                  +resm*(rr1(k,i,j)+rr1(k-1,i,j)))  &
     &                + cofwt(k  ,i,j)*(ts (k  ,i,j)+resm*rt1(k  ,i,j))  &
     &                + cofwt(k-1,i,j)*(ts (k-1,i,j)+resm*rt1(k-1,i,j))  &
     &                       + .25*(gu1(k  ,i  ,j  )*du1(k  ,i  ,j  )  &
     &                             +gu1(k  ,iu1,jpj)*du1(k  ,iu1,jpj)  &
     &                             +gu1(k-1,i  ,j  )*du1(k-1,i  ,j  )  &
     &                             +gu1(k-1,iu1,jpj)*du1(k-1,iu1,jpj))  &
     &                       + .25*(gu3(k  ,i  ,j  )*du3(k  ,i  ,j  )  &
     &                             +gu3(k  ,iu1,jpm)*du3(k  ,iu1,jpm)  &
     &                             +gu3(k-1,i  ,j  )*du3(k-1,i  ,j  )  &
     &                             +gu3(k-1,iu1,jpm)*du3(k-1,iu1,jpm))  &
     &                       + .25*(gu2(k  ,i  ,j  )*du2(k  ,i  ,j  )  &
     &                             +gu2(k  ,i  ,jv1)*du2(k  ,i,jv1  )  &
     &                             +gu2(k-1,i  ,j  )*du2(k-1,i  ,j  )  &
     &                             +gu2(k-1,i  ,jv1)*du2(k-1,i  ,jv1))  &
     &              + dts*smdivz*hh(i,j)*rdz*(div(k,i,j)-div(k-1,i,j))/  &
     &                                  (0.5*(t  (k,i,j)+t  (k-1,i,j)))
                  end do
               end do              
            end do
!
!           implicit solution for rw1, rr1, and rt1 
!
            do j=1,ny
               do i=1,nx
                  do k=2,nz1
                     rw1(k,i,j) = (rw1(k,i,j)-a(k,i,j)*rw1(k-1,i,j))  &
     &                            *alpha(k,i,j)
                  end do
                  do k=nz1,1,-1
                     rw1(k,i,j) = rw1(k,i,j)-gamma(k,i,j)*rw1(k+1,i,j)
                     rr1(k,i,j) = rs (k,i,j)  &
     &                              - cofrz*(rw1(k+1,i,j)-rw1(k  ,i,j))
                     rt1(k,i,j) = ts (k,i,j)  &
     &                                   -(coftz(k+1,i,j)*rw1(k+1,i,j)  &
     &                                   - coftz(k  ,i,j)*rw1(k  ,i,j))
                  end do
               end do
            end do

         end if

!        end of small time step

      end do
!
!     add large time step rw back into rw1 
!
      do j=1,ny
         jv1 = j-1
         if(jper*j.eq.1)  jv1 = ny1
         jm1 = max(j-1,1)
         if(jper*j.eq.1)  jm1 = ny1
         jp1 = min(j+1,ny)
         if(jper*j.eq.ny)  jp1=2
         do i=1,nx
            if(mod(i,2).eq.0)  then
               jpj=j
               jpm=jm1
            else
               jpj=jp1
               jpm=j
            end if
            im1 = i-1
            if(iper*i.eq.1)  im1 = nx1
            do k=2,nz1
               rw1(k,i,j) = rw(k,i,j) + 1./hh(i,j)*(rw1(k  ,i  ,j)  &
     &                    - .25*(gu1(k  ,i  ,j  )*ru11(k  ,i  ,j  )  &
     &                          +gu1(k  ,im1,jpj)*ru11(k  ,im1,jpj)  &
     &                          +gu1(k-1,i  ,j  )*ru11(k-1,i  ,j )  &
     &                          +gu1(k-1,im1,jpj)*ru11(k-1,im1,jpj))  &
     &                    - .25*(gu3(k  ,i  ,j  )*ru31(k  ,i  ,j )  &
     &                          +gu3(k  ,im1,jpm)*ru31(k  ,im1,jpm)  &
     &                          +gu3(k-1,i  ,j  )*ru31(k-1,i  ,j )  &
     &                          +gu3(k-1,im1,jpm)*ru31(k-1,im1,jpm))  &
     &                    - .25*(gu2(k  ,i  ,j  )*ru21(k  ,i  ,j  )  &
     &                          +gu2(k  ,i  ,jv1)*ru21(k  ,i  ,jv1)  &
     &                          +gu2(k-1,i  ,j  )*ru21(k-1,i  ,j  )  &
     &                          +gu2(k-1,i  ,jv1)*ru21(k-1,  i,jv1)))
            end do
         end do
      end do
      do j=1,ny
         do i=iper,nx
            do k=1,nz1
               ru11(k,i,j) = ru11(k,i,j)+ru1(k,i,j)
               ru31(k,i,j) = ru31(k,i,j)+ru3(k,i,j)
            end do
         end do
      end do
      do i=1,nx
         do j=jper,ny
            do k=1,nz1
               ru21(k,i,j) = ru21(k,i,j)+ru2(k,i,j)
            end do
         end do
      end do

      return
      end
