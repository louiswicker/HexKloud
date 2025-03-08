!                                                                     72
!-----------------------------------------------------------------------
!
      subroutine kessler_joe( t1t, qv1t, qc1, qr1,  &
     &                        rho, pii, dt, dz, nz1, ny, nx         )

      implicit none
      integer nx, ny, nz1
      real dt, dz
      real t1t (nz1,nx,ny), qv1t(nz1,nx,ny),    &
     &     qc1 (nz1,nx,ny), qr1 (nz1,nx,ny),  &
     &     rho (nz1,nx,ny), pii (nz1,nx,ny)
      integer mz
      parameter( mz=200 )
      real qrprod(mz), prod (mz), rcgs( mz), rcgsi (mz),qrr  &
     &    ,ern   (mz), vt   (mz), vtden(mz), gam   (mz)  &
     &    ,r     (mz), rhalf(mz), velqr(mz), buoycy(mz)  &
     &    ,pk    (mz), pc   (mz), f0   (mz), qvs   (mz)
!-----------------------------------------------------------------------
      real c1, c2, c3, c4, f5, mxfall, dtfall, fudge
      real cp, rdz, product, ackess, artemp, artot, ckess, dtl, dtlc
      real f2x, fvel, psl, veld, velu, xk, xki
      integer nfall
      integer  i,j,k,n,nz2

      ackess = 0.001
      ckess  = 2.2
      fvel   = 36.34
      f2x    = 17.27
      f5     = 237.3*F2X*2.5E6/1003.
      XK     = .2875          
      XKI    = 1./XK         
      psl    = 1000.

      do k=1,nz1
         r(k)     = 0.001*rho(k,1,1)
         rhalf(k) = sqrt(rho(1,1,1)/rho(k,1,1))
         pk(k)    = pii(k,1,1)
         pc(k)    = 3.8/(pk(k)**xki*psl)
         f0(k)    = 2.5E6/(1003.*pk(k))
      enddo
!
      dtl  = dt
      dtlc = dt
      nz2  = nz1 - 1

      do j=1,ny
         do i=1,nx
!
            DO K=1,NZ1
               QRPROD(K) = qc1(K,i,j)  &
     &               -(qc1(K,i,j)-DTL*AMAX1(ACKESS*(QC1(K,i,j)-.001),  &
     &                     0.))/(1.+DTL*CKESS*QR1(K,i,j)**.875)
!cc???         VELQR(K)  = (QR(K,i,j)*R(K))**1.1364*RHALF(K)
               VELQR(K)  = (QR1(K,i,j)*R(K))**1.1364*RHALF(K)
               QVS(K)    = PC(K)*EXP(F2X*(PK(K)*T1T(K,i,j)-273.)  &
     &                            /(PK(K)*T1T(K,i,j)- 36.))
            END DO
            VELU        = (QR1(2,i,j)*R(2))**1.1364*RHALF(2)
            VELD        = (QR1(1,i,j)*R(1))**1.1364*RHALF(1)
            qr1(1,i,j) = qr1(1,i,j)+DTL*(VELU-VELD)*FVEL/(R(1)*DZ)
            DO K=2,NZ2
               qr1(K,i,j) =qr1(K,i,j)+DTL*FVEL*(VELQR(K+1)-VELQR(K-1))  &
     &                                 /(R(K  )*DZ*(1.+1.))
!c             qr1(K,i,j) = qr1(K,i)+DTL*FVEL*(VELQR(K+1)-VELQR(K))
!c   &                                 /(R(K  )*DZ)
            END DO
            qr1(NZ1,i,j)  = qr1(NZ1,i,j)-DTL*FVEL*VELQR(NZ2)  &
     &                                 /(R(NZ1)*DZ*(1.+1.))
            ARTEMP     = 36340.*(.5*(VELQR(2)+VELQR(1))+VELD-VELU)
           ! ARTOT      = ARTOT+DTD*ARTEMP ! dtd is not defined, but artot is not used
            DO K=1,NZ1
               qc1(K,i,j) = AMAX1(qc1(K,i,j)-QRPROD(K),0.)
               qr1(K,i,j) = AMAX1(qr1(K,i,j)+QRPROD(K),0.)
               PROD(K)     = (QV1T(K,i,j)-QVS(K))/(1.+QVS(K)*F5  &
     &                             /(PK(K)*T1T(K,i,j)-36.)**2)
            END DO
            DO K=1,NZ1
               ERN(K)=AMIN1(DTLC*(((1.6+124.9*(R(K)*qr1(K,i,j))**.2046)  &
     &               *(R(K)*qr1(K,i,j))**.525)/(2.55E6*PC(K)  &
     &               /(3.8 *QVS(K))+5.4E5))*(DIM(QVS(K),QV1T(K,i,j))  &
     &               /(R(K)*QVS(K))),  &
     &                AMAX1(-PROD(K)-qc1(K,i,j),0.),qr1(K,i,j))
            END DO
            DO K=1,NZ1
               BUOYCY(K)   = F0(K)*(AMAX1(PROD(K),-qc1(K,i,j))-ERN(K))
               QV1T(K,i,j) = AMAX1(QV1T(K,i,j)  &
     &                        -AMAX1(PROD(K),-qc1(K,i,j))+ERN(K),0.)
               qc1(K,i,j) = qc1(K,i,j)+AMAX1(PROD(K),-qc1(K,i,j))
               qr1(K,i,j) = qr1(K,i,j)-ERN(K)
               T1T (K,i,j) = T1T (K,i,j)+BUOYCY(K)
            END DO
         end do
      end do

      RETURN           
      END              
!
!***********************************************************************
!
