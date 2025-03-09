 MODULE WENO

    implicit none

!
    !-----------------------------------------------------!
    !       WENO only:                                    !
    ! formulation for weno "smoothness indicators"        !
    !  1 = original (eg, Jiang and Shu, 1996, JCP)
    !  2 = Borges et al. (2008, JCP)

    integer, parameter :: siform = 2

    real(kind=8), parameter :: weps = 1.0d-20

    !-----------------------------------------------------!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 CONTAINS

 !dir$ attributes forceinline :: weno3

      REAL FUNCTION WENO3(s1,s2,s3)

      implicit none

      real(kind=4), intent(in) :: s1,s2,s3

      real(kind=8) :: b1,b2
      real(kind=8) :: w1,w2

      ! 3rd-order weighted essentially non-oscillatory (weno)
      ! Jiang and Shu, 1996, JCP

      b1 = (s1-s2)**2
      b2 = (s2-s3)**2

      if( siform.eq.1 )then
        ! original WENO (eg, Jiang and Shu, 1996, JCP)
        w1 = (1.0/3.0)/(weps+b1)**2
        w2 = (2.0/3.0)/(weps+b2)**2
      elseif( siform.eq.2 )then
        ! improved smoothness indicators (Borges et al, 2008, JCP)
        w1 = (1.0/3.0)*(1.0+min(1.0d30,abs(b1-b2)/(b1+weps))**2)
        w2 = (2.0/3.0)*(1.0+min(1.0d30,abs(b1-b2)/(b2+weps))**2)
      endif

      weno3 = ( w1*( (-1.0/2.0)*s1 + ( 3.0/2.0)*s2 )  &
               +w2*( ( 1.0/2.0)*s2 + ( 1.0/2.0)*s3 )  &
              )/( w1+w2 )

      END FUNCTION WENO3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 !dir$ attributes forceinline :: weno5

      REAL FUNCTION WENO5(s1,s2,s3,s4,s5)

      implicit none

      real(kind=4), intent(in) :: s1,s2,s3,s4,s5

      real(kind=8) :: b1,b2,b3
      real(kind=8) :: w1,w2,w3

      ! 5th-order weighted essentially non-oscillatory (weno)
      ! Jiang and Shu, 1996, JCP

      b1 = (13.0/12.0)*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
      b2 = (13.0/12.0)*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
      b3 = (13.0/12.0)*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

      if( siform.eq.1 )then
        ! original WENO (eg, Jiang and Shu, 1996, JCP)
        w1 = 0.1/(weps+b1)**2
        w2 = 0.6/(weps+b2)**2
        w3 = 0.3/(weps+b3)**2
      elseif( siform.eq.2 )then
        ! improved smoothness indicators (Borges et al, 2008, JCP)
        w1 = 0.1*(1.0+min(1.0d30,abs(b1-b3)/(b1+weps))**2)
        w2 = 0.6*(1.0+min(1.0d30,abs(b1-b3)/(b2+weps))**2)
        w3 = 0.3*(1.0+min(1.0d30,abs(b1-b3)/(b3+weps))**2)
   !    w1 = 0.1
   !    w2 = 0.6
   !    w3 = 0.3
      endif

      weno5 = ( w1*( ( 2.0/6.0)*s1 + (-7.0/6.0)*s2 + (11.0/6.0)*s3 )  &
               +w2*( (-1.0/6.0)*s2 + ( 5.0/6.0)*s3 + ( 2.0/6.0)*s4 )  &
               +w3*( ( 2.0/6.0)*s3 + ( 5.0/6.0)*s4 + (-1.0/6.0)*s5 )  &
              )/( w1+w2+w3 )

      END FUNCTION WENO5

END MODULE WENO
