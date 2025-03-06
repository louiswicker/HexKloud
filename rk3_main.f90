!===============================================================================
!********** nonlinear nonhydrostatic time-split model on hexagonal grid
!********** mass conserving, flux-form variables
!********** open or periodic boundaries.
!********** transformed terrain-following vertical coordinate.
!********** includes atmospheric mean state. constant stability
!********** reference code for terrain-folling height coordinate model
!
!  This code incorporates the wicker and skamarock RK3 timesplitting
!  algorithm
!                                                                     
!  This code also includes kessler microphysics
!  This code also includes NSSL-2M scheme microphysics (Ted Mansell)
!===============================================================================

      PROGRAM HEXKLOUD   ! In honor of Joe Klemp

      USE module_mp_nssl_2mom

      implicit none

      integer, parameter :: nz= 41, nx= 91, ny= 79, nz1=nz-1, nx1=nx-1, ny1=ny-1

!     parameter (nz= 41, nx= 61, ny= 53, nz1=nz-1, nx1=nx-1, ny1=ny-1)
!     parameter (nz= 41, nx=181, ny=157, nz1=nz-1, nx1=nx-1, ny1=ny-1)
!     parameter (nz= 41, nx= 47, ny= 40, nz1=nz-1, nx1=nx-1, ny1=ny-1)
!     parameter (nz= 41, nx= 121, ny=105, nz1=nz-1, nx1=nx-1, ny1=ny-1)
!     parameter (nz= 41, nx= 181, ny= 53, nz1=nz-1, nx1=nx-1, ny1=ny-1)
!     parameter (nz= 41, nx= 101, ny= 5, nz1=nz-1, nx1=nx-1, ny1=ny-1)
!     parameter (nz= 41, nx= 5, ny=101, nz1=nz-1, nx1=nx-1, ny1=ny-1)

      real u1  (nz1,0:nx,ny), u11 (nz1,0:nx,ny), ru1 (nz1,0:nx,ny)  &
     &    ,ru11(nz1,0:nx,ny), fu1 (nz1,0:nx,ny)  &
     &    ,ru1i(nz1,0:nx,ny), gu1 (nz1,0:nx,ny), du1 (nz1,0:nx,ny)  &
     &    ,u3  (nz1,0:nx,ny), u31 (nz1,0:nx,ny), ru3 (nz1,0:nx,ny)  &
     &    ,ru31(nz1,0:nx,ny), fu3 (nz1,0:nx,ny)  &
     &    ,ru3i(nz1,0:nx,ny), gu3 (nz1,0:nx,ny), du3 (nz1,0:nx,ny)  &
     &    ,u2  (nz1,nx,0:ny), u21 (nz1,nx,0:ny), ru2 (nz1,nx,0:ny)  &
     &    ,ru21(nz1,nx,0:ny), fu2 (nz1,nx,0:ny)  &
     &    ,ru2i(nz1,nx,0:ny), gu2 (nz1,nx,0:ny), du2 (nz1,nx,0:ny)  &
     &    ,w   (nz ,nx,  ny), w1  (nz ,nx  ,ny), rw  (nz ,nx  ,ny)  &
     &    ,rw1 (nz ,nx  ,ny), fw  (nz1,nx  ,ny)  &
     &    ,t   (nz1,nx  ,ny), t1  (nz1,nx  ,ny), rt  (nz1,nx  ,ny)  &
     &    ,rt1 (nz1,nx  ,ny), ft  (nz1,nx  ,ny)  &
     &    ,ti  (nz1,nx  ,ny), rti (nz1,nx  ,ny), rtb (nz1,nx  ,ny)  &
     &    ,rr  (nz1,nx  ,ny), rr1 (nz1,nx  ,ny)  &
     &    ,rri (nz1,nx  ,ny), rb  (nz1,nx  ,ny), fr  (nz1,nx  ,ny)  &
     &    ,p   (nz1,nx  ,ny), pb  (nz1,nx  ,ny), pii (nz1,nx  ,ny)  &
     &    ,ww  (nz ,nx  ,ny), rho (nz1,nx  ,ny), tb  (nz1,nx  ,ny)  &
     &    ,rs  (nz1,nx  ,ny), ts  (nz1,nx  ,ny), div (nz1,nx  ,ny)  &
     &    ,a   (nz1,nx  ,ny), b   (nz1,nx  ,ny), c   (nz1,nx  ,ny)  &
     &    ,cofwz (nz1,nx,ny), coftz (nz ,nx,ny), cofwt (nz1,nx,ny)  &
     &    ,alpha (nz1,nx,ny), gamma (nz1,nx,ny)                     &
     &    ,flux1(nz1,0:nx,ny),flux2(nz1,nx,0:ny),flux3(nz1,0:nx,ny)  &
     &    ,cofwrr(nz1,nx,ny), cofwr     (nx,ny), dhh1      (nx,ny)  &
     &    ,dhh2      (nx,ny), dhh3      (nx,ny), hh        (nx,ny)  &
     &    ,hs        (nx,ny), wdtz(nz), zu(nz1), zw(nz),   ds(nz1)  &
     &    , u1z(nz1),u2z(nz1), u3z(nz1), tz(nz1), fluxz(0:nz1,nx,ny)  &
     &    ,ax(nz), tzv (nz1), rqvb(nz1),  rel_hum (nz1), qvzv(nz1)

      real ru1_save (nz1,0:nx,ny), ru3_save (nz1,0:nx,ny)  &
     &    ,ru2_save (nz1,nx,0:ny), rw_save  (nz ,  nx,ny)  &
     &    ,rt_save  (nz1,nx  ,ny), rr_save  (nz1,  nx,ny)  &
     &    ,t_d_tend (nz1,nx  ,ny)

! rqx,rqx1 are mass content (qx*rho)
! qx, qx1 are mass mixing ratio
! fqx are tendencies?
! Now 4D to accomadate NSSL-2M scheme

      real, allocatable :: rqx(:,:,:,:),rqx1(:,:,:,:), qx(:,:,:,:), qx1(:,:,:,:), fqx(:,:,:,:)
      real, allocatable :: rsx(:,:,:,:),rsx1(:,:,:,:), sx(:,:,:,:), sx1(:,:,:,:), fsx(:,:,:,:)
      real, allocatable :: dz3d(:,:,:), dbz(:,:,:),ws(:,:,:), pres(:,:,:)
      real, allocatable :: rainnc(:,:), rainncv(:,:)

      integer :: nmoist, nscalar

      real*4 plt (nx,ny), pltx(nx,nz), plty(ny,nz), hxpl  (nx)  &
     &      ,xh  (nx,ny), xu1 (nx,ny), xu2 (nx,ny), xu3(nx,ny)  &
     &      ,yh  (nx,ny), yu1 (nx,ny), yu2 (nx,ny), yu3(nx,ny)  &
     &      ,x(nx),y(ny),time,pxl,pxr,pyl,pyr,pzl,zptop,xpll,xplr  &
     &      ,ypll,yplr,zplb,zplt,dxp,dyp,dzp,wmax(2401),waxis(2401)  &
     &      ,wmplt

      real*4 Azero(1)

      common /grid/ xh, xu1, xu2, xu3, yh, yu1, yu2, yu3

      integer imass, rk_step, ns_rk, total_steps
      character*3 slice(2)
      character*6 plane
      equivalence (plane,slice)

      integer, PARAMETER :: IERF=6,LUNI=2,IWID=1  
      integer :: IWTY = 20 !  1=ncgm/gmeta; 20=PostScript
      real :: ampl, angle, area, c1f, c2, c2f, cb
      real :: cofrz, cp, cti, d, dtl, dts, dtsa, dtsd
      real :: dtseps, dtsf, dtsg
      real :: dx, dy, dz
      real :: epssm, f, fac1, fura, g, hm
      integer :: i, ii, im1, ip, ip1, iper, ipf, ipi, iplt, ipp, n
      integer :: itr, ittrm, iwmax
      integer :: j, j1, jj, jm1, jn, jp1, jper, jpf, jpi, jpj, jpm, jpp, jv1, jwmax
      integer :: k, kk, kkk, km1, kwmax, nit, npl, npr, ns, ns0
      integer :: nxc, nxpl, nyc, nypl, nz2, nzpl
      real :: p0, pi, pitop, pressure, qvs
      real :: r, rad, radx, rady, radz, rcv, rd
      real :: rdx, rdy, rdz
      real :: resm, ritot, rrtot, rtot, rttop, rula
      real :: side, smdiv, smdivx, smdivz, sum
      real :: t0, tdiff, temp, thetak, tinit, tip, tk, tkm1, tkp1
      real :: tmax, tstp, u1m, u2m, u3m, um, ub, ur, vm, vnu, xa, xc
      real :: xht, xl, xn, xn2, xn2l, xn2m, xnu, xnus, xnus0, xnusz, xnusz0, xnut
      real :: ya, yc, yl, yht
      real :: zcent, zd, zinv, zt, ztemp
      integer, parameter :: lv = 1, lc = 2, lr = 3
      integer :: li = 4, ls = 5, lh = 6, lhl = 7
      integer :: lnc = 1, lnr = 2, lni = 3, lns = 4, lnh = 5, lnhl = 6, lccn = 7
      integer :: lvh = 8, lvhl = 9
      real    :: tmp
      real, dimension(20) :: nssl_params

      integer :: IDS=1,IDE=nx, JDS=1,JDE=ny, KDS=1,KDE=nz1, &
                 IMS=1,IME=nx, JMS=1,JME=ny, KMS=1,KME=nz1, &
                 ITS=1,ITE=nx, JTS=1,JTE=ny, KTS=1,KTE=nz1

      real    :: nssl_cccn = 6.e8, nssl_alphah=0, nssl_alphahl=1,  &
                 nssl_cnoh=4.e4, nssl_cnohl=4.e3, nssl_cnor=8.e6, nssl_cnos=3.0e6, &
                 nssl_rho_qh=600., nssl_rho_qhl=800., nssl_rho_qs=100.

      integer           :: nssl_ccn_is_ccna=1, nssl_2moment_on=1
      integer           :: mp_physics = 1 ! microphysics: 1=kessler; 18= NSSL 2-moment
      integer           :: iadvord = 5 ! advection order
      character(len=6)  :: order

      real              :: delt  = 3.     ! bubble temp
      real              :: dt    = 6.0    ! time step
      logical           :: debug = .false.

! Arrays for netCDF

      character(len=3),  dimension(20) :: varlabel
      character(len=20)                :: ncdf_file

      real, allocatable :: ncdf_var(:,:,:,:)

! Namelist declarations

      character(LEN=50) :: filename = 'namelist.input'
      logical           :: if_exist
      integer           :: iunit

      namelist /main/ mp_physics, iadvord, nssl_2moment_on, nssl_cccn, delt, dt, iwty, debug

! Start here and read namelist

      INQUIRE(file=trim(filename), exist=if_exist)

      if (  if_exist ) then
  
       iunit = 15
       open(15,file=trim(filename),status='old',form='formatted')
       rewind(15)
       read(15,NML=main)
      endif

      if ( mp_physics == 1 ) then
        nmoist = 3
        nscalar = 0
        allocate( dz3d(1,1,1), dbz(1,1,1), ws(1,1,1), pres(1,1,1) )
        dz3d(:,:,:) = dz
      elseif ( mp_physics == 18 ) then
        nmoist = 7
         if ( nssl_2moment_on == 1 ) then
            i = 5
            nscalar = 9
         elseif ( nssl_2moment_on == 0 ) then
            i = 0
            nscalar = 1
            lnc = 1; lnr = 1; lni = 1; lns = 1; lnh = 1; lnhl = 1; lccn = 1
            lvh = 1; lvhl = 1
            nssl_ccn_is_ccna = 0
         endif
         
        allocate( dz3d(nz1,nx,ny), dbz(nz1,nx,ny), ws(nz1,nx,ny), pres(nz1,nx,ny) )
        allocate( rainnc(nx,ny), rainncv(nx,ny) )

! call init?
       nssl_params(:)  = 0
       nssl_params(1)  = nssl_cccn
       nssl_params(2)  = nssl_alphah
       nssl_params(3)  = nssl_alphahl
       nssl_params(4)  = nssl_cnoh
       nssl_params(5)  = nssl_cnohl
       nssl_params(6)  = nssl_cnor
       nssl_params(7)  = nssl_cnos
       nssl_params(8)  = nssl_rho_qh
       nssl_params(9)  = nssl_rho_qhl
       nssl_params(10) = nssl_rho_qs
       nssl_params(11) = 0 ! nssl_ipelec_tmp
       nssl_params(12) = 0 ! config_flags%nssl_isaund
       nssl_params(13) = 0 ! reserved
       nssl_params(14) = 0 ! reserved
       nssl_params(15) = 0 ! reserved
       CALL nssl_2mom_init(nssl_params=nssl_params,ipctmp=i,mixphase=0,        &
                           nssl_density_on=.true.,                             &
                           nssl_hail_on=.true.,                                &
                           nssl_ccn_on= ( i >= 5 ),                            &
                           nssl_icdx=6,                                        &
                           nssl_icdxhl=6,                                      &
                           ccn_is_ccna=nssl_ccn_is_ccna)
      else
        write(0,*) 'unsupported value of mp_physics: ', mp_physics
        stop
      endif
      
      if ( iadvord == 2 ) then
        order = 'second'
      elseif ( iadvord == 3 ) then
        order = 'third '
      elseif ( iadvord == 4 ) then
        order = 'fourth'
      elseif ( iadvord == 5 ) then
        order = 'fifth '
      elseif ( iadvord == 6 ) then
        order = 'sixth '
      else
        write(0,*) 'invalid value of iadvord: ',iadvord, 'resetting to 5'
        iadvord = 5
        order = 'fifth '
      endif
      
      allocate( rqx(nz1,nx,ny,nmoist),  &
                rqx1(nz1,nx,ny,nmoist), &
                qx(nz1,nx,ny,nmoist),   &
                qx1(nz1,nx,ny,nmoist),  &
                fqx(nz1,nx,ny,nmoist) )

      if ( nscalar > 0 ) then
        allocate( rsx(nz1,nx,ny,nscalar),  &
                 rsx1(nz1,nx,ny,nscalar), &
                 sx(nz1,nx,ny,nscalar),   &
                 sx1(nz1,nx,ny,nscalar),  &
                 fsx(nz1,nx,ny,nscalar) )
      else
        allocate( rsx(1,1,1,1),  &
                 rsx1(1,1,1,1), &
                 sx(1,1,1,1),   &
                 sx1(1,1,1,1),  &
                 fsx(1,1,1,1) )
      endif

!--------------
!
      include "initialize.inc.f90"
!
!--------------

      dz3d(:,:,:) = dz

      Azero(1) = 0.0

!***********************************************************************
!     timestep loop
!***********************************************************************
!
!*****Large time step calculations
!
      kkk = ip
      do nit = 0, total_steps

      if(nit .ne. 0) then

      kkk   = kkk+1
      npr   = npr+1
      time  = nit*dt
      tinit =.05*xa
      if(npr.eq.1)  then
         write(6,*) 't,wmax= ',time, wmax(nit)
!        write(0,*) 't,wmax= ',time, wmax(nit)
         npr=0
      end if

      do j=1,ny
         do i=0,nx
            do k=1,nz1
               ru1_save(k,i,j) = ru1(k,i,j)
               ru3_save(k,i,j) = ru3(k,i,j)
            end do
         end do
      end do
      do j=0,ny
         do i=1,nx
            do k=1,nz1
               ru2_save(k,i,j) = ru2(k,i,j)
            end do
         end do
      end do
      do j=1,ny
         do i=1,nx
            do k=1,nz
               rw_save(k,i,j) = rw(k,i,j)
            end do
            do k=1,nz1
               rt_save(k,i,j) = rt(k,i,j)
               rr_save(k,i,j) = rr(k,i,j)
            end do
         end do
      end do
!
!*****Beginning of Runge Kutta time steps
!
      do rk_step = 1,3
!**********
!      do rk_step = 3,3

      if(rk_step .eq. 1) then
!         ns_rk = ns0/3
         ns_rk = 1
         dts = dt/float(ns_rk)/3.
      else if(rk_step .eq. 2) then
         ns_rk = ns0/2
         dts = dt/float(ns_rk)/2.
      else if(rk_step .eq. 3) then
         ns_rk = ns0
         dts = dt/float(ns_rk)
      end if
!**********
!      ns_rk = 1
!     dts = dt
!**********
      dtseps = .25*dts*(1.+epssm)
      cofrz = 2.*dtseps*rdz
      dtsa   = dts*side/area
      dtsd   = dts/d
      dtsf   = dts*sqrt(3.)/6.*f
      dtsg   = dts*.5/d*g
      xnus   = dts*xnus0
      xnusz  = dts*xnusz0
      do j=1,ny
         do i=1,nx
            cofwr(i,j) = dtseps*g*hh(i,j)
         end do
      end do
!
!        calculation of omega, ww = gu1*ru1 + gu2*ru2 + gu3*ru3 + hh* rw
!
         do j=1,ny
            jm1 = j-1
            if(jper*j.eq.1) jm1=ny1
            do i=1,nx
               if(mod(i,2).eq.0)  then
                  jpj=j
                  jpm=jm1
               else
                  jpj=jp1
                  jpm=j
               end if
               im1 = i-1
               if(iper*i.eq.1) im1=nx1
               do k=2,nz1
                  ww(k,i,j) = .25*(gu1(k  ,i  ,j  )*ru1(k  ,i  ,j  )  &
     &                            +gu1(k  ,im1,jpj)*ru1(k  ,im1,jpj)  &
     &                            +gu1(k-1,i  ,j  )*ru1(k-1,i  ,j  )  &
     &                            +gu1(k-1,im1,jpj)*ru1(k-1,im1,jpj))  &
     &                       +.25*(gu3(k  ,i  ,j  )*ru3(k  ,i  ,j  )  &
     &                            +gu3(k  ,im1,jpm)*ru3(k  ,im1,jpm)  &
     &                            +gu3(k-1,i  ,j  )*ru3(k-1,i  ,j  )  &
     &                            +gu3(k-1,im1,jpm)*ru3(k-1,im1,jpm))  &
     &                       +.25*(gu2(k  ,i  ,j  )*ru2(k  ,i  ,j  )  &
     &                            +gu2(k  ,i  ,jm1)*ru2(k  ,i  ,jm1)  &
     &                            +gu2(k-1,i  ,j  )*ru2(k-1,i  ,j  )  &
     &                            +gu2(k-1,i  ,jm1)*ru2(k-1,i  ,jm1))  &
     &                        +hh(i,j)*rw(k,i,j)
               end do
            end do
         end do

         call rhs_u1(u1,u11,ru1,fu1,ww,rho,ru2,ru3,u1z,u2z,u3z,u1m,u2m,  &
     &             u3m,ds,dtsa,dtsd,dtsf,dts,c1f,c2f,rdz,xnus,xnusz,  &
     &             nz1,nx,ny,iper,jper,flux1,flux2,flux3,fluxz,order)

         call rhs_u3(u3,u31,ru3,fu3,ww,rho,ru1,ru2,u1z,u2z,u3z,u1m,u2m,  &
     &             u3m,ds,dtsa,dtsd,dtsf,dts,c1f,c2f,rdz,xnus,xnusz,  &
     &             nz1,nx,ny,iper,jper,flux1,flux2,flux3,fluxz,order)

         call rhs_u2(u2,u21,ru2,fu2,ww,rho,ru1,ru3,u1z,u2z,u3z,u1m,u2m,  &
     &             u3m,ds,dtsa,dtsd,dtsf,dts,c1f,c2f,rdz,xnus,xnusz,  &
     &             nz1,nx,ny,iper,jper,flux1,flux2,flux3,fluxz,order)

!         call rhs_w( w,w1,fw,ww,p,pb,rt,rtb,rho,ru1,ru2,ru3,rcv,rb,rqv,  &
!     &               rqc,rqr,rqvb,dtsa,g,ds,dts,rdz,f,xnus,xnusz,nz1,  &
!     &               nx,ny,iper,jper,flux1,flux2,flux3,fluxz,order)
         call rhs_w( w,w1,fw,ww,p,pb,rt,rtb,rho,ru1,ru2,ru3,rcv,rb,rqx,  &
     &               nmoist,rqvb,dtsa,g,ds,dts,rdz,f,xnus,xnusz,nz1,  &
     &               nx,ny,iper,jper,flux1,flux2,flux3,fluxz,order)

         call rhs_s( t ,t1 ,ft ,ww,ru1,ru2,ru3,rho,ds,dts,dtsa,rdz,  &
     &               xnus,xnusz,nz1,nx,ny,iper,jper,  &
     &               ti,nz1,nx,ny,flux1,flux2,flux3,fluxz,order)

! qv
         call rhs_s( qx(1,1,1,lv),qx1(1,1,1,lv),fqx(1,1,1,lv),ww,ru1,ru2,ru3,rho,ds,dts,dtsa,rdz,  &
     &               xnus,xnusz,nz1,nx,ny,iper,jper,  &
     &               qvzv,nz1,1,1,flux1,flux2,flux3,fluxz,order)

! other mixing ratios
         do n = 2,nmoist
           call rhs_s( qx(1,1,1,n),qx1(1,1,1,n),fqx(1,1,1,n),ww,ru1,ru2,ru3,rho,ds,dts,dtsa,rdz,  &
     &               xnus,xnusz,nz1,nx,ny,iper,jper,  &
     &               Azero, 1  ,1,1,flux1,flux2,flux3,fluxz,order)
         enddo
!          call rhs_s( qc,qc1,fqc,ww,ru1,ru2,ru3,rho,ds,dts,dtsa,rdz,  &
!      &               xnus,xnusz,nz1,nx,ny,iper,jper,  &
!      &               Azero, 1  ,1,1,flux1,flux2,flux3,fluxz,order)
! 
!          call rhs_s( qr,qr1,fqr,ww,ru1,ru2,ru3,rho,ds,dts,dtsa,rdz,  &
!      &               xnus,xnusz,nz1,nx,ny,iper,jper,  &
!      &               Azero, 1  ,1,1,flux1,flux2,flux3,fluxz,order)


! other scalars
         do n = 1,nscalar
           call rhs_s( sx(1,1,1,n),sx1(1,1,1,n),fsx(1,1,1,n),ww,ru1,ru2,ru3,rho,ds,dts,dtsa,rdz,  &
     &               xnus,xnusz,nz1,nx,ny,iper,jper,  &
     &               Azero, 1  ,1,1,flux1,flux2,flux3,fluxz,order)
         enddo

         call rhs_rho( fr,ru1,ru2,ru3,ww,dts,dtsa,rdz,  &
     &                 nz1,nx,ny,iper,jper      )

!
!--------------

      include "boundaries.inc.f90"

!--------------
!
!        advance moisture variables over interval
!
         do j=1,ny
            do i=1,nx
               do k=1,nz1
                  do n = 1,nmoist
                    rqx(k,i,j,n) = amax1(rqx1(k,i,j,n) + ns_rk*fqx(k,i,j,n),0.0)
                  enddo
                  do n = 1,nscalar
                    rsx(k,i,j,n) = amax1(rsx1(k,i,j,n) + ns_rk*fsx(k,i,j,n),0.0)
                  enddo
!                   rqv(k,i,j) = amax1(rqv1(k,i,j) + ns_rk*fqv(k,i,j),0.0)
!                   rqc(k,i,j) = amax1(rqc1(k,i,j) + ns_rk*fqc(k,i,j),0.0)
!                   rqr(k,i,j) = amax1(rqr1(k,i,j) + ns_rk*fqr(k,i,j),0.0)
               end do
            end do
         end do
!
!        add in the diabatic theta_v tendency from the last timestep
!
         do j=1,ny
            do i=1,nx
               do k=1,nz1
                  ft(k,i,j) = ft(k,i,j) + dts*rho(k,i,j)*t_d_tend(k,i,j)
               end do
            end do
         end do
!
!********small time step calculations
!
!        coefficients for tri-diagonal matrix and do small steps
!
         call calc_scoef( dtseps, c2, hh, rdz, t, p, tb,   &
     &                       rho, rb, rqx, nmoist, rqvb,  &
     &                       g, rcv, cofwz, coftz, cofwt,  &
     &                       cofwr, cofwrr, cofrz,  &
     &                       a, b, c, alpha, gamma, nx,ny,nz1    )
!
         call smlstep( ru1,ru11,fu1,ru2,ru21,fu2,ru3,ru31,fu3,   &
     &                    rw,rw1,fw,t,ts,rt1,ft,rs,rr1,fr,  &
     &                    p,ww,div,du1,du2,du3,hh,gu1,gu2,gu3,  &
     &                    a,alpha,gamma,dhh1,dhh2,dhh3,  &
     &                    cofrz,coftz,cofwz,cofwr,cofwrr,cofwt,  &
     &                    rdz,dts,dtsa,dtsd,c2,smdivx,smdivz,resm,  &
     &                    nx,ny,nz1,ns_rk,iper,jper  )

         if(rk_step .le. 2) then
!
!********set levels for full step
!
            do j=1,ny
               do i=0,nx
                  do k=1,nz1
                     ru1 (k,i,j) = ru11    (k,i,j)
                     ru11(k,i,j) = ru1_save(k,i,j)
                     ru3 (k,i,j) = ru31    (k,i,j)
                     ru31(k,i,j) = ru3_save(k,i,j)
                  end do
               end do
            end do
            do j=0,ny
               do i=1,nx
                  do k=1,nz1
                     ru2 (k,i,j) = ru21    (k,i,j)
                     ru21(k,i,j) = ru2_save(k,i,j)
                  end do
               end do
            end do
            do j=1,ny
               do i=1,nx
                  do k=1,nz
                     rw (k,i,j) = rw1    (k,i,j)
                     rw1(k,i,j) = rw_save(k,i,j)
                  end do
                  do k=1,nz1
                     rt (k,i,j) = rt1    (k,i,j)
                     rt1(k,i,j) = rt_save(k,i,j)
                     rr (k,i,j) = rr1    (k,i,j)
                     rr1(k,i,j) = rr_save(k,i,j)
                     rho(k,i,j) = rb(k,i,j) + rr(k,i,j)
                  end do
               end do
            end do

            do j=1,ny
               jp1 = min(j+1,ny)
               if(jper*j.eq.ny)  jp1 = 2
               jv1 = j-1
               if(jper*j.eq.1 )  jv1 = ny1
               jm1 = max(j-1,1)
               if(jper*j.eq.1 )  jm1 = ny1
               do i=1,nx
                  if(mod(i,2).eq.0)  then
                     jpj = j
                     jpm = jm1
                  else
                     jpj = jp1
                     jpm = j
                  end if
                  ip1 = min(i+1,nx)
                  if(iper*i.eq.nx)  ip1 = 2
                  im1 = i-1
                  if(iper*i.eq.1 )  im1 = nx1
                  do k=1,nz1
                     u1(k,i,j) = ru1(k,i,j)  &
     &                              /(.5*(rho(k,i,j)+rho(k,ip1,jpm)))
                     u3(k,i,j) = ru3(k,i,j)  &
     &                              /(.5*(rho(k,i,j)+rho(k,ip1,jpj)))
                     u2(k,i,j) = ru2(k,i,j)  &
     &                              /(.5*(rho(k,i,j)+rho(k,i  ,jp1)))
                  end do

                  do k=1,nz1
                     t (k,i,j) = (rtb(k,i,j)+rt(k,i,j))/rho(k,i,j)
                     p (k,i,j) = (hh(i,j)*(rtb(k,i,j)+rt(k,i,j))  &
     &                                   /(100000./287./300.) )**rcv
                     do n = 1,nmoist
                       qx(k,i,j,n) = rqx(k,i,j,n)/rho(k,i,j)
                     enddo

                     if ( nscalar > 1 ) then
                     do n = 1,nscalar
                       sx(k,i,j,n) = rsx(k,i,j,n)/rho(k,i,j)
                     enddo
                     endif

!                      qv(k,i,j) = rqv(k,i,j)/rho(k,i,j)
!                      qc(k,i,j) = rqc(k,i,j)/rho(k,i,j)
!                      qr(k,i,j) = rqr(k,i,j)/rho(k,i,j)

                  end do

                  k=1
                     rw(k,i,j) = - .25/hh(i,j)*  &
     &                         (3.*(gu1(k  ,i  ,j  )*ru1(k  ,i  ,j  )  &
     &                             +gu1(k  ,im1,jpj)*ru1(k  ,im1,jpj))  &
     &                            -(gu1(k+1,i  ,j  )*ru1(k+1,i  ,j  )  &
     &                             +gu1(k+1,im1,jpj)*ru1(k+1,im1,jpj))  &
     &                         +3.*(gu3(k  ,i  ,j  )*ru3(k  ,i  ,j  )  &
     &                             +gu3(k  ,im1,jpm)*ru3(k  ,im1,jpm))  &
     &                            -(gu3(k+1,i  ,j  )*ru3(k+1,i  ,j  )  &
     &                             +gu3(k+1,im1,jpm)*ru3(k+1,im1,jpm))  &
     &                         +3.*(gu2(k  ,i  ,j  )*ru2(k  ,i  ,j  )  &
     &                             +gu2(k  ,i  ,jv1)*ru2(k  ,i  ,jv1))  &
     &                            -(gu2(k+1,i  ,j  )*ru2(k+1,i  ,j)  &
     &                             +gu2(k+1,i  ,jv1)*ru2(k+1,i  ,jv1)))
                     w (k,i,j)=2.*rw(k,i,j)/(3.*rho(k,i,j)-rho(k+1,i,j))
                  do k=2,nz1
                     w (k,i,j) =rw(k,i,j)/(.5*(rho(k,i,j)+rho(k-1,i,j)))
                  end do
               end do
            end do
            if(iper.eq.0)  then
               do j=1,ny
                  do k=1,nz1
                     u1  (k,0,j) = ru1(k,0,j)/rho(k,1,j)
                     u3  (k,0,j) = ru3(k,0,j)/rho(k,1,j)
                  end do
               end do
            end if
            if(jper.eq.0)  then
               do i=1,nx
                  do k=1,nz1
                     u2  (k,i,0) = ru2(k,i,0)/rho(k,i,1)
                  end do
               end do
            end if
         else if ( rk_step .eq. 3) then
!
!********reset levels for full step
!
!        subtract out diabatic theta_v tendency from the last timestep
!
            do j=1,ny
               do i=1,nx
                  do k=1,nz1
                     rt1(k,i,j) = rt1(k,i,j)   &
     &                           - ns_rk*dts*rho(k,i,j)*t_d_tend(k,i,j)
                  end do
               end do
            end do

            do j=1,ny
               do i=1,nx
                  do k=1,nz1
                     ru1(k,i,j) = ru11(k,i,j)
                     ru3(k,i,j) = ru31(k,i,j)
                     ru2(k,i,j) = ru21(k,i,j)
                     rw (k,i,j) = rw1 (k,i,j)
                     rt (k,i,j) = rt1 (k,i,j)
                     rr (k,i,j) = rr1 (k,i,j)    
                     rho(k,i,j) = rb  (k,i,j)+rr(k,i,j)
                  end do
               end do
            end do
            do j=1,ny
               jp1 = min(j+1,ny)
               if(jper*j.eq.ny) jp1 = 2
               jv1 = j-1
               if(jper*j.eq.1)  jv1 = ny1
               jm1 = max(j-1,1)
               if(jper*j.eq.1)  jm1 = ny1
               do i=1,nx
                  if(mod(i,2).eq.0)  then
                     jpj = j
                     jpm = jm1
                  else
                     jpj = jp1
                     jpm = j
                  end if
                  ip1 = min(i+1,nx)
                  if(iper*i.eq.nx)  ip1 = 2
                  im1 = i-1
                  if(iper*i.eq.1 )  im1 = nx1
                  do k=1,nz1
                     u1(k,i,j)  = ru1(k,i,j)  &
     &                              /(.5*(rho(k,i,j)+rho(k,ip1,jpm)))
                     u3(k,i,j)  = ru3(k,i,j)  &
     &                              /(.5*(rho(k,i,j)+rho(k,ip1,jpj)))
                     u2(k,i,j)  = ru2(k,i,j)  &
     &                              /(.5*(rho(k,i,j)+rho(k,i  ,jp1)))
                     u11(k,i,j) = u1 (k,i,j)
                     u31(k,i,j) = u3 (k,i,j)
                     u21(k,i,j) = u2 (k,i,j)
                  end do
                  do k=1,nz1
                     t  (k,i,j) = (rtb(k,i,j)+rt(k,i,j))/rho(k,i,j)
                     t1 (k,i,j) = t(k,i,j)
                     p  (k,i,j) = (hh(i,j)*(rtb(k,i,j)+rt(k,i,j))  &
     &                                   /(100000./287./300.) )**rcv
                     do n = 1,nmoist
                       qx(k,i,j,n) = rqx(k,i,j,n)/rho(k,i,j)
                     enddo

                     if ( nscalar > 1 ) then
                     do n = 1,nscalar
                       sx(k,i,j,n) = rsx(k,i,j,n)/rho(k,i,j)
                     enddo
                     endif

!                      qv (k,i,j) = rqv(k,i,j)/rho(k,i,j)
!                      qc (k,i,j) = rqc(k,i,j)/rho(k,i,j)
!                      qr (k,i,j) = rqr(k,i,j)/rho(k,i,j)

                     qx1(k,i,j,1) = qx(k,i,j,1) ! qv
!                     qv1(k,i,j) = qv(k,i,j)
!                    qc1(k,i,j) = qc(k,i,j)
!                    qr1(k,i,j) = qr(k,i,j)

                     do n = 1,nmoist
                       rqx1(k,i,j,n) = rqx(k,i,j,n)
                     enddo

                     if ( nscalar > 1 ) then
                     do n = 1,nscalar
                       rsx1(k,i,j,n) = rsx(k,i,j,n)
                     enddo
                     endif
!                      rqv1(k,i,j) = rqv(k,i,j)
!                      rqc1(k,i,j) = rqc(k,i,j)
!                      rqr1(k,i,j) = rqr(k,i,j)
               
                  end do

                  k=1
                     rw(k,i,j) = - .25/hh(i,j)*  &
     &                         (3.*(gu1(k  ,i  ,j  )*ru1(k  ,i  ,j  )  &
     &                             +gu1(k  ,im1,jpj)*ru1(k  ,im1,jpj))  &
     &                            -(gu1(k+1,i  ,j  )*ru1(k+1,i  ,j  )  &
     &                             +gu1(k+1,im1,jpj)*ru1(k+1,im1,jpj))  &
     &                         +3.*(gu3(k  ,i  ,j  )*ru3(k  ,i  ,j  )  &
     &                             +gu3(k  ,im1,jpm)*ru3(k  ,im1,jpm))  &
     &                            -(gu3(k+1,i  ,j  )*ru3(k+1,i  ,j  )  &
     &                             +gu3(k+1,im1,jpm)*ru3(k+1,im1,jpm))  &
     &                         +3.*(gu2(k  ,i  ,j  )*ru2(k  ,i  ,j  )  &
     &                             +gu2(k  ,i  ,jv1)*ru2(k  ,i  ,jv1))  &
     &                            -(gu2(k+1,i  ,j  )*ru2(k+1,i  ,j)  &
     &                             +gu2(k+1,i  ,jv1)*ru2(k+1,i  ,jv1)))
                     w (k,i,j)=2.*rw(k,i,j)/(3.*rho(k,i,j)-rho(k+1,i,j))
                     w1(k,i,j)=   w (k,i,j)
                  do k=2,nz1
                     w (k,i,j) =rw(k,i,j)/(.5*(rho(k,i,j)+rho(k-1,i,j)))
                     w1(k,i,j) =w (k,i,j)
                  end do
               end do
            end do
            if(iper.eq.0)  then
               do j=1,ny
                  do k=1,nz1
                     ru1(k,0,j) = ru11(k,0,j)
                     u1 (k,0,j) = ru1 (k,0,j)/rho(k,1,j)
                     u11(k,0,j) = u1  (k,0,j)
                     ru3(k,0,j) = ru31(k,0,j)
                     u3 (k,0,j) = ru3 (k,0,j)/rho(k,1,j)
                     u31(k,0,j) = u3  (k,0,j)
                  end do
               end do
            end if
            if(jper.eq.0)  then
               do i=1,nx
                  do k=1,nz1
                     ru2(k,i,0) = ru21(k,i,0)
                     u2 (k,i,0) = ru2 (k,i,0)/rho(k,i,1)
                     u21(k,i,0) = u2  (k,i,0)
                  end do
               end do
            end if
         end if ! if rk_step = 2

      end do  ! rk step loop

!-----------------------------------------------------------------------
!
!  here is the microphysics (warm rain)
!
      do j=1,ny
         do i=1,nx
            do k=1,nz1
               t_d_tend(k,i,j) = t(k,i,j)
               t       (k,i,j) = t0*t(k,i,j)/(1.+1.61*qx(k,i,j,1))
            end do
         end do
      end do

      if ( mp_physics == 1 ) then
!      call kessler( t, qv, qc, qc1, qr, qr1, rho, p,
!     *              dt, dz, nz1, nx, ny                  )
!      call kessler( t, qv, qc, qc1, qr, qr1, rb, pb,
!     *              dt, dz, nz1, nx, ny                  )
!c    call kessler_joe( t, qv, qc, qc1, qr, qr1, rho, pb,
      call kessler_joe( t, qx(1,1,1,lv), qx(1,1,1,lc), qx(1,1,1,lr), rho, pb,  & ! 1=qv, 2=qc, 3=qr
     &              dt, dz, nz1, nx, ny                  )

      elseif ( mp_physics == 18 ) then

        ! ws(1:nz1,:,:) = w(1:nz1,:,:)
         iwmax = nxc
         jwmax = nyc
         kwmax = 1
         tmp = 1.0
         do j=1,ny
            do i=1,nx
               do k=1,nz1
                  ws(k,i,j) = w(k,i,j)
                  if(w(k,i,j).gt.tmp)  then
                     tmp = w(k,i,j)
                     iwmax = i
                     jwmax = j
                     kwmax = k
                  end if
               end do
            end do
         end do

      if ( debug ) then
      if ( nscalar > 1 ) then
        write(6,*) 'k,pres,w,qc,nc,qi,qh,t = '
      else
        write(6,*) 'k,pres,w,qc,qi,qh,t = '
      endif
      endif
      do j=1,ny
         do i=1,nx
            do k=1,nz1
              pres(k,i,j) = 1.e5*p(k,i,j)**(1./.2875)
              if ( debug .and. i == iwmax .and. j == jwmax ) then
               if ( nscalar > 1 ) then
                write(6,*) k,pres(k,i,j),ws(k,i,j),qx(k,i,j,lc)*1000., sx(k,i,j,lnc), &
                  qx(k,i,j,li)*1000.,qx(k,i,j,lh)*1000., pb(k,i,j)*t(k,i,j)-273.15
               else
                write(6,*) k,pres(k,i,j),ws(k,i,j),qx(k,i,j,lc)*1000., &
                  qx(k,i,j,li)*1000.,qx(k,i,j,lh)*1000., pb(k,i,j)*t(k,i,j)-273.15
               endif
              endif
            enddo
         enddo
      enddo

!         pres(1:nz1,1:nx,1:ny) = 1.e5*p(1:nz1,1:nx,1:ny)**(1./.2875)
         
         CALL nssl_2mom_driver(                          &
                     ITIMESTEP=nit,                      &
                     TH=t,                              &
                     QV=qx(1,1,1,lv),                         &
                     QC=qx(1,1,1,lc),                         &
                     QR=qx(1,1,1,lr),                         &
                     QI=qx(1,1,1,li),                         &
                     QS=qx(1,1,1,ls),                         &
                     QH=qx(1,1,1,lh),                         &
                     QHL=qx(1,1,1,lhl),                        &
 !                    CCW=qnc_curr,                       &
                     CCW=sx(1,1,1,lnc),                    &
                     CRW=sx(1,1,1,lnr),                       &
                     CCI=sx(1,1,1,lni),                       &
                     CSW=sx(1,1,1,lns),                       &
                     CHW=sx(1,1,1,lnh),                       &
                     CHL=sx(1,1,1,lnhl),                       &
                     VHW=sx(1,1,1,lvh), f_vhw=(lvh > 1),      &
                     VHL=sx(1,1,1,lvhl), f_vhl=(lvhl > 1),      &
!                      ZRW=qzr_curr,  f_zrw = f_qzr,       &
!                      ZHW=qzg_curr,  f_zhw = f_qzg,       &
!                      ZHL=qzh_curr,  f_zhl = f_qzh,       &
                     cn=sx(1,1,1,lccn),  f_cn=(lccn > 1),    &
                     PII=pb,                               &
                     P=pres,                                &
                     W=ws,                               &
                     DZ=dz3d,                            &
                     DTP=dt,                             &
                     DN=rho,                             &
                      RAINNC   = RAINNC,                  &
                      RAINNCV  = RAINNCV,                 &
!                      SNOWNC   = SNOWNC,                  &
!                      SNOWNCV  = SNOWNCV,                 &
!                      HAILNC   = HAILNC,                  &
!                      HAILNCV  = HAILNCV,                 &
!                      GRPLNC   = GRAUPELNC,               &
!                      GRPLNCV  = GRAUPELNCV,              &
!                      SR=SR,                              &
                     dbz      = dbz      ,               &
!                     ssat3d   = ssat,  f_ssat=f_ssat,    &
!                     ssati    = ssati, f_ssati=f_ssati,  &
!#if ( WRF_CHEM == 1 )
!                    WETSCAV_ON = config_flags%wetscav_onoff == 1, &
!                    EVAPPROD=evapprod,RAINPROD=rainprod, &
!#endif
                     nssl_progn=.false.,              &
                     diagflag = .true.,                &
                     ke_diag = nz1,                &
!                      cu_used=cu_used,                    &
!                      qrcuten=qrcuten,                    &  ! hm
!                      qscuten=qscuten,                    &  ! hm
!                      qicuten=qicuten,                    &  ! hm
!                      qccuten=qccuten,                    &  ! hm
!                      re_cloud=re_cloud,                  &
!                      re_ice=re_ice,                      &
!                      re_snow=re_snow,                    &
!                      has_reqc=has_reqc,                  & ! ala G. Thompson
!                      has_reqi=has_reqi,                  & ! ala G. Thompson
!                      has_reqs=has_reqs,                  & ! ala G. Thompson
!                      hail_maxk1=hail_maxk1,              &
!                      hail_max2d=hail_max2d,              &
                     nwp_diagnostics=0, &
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde, &
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme, &
                  ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte  &
                                                                    )
      
      endif

      do j=1,ny
         do i=1,nx
            do k=1,nz1

               t   (k,i,j) = t(k,i,j)*(1.+1.61*qx(k,i,j,1))/t0 
               t_d_tend(k,i,j) = (t(k,i,j) - t_d_tend(k,i,j))/(ns*dts)
!c             t_d_tend(k,i,j) = 0.
               t1  (k,i,j) = t(k,i,j)
               rt  (k,i,j) = t(k,i,j)*rho(k,i,j) - rtb(k,i,j)
               rt1 (k,i,j) = rt(k,i,j)
        
               do n = 1,nmoist
                 rqx (k,i,j,n) = qx(k,i,j,n)*rho(k,i,j)
                 rqx1(k,i,j,n) = rqx(k,i,j,n)
                 qx1 (k,i,j,n) = qx(k,i,j,n)
               enddo

               if ( nscalar > 1 ) then
               do n = 1,nscalar
                 rsx (k,i,j,n) = sx(k,i,j,n)*rho(k,i,j)
                 rsx1(k,i,j,n) = rsx(k,i,j,n)
                 sx1 (k,i,j,n) = sx(k,i,j,n)
               enddo
               endif

!                rqv (k,i,j) = qv(k,i,j)*rho(k,i,j)
!                rqv1(k,i,j) = rqv(k,i,j)
!                qv1 (k,i,j) = qv(k,i,j)
! 
!                rqc (k,i,j) = qc(k,i,j)*rho(k,i,j)
!                rqc1(k,i,j) = rqc(k,i,j)
!                qc1 (k,i,j) = qc(k,i,j)
! 
!                rqr (k,i,j) = qr(k,i,j)*rho(k,i,j)
!                rqr1(k,i,j) = rqr(k,i,j)
!                qr1 (k,i,j) = qr(k,i,j)

               p   (k,i,j) = (hh(i,j)*(rtb(k,i,j)+rt(k,i,j))  &
     &                               /(100000./287./300.))**rcv
            end do
         end do
      end do

!      wmax = 0.
!      do j=1,ny
!         do i=1,nx
!            do k=2,nz1
!               wmax = amax1(wmax,w(k,i,j))
!            end do
!         end do
!      end do

!      vdiffm = 0.
!      do j=1,ny
!         do i=1,nx
!            do k=1,nz1
!              if(mod(i,2).eq.0.)  then
!                 nyj=ny+1-j
!                 nyj=nyc+1-j
!                 if(nyj.lt.1)  nyj=nyj+ny1
!                  vdiff=abs(u2(k,i,j)+u2(k,i,nyj))
!              else
!                 nyj=ny-j
!                 nyj=nyc-j
!                 if(nyj.lt.1)  nyj=nyj+ny1
!                  vdiff=abs(u2(k,i,j)+u2(k,i,nyj))
!              end if
!              if(vdiff.gt.vdiffm)  then
!                  vdiffm = vdiff
!                 ivm=i
!                 jvm=j
!                 kvm=k
!              end if
!            end do
!         end do
!         write(6,*) j,yu2(91,j),u2(1,91,j),yu2(92,j),u2(1,92,j)
!      end do
!     write(6,*) vdiffm,u2(kvm,ivm,jvm),ivm,jvm,kvm

      end if !  take step only after plotting first
!
!
!**** processing for plotting
!
!--------------
!
      include "plotting.inc.f90"
!
!--------------
!
      end do  ! for timestep loop

      call clsgks()

      stop
      end
     
