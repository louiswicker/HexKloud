c                                                                     72
      pi    = 4.*atan(1.)
cc      imass = 1
      imass = 0
c      iper  = 0
      iper  = 1
c      jper  = 0
      jper  = 1

c      ur   = 0.
cc      ur   = -6.
      ur   = -15.
      angle= 0.
c      angle= 45.
      um    =ur*cos(angle/180.*pi)
      vm    =ur*sin(angle/180.*pi)
      u1m  = .5*sqrt(3.)*um-.5*vm
      u3m  = .5*sqrt(3.)*um+.5*vm
      u2m  = vm
      write(6,*) 'U1M = ',u1m,'  U3M = ',u3m,' U2M = ',u2m

c      xl    = 20000.
c      xl    = 28000.
      xl    = 84000.
cc      xl    = 56000.
      side  = 2.*xl/(3.*float(nx1))
      d     = sqrt(3.)*side
      yl    = d*float(ny1)

c     initialization for 2-d y-z simulation
c      yl    = 50000.
c      d     = yl/float(ny1)
c      side  = d/sqrt(3.)
c      xl    = 1.5*side*float(nx1)

      write(6,*) 'XL = ',xl,'   YL = ',yl,'   D = ',d

      xn2   = 0.0001
      xn2m  = 0.0001
      xn2l  = 0.0001
      zinv  = 10000.
      xn    = sqrt(xn2)
c      f     = .0001
      f     = 0.

      hm     = 0.
      ampl  = 1.
      xa    = 10000.
      ya    = 10000.

c      dt    = 6.0
      dt    = 12.0
c      dt    = 24.0
	  
c      ns0   = 6
      ns0   = 8
c      ns0   = 4
	  
      delt= 3.
c      delt= 2.


      vnu   = 500.
      xnu   = vnu*dt/d**2
	  
cc      xnu   = .015
cc      xnu   = .008
cc      xnu   = .012
cc      xnu   = .006 ! original value
cc      xnu   = .003
c      xnu   = .0

      smdiv = .1
c      smdiv = .0
      epssm = .1

      tip   = 300. ! plotting interval

      ip    = nint(tip/dt)
	  
c	  ip = 1
	  
      tstp  = 12.*ip*dt ! total time is 12 plotting intervals
      nz2   = nz1-1
      t0    = 300.
      r     = 287.
      cp    = 1003.
      rcv   = r/(cp-r)
      p0    = 100000.
      cti   = 1./(cp*t0)
      g     = 9.81
      c2    = cp*rcv*t0
      cb    = 25.
	  delt = delt/t0

      zt    = 20000.
      zd    = 10000.
c      xnut = .025
c      xnut = .05
      xnut = .0

c

c      xpll   = 0.
c      xplr   = xl
cc      xpll   = 0.25*xl
cc      xplr   = 0.75*xl
      xpll   =    xl/6.
      xplr   = 5.*xl/6.
c      ypll   = 0.
c      yplr   = yl
cc      ypll   = 0.25*yl
cc      yplr   = 0.75*yl
      ypll   =    yl/6.
      yplr   = 5.*yl/6.
      zplb   =  0.
c      zplt   =  zd
      zplt   =  zt
	  wmplt  = 50.

      xc     = .5*xl
      yc     = .5*yl
      dx     = xl/float(nx1)
      dy     = yl/float(ny1)
      dz     = zt/float(nz1)
      nxc    = nx1/2+1
      if(mod(nxc,2).eq.0) nxc=nxc-1
      nyc    = ny1/2+1
      dxp    = 100.*nint(sqrt(xl*yl/float(nx1*ny1))/100.)
	  dyp    = dxp
	  dzp    = dz

      dts    = dt/float(ns0)
      dtseps = .25*dts*(1.+epssm)
      ns     = ns0
      rd     = 1./d
      rdz    = 1./dz
      area   = 1.5*d*side
      dtsa   = dts*side/area
      dtsd   = dts/d
      dtsf   = dts*sqrt(3.)/6.*f
      dtsg   = dts*.5/d*g
      c1f    = 2./3.
      c2f    = 1./3.

      resm   = (1.-epssm)/(1.+epssm)
      smdivx = smdiv*d **2/dts
      smdivz = smdiv*dz**2/dts

      xnus0   = xnu/dt
      xnusz0  = xnus0*rdz**2*d**2

cc      xnusz0 = .25*xnusz0

      xnus0   = 2./3.*xnus0      

c*****adjustment to make mixing the same as xyz code
c      deltax = sqrt(xl*yl/float(nx1*ny1))
c	  write(6,*) 'deltax = ',deltax,' dz = ',dz
c	  xnusz0 = xnus0*deltax**2/dz**2
c	  xnus0  = xnus0*2./3.*deltax**2/d**2
      write(6,*) 'XNU = ',xnu/dt,' XNUS = ',xnus0,' XNUSZ = ',xnusz0
	  
      iplt   = 1
      ipi    = nint(xpll/dx)+1
      ipf    = nint(xplr/dx)+1
      jpi    = nint(ypll/dy)+1
      jpf    = nint(yplr/dy)+1
      nxpl   = ipf-ipi+1
      nypl   = jpf-jpi+1
      npl    = nxpl
      nzpl   = nint(nz1*zplt/zt)
      dtl    = dt
      zw(1)  = 0.
      ax(1)  = 0.
	  write(6,*) 'NX = ',nx,' NXC = ',nxc,' NY = ',ny,' NYC = ',nyc
      write(6,*) 'xa =',xa,' xl =',xl,' yl =',yl,' dt   =',dt
      write(6,*) 'zt =',zt,' zd =',zd,' xnut =',xnut
c
      do j=1,ny
         do i=1,nx
            xh(i,j) = (i-1)*1.5*side
            x (i)   = xh(i,1)
            if(mod(i,2).eq.0)  then
               yh(i,j) = (j-1.5)*d
            else
               yh(i,j) = (j- 1.)*d
            end if
            xu3(i,j) = xh(i,j)+.75*side
            yu3(i,j) = yh(i,j)+d/4.
            xu1(i,j) = xh(i,j)+.75*side
            yu1(i,j) = yh(i,j)-d/4.
            xu2(i,j) = xh(i,j)
            yu2(i,j) = yh(i,j)+d/2.
         end do
         y (j)   = yh(1,j)
      end do
c
      do k=1,nz1
         zu(k  )  = dz*(k-.5)
         zw(k+1)  = dz*k
         ax(k+1)  = 0.
         ds(k)    = 0.
         if(zu(k).le.zinv)  then
            tz(k)=1.+xn2l/g*zu(k)
         else
            tz(k)=1.+xn2l/g*zinv+xn2/g*(zu(k)-zinv)
         end if
         if(zu(k).gt.zd+.1)  then
            ds(k) = ds(k)+xnut*sin(.5*pi*(zu(k)-zd)/(zt-zd))**2
         end if
         ds(k) = ds(k)/float(2*ns0)
      end do
c
c----------------
c
      call init_sound( tz,rel_hum,u1z,u2z,u3z,zu,nz1 )
c
c----------------
c
c  nondimensionalize the sounding
c
      do k=1,nz1
        tz(k) = (tz(k))/t0
        tzv(k) = tz(k)
cc        rel_hum(k) = 0.
        qvzv(k) = 0.
      enddo
c
      do j=1,ny
         do i=1,nx
            hs(i,j) = hm/(1.+((xh(i,j)-xc)/xa)**2+((yh(i,j)-yc)/ya)**2)
            if(mod(i,2).eq.0)  then
               hxpl(i) = ampl*.5*(hs(i,nyc)+hs(i,nyc+1))
            ELSE
               hxpl(i) = ampl*hs(i,nyc)
            END IF
            hh(i,j) = zt/(zt-hs(i,j))
         end do
      end do
c
      if(iper.eq.0)  then
         do j=1,ny
            jp1 = min(j+1,ny)
            if(jper*j.eq.ny) jp1=2
            jm1 = max(j-1,1)
            if(jper*j.eq.1) jm1=ny1
            do k=1,nz1
               gu1(k,0,j) = (zu(k)-zt)/(zt-hs(1,jm1))
     &                                *rd*(hs(2,jm1)-hs(1,jm1))
               gu3(k,0,j) = (zu(k)-zt)/(zt-hs(1,j  ))
     &                                *rd*(hs(2,jp1)-hs(1,jp1-1))
               u1(k,0,j)  = u1m + u1z(k)
               u3(k,0,j)  = u3m + u3z(k)
            end do
         end do
      end if
      if(jper.eq.0)  then
         do i=1,nx
            do k=1,nz1
               gu2(k,i,0) = (zu(k)-zt)/(zt-hs(i,1))*rd*(hs(i,2)-hs(i,1))
               u2(k,i,0)  = u2m + u2z(k)
            end do
         end do
      end if

c--------------------------------------------------------------
      do ittrm = 1,10
c--------------------------------------------------------------


      do j=1,ny
         jp1 = min(j+1,ny)
         if(jper*j.eq.ny) jp1=2
         jm1 = max(j-1,1)
         if(jper*j.eq.1  )  jm1=ny1
         do i=1,nx
	        j1  = j
	        jn  = j
            if(mod(i,2).eq.0)  then
               jpj = j
               jpm = jm1
	       if(jper.eq.0.and.j.eq.1 )  j1 = 2
            else
               jpj=jp1
               jpm=j
	       if(jper.eq.0.and.j.eq.ny)  jn = ny1
            end if
            ip1 = min(i+1,nx)
            if(iper*i.eq.nx) ip1=2
            dhh1(i,j) = (hh(ip1,jpm)-hh(i,j1))/(hh(ip1,jpm)+hh(i,j1))
            dhh3(i,j) = (hh(ip1,jpj)-hh(i,jn))/(hh(ip1,jpj)+hh(i,jn))
            dhh2(i,j) = (hh(i  ,jp1)-hh(i,j ))/(hh(i  ,jp1)+hh(i,j ))
            do k=1,nz1
               ztemp    = zu(k)+(1.-zu(k)/zt)*hs(i,j)
               gu1(k,i,j) = (zu(k)-zt)/(zt-.5*(hs(ip1,jpm)+hs(i,j)))
     &                                    *rd*(hs(ip1,jpm)-hs(i,j1))
               gu3(k,i,j) = (zu(k)-zt)/(zt-.5*(hs(ip1,jpj)+hs(i,j)))
     &                                    *rd*(hs(ip1,jpj)-hs(i,jn))
               gu2(k,i,j) = (zu(k)-zt)/(zt-.5*(hs(i  ,jp1)+hs(i,j)))
     &                                    *rd*(hs(i  ,jp1)-hs(i,jp1-1))
               u1(k,i,j)  = u1m + u1z(k)
               u3(k,i,j)  = u3m + u3z(k)
               u2(k,i,j)  = u2m + u2z(k)
cc               tb (k,i) = tzv(k)
c               if(ztemp.le.zinv)  then
c                  t(k,i,j)=1.+xn2l/g*ztemp
c               else
c                  t(k,i,j)=1.+xn2l/g*zinv+xn2/g*(ztemp-zinv)
c               end if
               t (k,i,j) = tz(k)
               ti(k,i,j) = t(k,i,j)
               tb(k,i,j) = t(k,i,j)
            end do
         end do
      end do
c
c  assume zero moisture to start, itterate hydrostatic equation
c
c      do ittrm = 1,10
c--------------------------------------------------------------

      fac1 = 0.61
      pitop = 1.-.5*dz*g/(cp*t (1,1,1)*(1.+fac1*qvzv(1))*t0*hh(1,1))
      do k=2,nz1
         tk = t(k,1,1)*(1.+fac1*qvzv(k))
         tkm1 = t(k-1,1,1)*(1.+fac1*qvzv(k-1))
         pitop = pitop-dz*g/(cp*.5*(tk+tkm1)*t0*hh(1,1))
      end do
      tk = t(nz1,1,1)*(1.+fac1*qvzv(nz1))
      pitop = pitop-.5*dz*g/(cp*t (nz1,1,1)*t0*hh(1,1))
      do j=1,ny
         do i=1,nx
            tk = tb(nz1,i,j)*(1.+fac1*qvzv(nz1))
            pb(nz1,i,j) = pitop+.5*dz*g/(cp*tk*t0*hh(i,j))
            tk = t(nz1,i,j)*(1.+fac1*qvzv(nz1))
            p (nz1,i,j) = pitop+.5*dz*g/(cp*tk*t0*hh(i,j))
            do k=nz2,1,-1
               tk = tb(k,i,j)*(1.+fac1*qvzv(k))
               tkp1 = tb(k+1,i,j)*(1.+fac1*qvzv(k+1))
               pb(k,i,j)  = pb(k+1,i,j)
     &                      + dz*g/(cp*.5*(tk+tkp1)*t0*hh(i,j)) 
               tk = t(k,i,j)*(1.+fac1*qvzv(k))
               tkp1 = t(k+1,i,j)*(1.+fac1*qvzv(k+1))
               p (k,i,j)  = p (k+1,i,j)
     &                      + dz*g/(cp*.5*(tk+tkp1)*t0*hh(i,j))
            end do
            do k=1,nz1
               rb(k,i,j)  = (100000./287./300.)*pb(k,i,j)**(1./rcv)/
     &                      ((tb(k,i,j)*(1.+1.61*qvzv(k)))*hh(i,j))
               rr(k,i,j)  = (100000./287./300.)*p (k,i,j)**(1./rcv)/
     &                      ((t (k,i,j)*(1.+1.61*qvzv(k)))*hh(i,j))
     &                      -rb(k,i,j)
               rtb(k,i,j) = rb(k,i,j)*tb(k,i,j)*(1.+1.61*qvzv(k))
               rti(k,i,j) = t(k,i,j)*rr(k,i,j)+rb(k,i,j)*
     &                      (t(k,i,j)-tb(k,i,j))*(1.+1.61*qvzv(k))
           end do
         end do
      end do
c
      rttop = 0.
      do itr=1,10
         tmax=0.
         do j=1,ny
            do i=1,nx
               rt(nz1,i,j) = rttop-.5*dz*g/(c2*p(nz1,i,j)*hh(i,j))*
     &                 (rb(nz1,i,j)*(p(nz1,i,j)-pb(nz1,i,j))/pb(nz1,i,j)
     &                 -rr(nz1,i,j))
               do k=nz2,1,-1
                  rt(k,i,j) = rt(k+1,i,j)-dz*g/(c2*(p(k,i,j)+p(k+1,i,j))       
     &                        *hh(i,j))*
     &                          (rb(k  ,i,j)*(p(k  ,i,j)-pb(k  ,i,j))
     &                          /pb(k  ,i,j)-rr(k  ,i,j)
     &                          +rb(k+1,i,j)*(p(k+1,i,j)-pb(k+1,i,j))
     &                          /pb(k+1,i,j)-rr(k+1,i,j))
               end do
               do k=1,nz1
                  tdiff=abs(rt(k,i,j)-rti(k,i,j))
                  if(tdiff.gt.tmax)  tmax = tdiff
                  p  (k,i,j) = (hh(i,j)*(rtb(k,i,j)+rt(k,i,j))
     &                                 /(100000./287./300.))**rcv
                  rr (k,i,j) = rt(k,i,j)-rb(k,i,j)*
     &                         (t(k,i,j)-tb(k,i,j))*(1.+1.61*qvzv(k))
     &                        /(t(k,i,j)*(1.+1.61*qvzv(k)))
                  rri(k,i,j) = rr(k,i,j)
                  rti(k,i,j) = rt(k,i,j)
                  pii(k,i,j) = p (k,i,j)
               end do
            end do
         end do
       end do
c
c     potential temperature perturbation
c
      radx  = 10000.
      rady  = 10000.
      radz  = 1500.
      zcent = 1500.
cc      radz = 2000.
cc      zcent = 3000.
      if(delt.ne.0.)  then
         do j=1,ny
            do i=1,nx
            xht = (i-1)*1.5*side
            if(mod(i,2).eq.0)  then
               yht = (j-1.5)*d
            else
               yht = (j- 1.)*d
            end if
               do k=1,nz1
                  ztemp      = zu(k)+(1.-zu(k)/zt)*hs(i,j)
c                  rad=sqrt(((xh(i,j)-xc)/radx)**2
c     &                    +((ztemp-zcent)/radz)**2)
c                  rad=sqrt(((yh(i,j)-yc)/rady)**2
c     &                    +((ztemp-zcent)/radz)**2)
c                  rad=sqrt(((xh(i,j)-xc)/radx)**2+((yh(i,j)-yc)/rady)**2
c     &                    +((ztemp-zcent)/radz)**2)
                  rad=sqrt(((xht-xc)/radx)**2+((yht-yc)/rady)**2
     &                    +((ztemp-zcent)/radz)**2)
                  if(rad.lt.1)  then
                     t(k,i,j)=t(k,i,j)+delt*cos(.5*pi*rad)**2
c equiv line
c from joe's code
c                  TTEMP = .5*DELT*(COS(PI*RAD)+1.)+TTEMP
c
                  end if
                  rt1(k,i,j)=rt(k,i,j)
               end do
            end do
         end do

c         do j=1,ny
c            do i=1,nx
cc			   ii=i+nxc-1
cc			   if(i.gt.nxc)  ii=nx+1-i
c			   ii =i
c	   	       jj=j+nyc-1
c			   if(j.gt.nyc)  jj=ny+2-mod(i,2)-j
c              do k=1,nz1
c                  t  (k,i,j)=t (k,ii,jj)
c               end do
c            end do
c         end do

         do itr=1,10
            tmax=0.
            do j=1,ny
               do i=1,nx
                  rt(nz1,i,j) = rttop-.5*dz*g/(c2*p(nz1,i,j)*hh(i,j))*
     &                         (rb(nz1,i,j)*(p(nz1,i,j)-pb(nz1,i,j))/
     &                          pb(nz1,i,j)-rr(nz1,i,j))
                  do k=nz2,1,-1
                     rt(k,i,j) = rt(k+1,i,j)-.5*dz*g/(c2*p(k,i,j)
     &                           *hh(i,j))*
     &                             (rb(k  ,i,j)*(p(k  ,i,j)-pb(k  ,i,j))
     &                             /pb(k  ,i,j)-rr(k  ,i,j)
     &                             +rb(k+1,i,j)*(p(k+1,i,j)-pb(k+1,i,j))
     &                             /pb(k+1,i,j)-rr(k+1,i,j))
                  end do
                  do k=1,nz1
                     tdiff=abs(rt(k,i,j)-rt1(k,i,j))
                     if(tdiff.gt.tmax)  tmax=tdiff
                     p  (k,i,j) = (hh(i,j)*(rtb(k,i,j)+rt(k,i,j))
     &                            /(100000./287./300.))**rcv
                     rr (k,i,j) = (rt(k,i,j)-rb(k,i,j)
     &                            *(t(k,i,j)-tb(k,i,j)))/t(k,i,j)
                     rt1(k,i,j) =  rt(k,i,j)
                  end do
               end do
            end do
         end do
      end if
c
      ritot=0.
      rrtot=0.
      do j=1,ny-jper
         do i=1,nx-iper
            do k=1,nz1
               ritot = ritot+rb(k,i,j)+rr(k,i,j)
               rrtot = rrtot+rr(k,i,j)
            end do 
         end do
      end do
c
      do j=1,ny
         do i=1,nx
            do k=1,nz1
               rho(k,i,j) = rb(k,i,j)+rr(k,i,j)
            end do
         end do
      end do

c
c  change humidity to mixing ratio
c
      do k=1,nz1
          temp      = (pii(k,1,1))*t0*tz(k)
          pressure = 1.000e+05 * ((pii(k,1,1))**(1003./287.))
          qvs       = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pressure
          qvzv(k) = amin1(0.014,rel_hum(k)*qvs)
      end do
c
c  recompute the virtual temperature
c
      do k=1,nz1
        tzv(k) = tz(k)*(1.+1.61*qvzv(k))
      end do
      do j=1,ny
         do i=1,nx
            do k=1,nz1
               ti(k,i,j) = tzv(k)
            end do
         end do
      end do
c
c  end of loop for moist iteration
c
      end do
c
c  write out the sounding
c
      write(6,*) ' sounding for the simulation '
      do k=1,nz1
         write(6,*) zu(k),t0*tz(k),1000.*qvzv(k),
     *              rb(k,1,1)*(1.+qvzv(k)),u1z(k),u2z(k),u3z(k)
      enddo
c
c reset the temperature to be the "klemp virtual temperature"
c
      do k=1,nz1
         rqvb(k) = qvzv(k)*rb(k,1,1)
      end do
      do j=1,ny
         do i=1,nx
            do k=1,nz1
               t (k,i,j) = t (k,i,j)*(1.+1.61*qvzv(k))
               tb(k,i,j) = tb(k,i,j)*(1.+1.61*qvzv(k))
            end do
         end do
      end do
c
c  note, qvzv is the
c  mixing ratio, both for the initial state.
c  rqvb(k) is rho_d * qv_init
c
c----------------------------------------------------------
      do j=1,ny
         jm1 = max(j-1,1)
         if(jper*j.eq.1 )  jm1 = ny1
         jp1 = min(j+1,ny)
         if(jper*j.eq.ny)  jp1 = 2
         do i=1,nx
            ip1 = min(i+1,nx)
            if(iper*i.eq.nx)  ip1 = 2
            if(mod(i,2).eq.0)  then
               jpj=j
               jpm=jm1
            else
               jpj=jp1
               jpm=j
            end if
            do k=1,nz1
               km1      = max(k-1,1)
               u11 (k,i,j) = u1(k,i,j)
               ru1 (k,i,j) = .5*(rho(k,i,j)+rho(k,ip1,jpm))*u1(k,i,j)    
               ru11(k,i,j) = ru1(k,i,j)
               ru1i(k,i,j) = ru1(k,i,j)
               u31 (k,i,j) = u3(k,i,j)
               ru3 (k,i,j) = .5*(rho(k,i,j)+rho(k,ip1,jpj))*u3(k,i,j)    
               ru31(k,i,j) = ru3(k,i,j)
               ru3i(k,i,j) = ru3(k,i,j)
               u21 (k,i,j) = u2(k,i,j)
               ru2 (k,i,j) = .5*(rho(k,i,j)+rho(k,i  ,jp1))*u2(k,i,j)    
               ru21(k,i,j) = ru2(k,i,j)
               ru2i(k,i,j) = ru2(k,i,j)
               w   (k,i,j) = 0.
               w1  (k,i,j) = 0.
               ww  (k,i,j) = 0.
               rw  (k,i,j) = .5*(rho(k,i,j)+rho(km1,i,j))*w(k,i,j)
               rw1 (k,i,j) = rw (k,i,j)
               t1  (k,i,j) = t (k,i,j)
               rt1 (k,i,j) = rt (k,i,j)
               rr1 (k,i,j) = rr(k,i,j)
               div (k,i,j) = 0.
               qv  (k,i,j) = qvzv(k)
               qv1 (k,i,j) = qvzv(k)
               rqv (k,i,j) = qv(k,i,j)*rho(k,i,j)
               rqv1(k,i,j) = qv(k,i,j)*rho(k,i,j)
               qc  (k,i,j) = 0.
               qr  (k,i,j) = 0.
               qc1 (k,i,j) = 0.
               qr1 (k,i,j) = 0.
               rqc (k,i,j) = 0.
               rqr (k,i,j) = 0.
               rqc1(k,i,j) = 0.
               rqr1(k,i,j) = 0.   
            end do
            w  (nz,i,j) = 0.
            w1 (nz,i,j) = 0.
            ww (nz,i,j) = 0.
            rw (nz,i,j) = 0.
            rw1(nz,i,j) = 0.
         end do 
      end do 
	  
      if(iper.eq.0)  then
         do j=1,ny
            do k=1,nz1
               u11 (k,0,j)  = u1 (k,0,j)
               ru1 (k,0,j)  = rho(k,1,j)*u1(k,0,j)
               ru11(k,0,j)  = ru1(k,0,j)
               ru1i(k,0,j)  = ru1(k,0,j)
               u31 (k,0,j)  = u3 (k,0,j)
               ru3 (k,0,j)  = rho(k,1,j)*u3(k,0,j)
               ru31(k,0,j)  = ru3(k,0,j)
               ru3i(k,0,j)  = ru3(k,0,j)
            end do
         end do
      end if
      if(jper.eq.0)  then
         do i=1,nx
            do k=1,nz1
               u21 (k,i,0)  = u2  (k,i,0)
               ru2 (k,i,0)  = rho(k,i,1)*u2(k,i,0)
               ru21(k,i,0)  = ru2 (k,i,0)
               ru2i(k,i,0)  = ru2 (k,i,0)
            end do
         end do
      end if
      do k=1,nz
         wdtz(k) = 0.
      end do
      cofrz = 2.*dtseps*rdz
      do j=1,ny
         do i=1,nx
            cofwr(i,j) = dtseps*g*hh(i,j)
            a    (1,i,j) = 0.
            b    (1,i,j) = 1.
            c    (1,i,j) = 0.
            gamma(1,i,j) = 0.
         end do
      end do
      if(iplt.eq.1)   then
        IF ( .true. ) THEN
C
C Open GKS.
C
        CALL GOPKS (IERF,0)
        CALL GOPWK (IWID,LUNI,IWTY)
        CALL GACWK (IWID)

         ELSE

         call opngks

         ENDIF
         call frame
      end if
      kkk  = 0
      nit  = 0
      npr  = 0
      time = 0.
      total_steps = nint(tstp/dt)
      ipp  = 1
      jpp  = 1
      write(6,*) ' sounding at 1 '
      write(6,*) ' k, z, theta, rho*thetam, pii, rho, qv0 '
      do k=nz1,1,-1
         thetak = t0*tb(k,ipp,jpp)/(1.+1.61*qvzv(k))
         write(6,666) 
     *      k,zu(k),thetak,t0*(rt(k,ipp,jpp)+rtb(k,ipp,jpp)),
     *      pb(k,ipp,jpp),rb(k,ipp,jpp),qvzv(k)
      enddo
cc      stop
 666  format(1x,i3,1x,f6.0,5(1x,e13.6))

         do j=1,ny
            do i=1,nx
               do k=1,nz1
                  t_d_tend(k,i,j) = 0.
               end do
            end do
         end do
