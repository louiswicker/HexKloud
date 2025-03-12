!                                                                     72
         iwmax = nxc
         jwmax = nyc
         kwmax = 1
         wmax (nit+1) = 0.
         waxis(nit+1) = nit 
         do j=1,ny
            do i=1,nx
               do k=1,nz1
                  if(w(k,i,j).gt.wmax(nit+1))  then
                     wmax(nit+1) = w(k,i,j)
                     iwmax = i
                     jwmax = j
                     kwmax = k
                  end if
               end do
            end do
         end do
         
!        if((kkk.ge.ip).or. (nit.eq.0)) then      
        IF ( outputflag ) THEN

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

            print 901 ,time
  901       format(1h ,'time =',f10.1)
            kkk=0
! 111        continue
            rtot=0.
            do j=1,ny-jper
               do i=1,nx-iper
                  do k=1,nz1
                     rtot = rtot+rr(k,i,j)
                  end do 
               end do
            end do

            rtot  = 100.*(rtot-rrtot)/ritot
            write(6,*) '       percent mass change1 =',rtot

            if(iplt.eq.1)  then       

            pzl = .05
            zptop = 0.05 + 0.9*0.8
            zptop = .95
        if(xplr-xpll.ge.yplr-ypll)  then
           pxl = .05
           pxr = .95
           pyl = .5 - .45*(yplr-ypll)/(xplr-xpll)
           pyr = .5 + .45*(yplr-ypll)/(xplr-xpll)
        else
           pyl = .05
           pyr = .95
           pxl = .5 - .45*(xplr-xpll)/(yplr-ypll)
           pxr = .5 + .45*(xplr-xpll)/(yplr-ypll)
        end if
!
!************  Maximum vertical velocity
!
        call wplot(wmax,waxis,wmplt,ip,nit+1,total_steps)
        write(6,*) 'i,j,k,wmax: ', IWMAX,JWMAX,KWMAX,WMAX(nit+1)
        
         if ( doplot ) then

!            go to 122

!
!************  x-z cross sections
!
!               j = nyc
               j = jwmax

               write(slice(2),300)  j
  300          format(i3)
               slice(1) = ' J='
               
               jp1 =min(j+1,ny)
               if(jper*j.eq.ny)  jp1 = 2
               jm1 = j-1
               if(jper*j.eq.1 )  jm1 = ny1
!
!************  x-z perturbation theta cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = .5*  &
     &                      (t(k,ii,j  )/(1.+1.61*qx(k,ii,j  ,lv))  &
     &                      +t(k,ii,jp1)/(1.+1.61*qx(k,ii,jp1,lv)))
                     else
                        pltx(i,k) = t(k,ii,j)/(1.+1.61*qx(k,ii,j,lv))
                     end if
                     pltx(i,k) = t0*(pltx(i,k)-tz(k)) 
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'t ',plane,'p',  &
     &           time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
               call setusv('LW',1000)
!               call curve(x,hxpl,nx)
               call frame
!
!************  x-z perturbation u cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  im1=ii-1
                  if(iper*ii.eq.1) im1 = nx1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = ampl*.5/sqrt(3.)  &
     &                                    *(u1(k,ii,jp1)+u1(k,im1,j)  &
     &                                     +u3(k,ii,j  )+u3(k,im1,j)  &
     &                                  -2.*(u1z(k)+u1m+u3z(k)+u3m))
                     else
                        pltx(i,k) = ampl*.5/sqrt(3.)  &
     &                                    *(u1(k,ii,j)+u1(k,im1,jp1)  &
     &                                     +u3(k,ii,j)+u3(k,im1,j )  &
     &                                  -2.*(u1z(k)+u1m+u3z(k)+u3m))
                     end if
                  end do
               end do
!               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'u ',plane','p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z perturbation v cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = ampl*(u2(k,ii,j)-u2z(k)-u2m) 
                     else
                        pltx(i,k) = ampl*(.5*(u2(k,ii,j)+u2(k,ii,jm1))  &
     &                                      -u2z(k)-u2m)   
                     end if
                  end do
               end do
!              call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'v ',plane,'p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z perturbation w cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k)=0.5*ampl*(w(k,ii,j)+w(k,ii,jp1))
                     else
                        pltx(i,k)=ampl*w(k,ii,j)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz,0.,0.,0.,'w ',plane,'w',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
               call setusv('LW',1000)
!               call curve(x,hxpl,nx)
               call frame
!
!************  x-z perturbation p cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k)=0.5*(p(k,ii,j  )-pii(k,ii,j  )  &
     &                                +p(k,ii,jp1)-pii(k,ii,jp1))
                     else
                        pltx(i,k) = p(k,ii,j)-pii(k,ii,j)
                     end if
                      pltx(i,k)=ampl*cp*t0*pltx(i,k)
                  end do
               end do
!               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'p ',plane,'p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z perturbation rr cross section
!                          
!               do i=1,nxpl
!                 ii = i+ipi-1
!                  do k=1,nz1
!                     if(mod(ii,2).eq.0.)  then
!                        pltx(i,k)=0.5
!     &                        ((rr(k,ii,j  )-rri(k,ii,j  ))*hh(ii,j  )
!     &                        +(rr(k,ii,jp1)-rri(k,ii,jp1))*hh(ii,jp1))
!                    else
!                        pltx(i,k)=(rr(k,ii,j)-rri(k,ii,j))*hh(ii,j)
!                    end if
!                  end do
!               end do
!               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.00,'r ',plane,'p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z div cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k)=ampl*0.5*(div(k,ii,j  )/rho(k,ii,j  )  &
     &                                     +div(k,ii,jp1)/rho(k,ii,jp1))
                     else
                        pltx(i,k)=ampl*div(k,ii,j)/rho(k,ii,j)  
                     end if
                  end do
               end do
!               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'d ',plane,'p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z perturbation qv cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*(.5*(qx(k,ii,j,lv)+qx(k,ii,jp1,lv))  &
     &                                                   -qvzv(k))
                     else
                        pltx(i,k) = 1000.*(qx(k,ii,j,lv)-qvzv(k))  
                     end if
                  end do
               end do
!               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qv',plane,'p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z perturbation qc cross section
!                          
               IF ( lc > 0 ) THEN
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,lc)+qx(k,ii,jp1,lc))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,lc)
                     end if
                  end do
               end do
              call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qc',plane,'p',   &
    &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
              call frame
              ENDIF
!
!************  x-z perturbation qr cross section
!                          
               IF ( lr > 0 ) THEN
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,lr)+qx(k,ii,jp1,lr))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,lr)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qr',plane,'p',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
               call frame
               ENDIF

!************  x-z perturbation qh cross section
!                          
             if ( nmoist >= lh ) then
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,lh)+qx(k,ii,jp1,lh))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,lh)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qh',plane,'p',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
               call frame

             if ( nmoist >= lhl ) then
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,lhl)+qx(k,ii,jp1,lhl))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,lhl)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'hl',plane,'p',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
               call frame
               endif

               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 0.5*(dbz(k,ii,j)+dbz(k,ii,jp1))
                     else
                        pltx(i,k) = dbz(k,ii,j)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'db',plane,'p',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
               call frame
             endif

!************  x-z perturbation qi cross section
!                          
             if ( nmoist >= li ) then
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,li)+qx(k,ii,jp1,li))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,li)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qi',plane,'p',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
               call frame
             endif

!  122          continue
!               go to 123
!
!************  y-z cross sections **********************************
!
!               i = nxc
               i = iwmax
               write(slice(2),300)  i
               slice(1) = ' I='
!
!************  y-z perturbation theta cross section
!                          
               do j=1,nypl
                  jj = j+jpi-1
                  do k=1,nz1
                     plty(j,k) = t0*t(k,i,jj)/(1.+1.61*qx(k,i,jj,lv))  &
     &                          -t0*tz(k)
                  end do
               end do
               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'t ',plane,'p',  &
     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
               call setusv('LW',1000)
!               call curve(x,hxpl,nx)
               call frame
!
!************  y-z perturbation u cross section
!                          
               do j=1,nypl
                  jj = j+jpi-1
                  jp1 = min(jj+1,ny)
                  if(jper*jj.eq.ny) jp1 = 2
                  jm1 = max(jj-1,1)
                  if(jper*jj.eq.1  )  jm1 = ny1
                  if(mod(i,2).eq.0)  then
                     jpj = j
                     jpm = jm1
                  else
                     jpj = jp1
                     jpm = jj
                  end if
                  do k=1,nz1
                     plty(j,k) = ampl*.5/sqrt(3.)  &
     &                                    *(u1(k,i,jj)+u1(k,i-1,jpj)  &
     &                                     +u3(k,i,jj)+u3(k,i-1,jpm)  &
     &                                  -2.*(u1z(k)+u1m+u3z(k)+u3m))
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'up',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation v cross section
!                          
               do j=1,nypl
                  jj = j+jpi-1
                  jm1 = jj-1
                  if(jper*jj.eq.1) jm1 = ny1
                  do k=1,nz1
                     plty(j,k) = ampl*(.5*(u2(k,i,jj)+u2(k,i,jm1))  &
     &                                    -u2z(k)-u2m) 
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'vp',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation w cross section
!                          
               do j=1,nypl
                  jj = j+jpi-1
                  do k=1,nz
                     plty(j,k)=ampl*w(k,i,jj)
                  end do
               end do
               call conplot(plty,ny,nypl,nz,0.,0.,0.,'w ',plane,'w',  &
     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
               call setusv('LW',1000)
!               call curve(x,hxpl,nx)
               call frame
!
!************  y-z perturbation p cross section
!                          
               do j=1,nypl
                  jj = j+jpi-1
                  do k=1,nz1
                     plty(j,k) = ampl*cp*t0*(p(k,i,jj)-pii(k,i,jj))
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'p ',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation rr cross section
!                          
!               do j=1,nypl
!                  jj = j+jpi-1
!                  do k=1,nz1
!                     plty(j,k)=ampl*(rr(k,i,jj)-rri(k,i,jj))*hh(i,jj)
!                  end do
!               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'r ',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z div cross section
!                          
               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=ampl*div(k,i,jj)/rho(k,i,jj)
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'d ',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation qv cross section
!                          
               do j=1,nypl
                  jj = j+jpi-1
                  do k=1,nz1
                     plty(j,k)=1000.*(qx(k,i,jj,lv)-qvzv(k))
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qv',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation qc cross section
!                          
               IF ( lc > 0 ) THEN
               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=1000.*qx(k,i,jj,lc)
                  end do
               end do
               ENDIF
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qc',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation qr cross section
!                          
               IF ( lr > 0 ) THEN
               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=1000.*qx(k,i,jj,lr)
                  end do
               end do
               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qr',plane,'p',  &
     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call curve(x,hxpl,nx)
               call frame
               ENDIF

!************  y-z perturbation qh cross section
!                          
               IF ( nmoist >= lh ) THEN
               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=1000.*qx(k,i,jj,lh)
                  end do
               end do
               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qh',plane,'p',  &
     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call curve(x,hxpl,nx)
               call frame

               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=dbz(k,i,jj)
                  end do
               end do
               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'db',plane,'p',  &
     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call curve(x,hxpl,nx)
               call frame
               
               ENDIF
!  123          continue
!cc
!               go to 124
!
!************  x-y cross sections **********************************
!
               do kk=1,3
               if(kk.eq.1)  k = 1
               if(kk.eq.2)  k = 5
               if(kk.eq.3)  k = 10
!              k = 5
               write(slice(2),300)  k
               slice(1) = ' K='
               
!
!************  x-y perturbation theta cross section
!                          
!               do j=1,ny
!                  do i=1,nx
!                     plt(i,j) = t0*t(k,i,j)/(1.+1.61*qv(k,i,j))
!     &                          -t0*tz(k)
!                  end do
!               end do

!               call trplotc(plt,nx,ny,0.,0.,0.,'t ',plane,'h ',time,
!     &              pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp,'C')
!c     &               .08,.93,.1 ,.95,xpll,xplr,ypll,yplr,dxp,dyp,'C')
!               call frame

               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = t0*.5*(t(k,ii,jj )  &
     &                                  /(1.+1.61*qx(k,ii,jj ,lv))  &
     &                                      +t(k,ii,jp1)  &
     &                                  /(1.+1.61*qx(k,ii,jp1,lv)))  &
     &                            -t0*tz(k)
                     else
                        plt(i,j) = t0*t(k,ii,jj)/(1.+1.61*qx(k,ii,jj,lv))  &
     &                            -t0*tz(k)
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'t ',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame
!
!************  x-y perturbation u1 cross section
!                          
!               do j=1,ny
!                  do i=1,nx
!                     plt(i,j) = u1(k,i,j)
!                  end do
!               end do
!
!               call trplot(plt,nx,ny,0.,0.,0.,'u1',plane,'u1',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y perturbation u3 cross section
!                          
!               do j=1,ny
!                  do i=1,nx
!                     plt(i,j) = u3(k,i,j)
!                  end do
!               end do
!
!               call trplot(plt,nx,ny,0.,0.,0.,'u3',plane,'u3',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y perturbation u cross section
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
                     im1=i-1
                     if(iper*i.eq.1) im1 = nx1
                     plt(i,j) = ampl*.5/sqrt(3.)  &
     &                                    *(u1(k,i,j)+u1(k,im1,jpj)  &
     &                                     +u3(k,i,j)+u3(k,im1,jpm)  &
     &                                  -2.*(u1z(k)+u1m+u3z(k)+u3m))
                  end do
               end do
!               call trplot(plt,nx,ny,0.,0.,0.,'u ',plane,'h ',time,
!     &                pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dzp)
!               call setusv('LW',1000)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-y perturbation u2 cross section
!                          
                do j=1,ny
                  do i=1,nx
                     plt(i,j) = u2(k,i,j)
                  end do
               end do

               call trplot(plt,nx,ny,0.,0.,0.,'u2',plane,'u2',time, &
                            pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
                call frame
!
!************  x-y perturbation w cross section
!                          
!               do j=1,ny
!                  do i=1,nx
!                     plt(i,j) = .5*(w(k,i,j)+w(k+1,i,j))
!                  end do
!               end do

!               call trplotc(plt,nx,ny,0.,0.,0.,'w field ','h ',time,
!c     &              pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp,'C')
!     &               .08,.93,.1 ,.95,xpll,xplr,ypll,yplr,dxp,dyp,'C')
!               call frame

               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = .25*(w(k,ii,jj )+w(k+1,ii,jj )  &
     &                                 +w(k,ii,jp1)+w(k+1,ii,jp1))
                     else
                        plt(i,j) = .5*(w(k,ii,jj)+w(k+1,ii,jj))
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'w ',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame
!
!************  x-y perturbation p cross section
!                          
               do j=1,ny
                  do i=1,nx
                     plt(i,j) = ampl*cp*t0*(p(k,i,j)-pii(k,i,j))
                  end do
               end do

!               call trplot(plt,nx,ny,0.,0.,0.,'p ',plane,'h ',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y div cross section
!                          
               do j=1,ny
                  do i=1,nx
                     plt(i,j)=ampl*div(k,i,j)/rho(k,i,j)
                  end do
               end do

!               call trplot(plt,nx,ny,0.,0.,0.,'div fld ','h ',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y div cross section
!                          
!               do j=1,ny
!              jp1 =j+1
!              if(j.eq.ny)  jp1 = 2
!                  do i=1,nx
!                     if(mod(i,2).eq.0)  then
!                     plt(i,j)=ampl*(div(k,i,j)+div(k,i,jp1))
!     &                            /(rho(k,i,j)+rho(k,i,jp1))
!                    else
!                     plt(i,j)=ampl*div(k,i,j)/rho(k,i,j)
!                    end if
!                  end do
!               end do

!               call conplot(plt,nx,ny,0.,0.,0.,'dv',plane,'p',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y perturbation qv cross section
!                          
               do j=1,ny
                  do i=1,nx
                     plt(i,j) = 1000.*(qx(k,i,j,lv)-qvzv(k))
                  end do
               end do

!               call trplot(plt,nx,ny,0.,0.,0.,'qv',plane,'h ',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y perturbation qc cross section
!                          
               IF ( lc > 0 ) THEN
               do j=1,ny
                  do i=1,nx
                     plt(i,j) = 1000.*qx(k,i,j,lc)
                  end do
               end do
               ENDIF

!               call trplot(plt,nx,ny,0.,0.,0.,'qc',plane,'h ',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y perturbation qr cross section
!                          
!               do j=1,ny
!                  do i=1,nx
!                     plt(i,j) = 1000.*qr(k,i,j)
!                  end do
!               end do

!               call trplotc(plt,nx,ny,0.,0.,0.,'qr',plane,'h ',time,
!     &              .08,.93,.1,.95,xpll,xplr,ypll,yplr,dxp,dyp,'C')
!               call frame
               
               IF ( lr > 0 ) THEN
               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = 1000.*.5*(qx(k,ii,jj,lr)+qx(k,ii,jp1,lr))
                     else
                        plt(i,j) = 1000.*qx(k,ii,jj,lr)
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'qr',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame
               ENDIF

              IF ( nmoist >= lh ) THEN
               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = 1000.*.5*(qx(k,ii,jj,lh)+qx(k,ii,jp1,lh))
                     else
                        plt(i,j) = 1000.*qx(k,ii,jj,lh)
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'qh',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame


               IF ( nmoist >= lhl ) THEN
               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = 1000.*.5*(qx(k,ii,jj,lhl)+qx(k,ii,jp1,lhl))
                     else
                        plt(i,j) = 1000.*qx(k,ii,jj,lhl)
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'hl',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame
               ENDIF

               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = 0.5*(dbz(k,ii,jj)+dbz(k,ii,jp1))
                     else
                        plt(i,j) = dbz(k,ii,jj)
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'db',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame
               
               ENDIF

               end do

!  124          continue

        endif ! doplot
        end if

     END IF

