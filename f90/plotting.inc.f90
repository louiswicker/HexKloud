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
         
         if((kkk.ge.ip).or. (nit.eq.0)) then      
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
            write(6,*) '       percent mass change =',rtot

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
     &                      (t(k,ii,j  )/(1.+1.61*qx(k,ii,j  ,1))  &
     &                      +t(k,ii,jp1)/(1.+1.61*qx(k,ii,jp1,1)))
                     else
                        pltx(i,k) = t(k,ii,j)/(1.+1.61*qx(k,ii,j,1))
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
                        pltx(i,k) = 1000.*(.5*(qx(k,ii,j,1)+qx(k,ii,jp1,1))  &
     &                                                   -qvzv(k))
                     else
                        pltx(i,k) = 1000.*(qx(k,ii,j,1)-qvzv(k))  
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
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,2)+qx(k,ii,jp1,2))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,2)
                     end if
                  end do
               end do
!               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qc',plane,'p',
!     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  x-z perturbation qr cross section
!                          
               do i=1,nxpl
                  ii = i+ipi-1
                  do k=1,nz1
                     if(mod(ii,2).eq.0.)  then
                        pltx(i,k) = 1000.*.5*(qx(k,ii,j,3)+qx(k,ii,jp1,3))
                     else
                        pltx(i,k) = 1000.*qx(k,ii,j,3)
                     end if
                  end do
               end do
               call conplot(pltx,nx,nxpl,nz1,0.,0.,0.,'qr',plane,'p',  &
     &              time,pxl,pxr,pzl,zptop,xpll,xplr,zplb,zplt,dxp,dzp)
!               call curve(x,hxpl,nx)
               call frame

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
                     plty(j,k) = t0*t(k,i,jj)/(1.+1.61*qx(k,i,jj,1))  &
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
                     plty(j,k)=1000.*(qx(k,i,jj,1)-qvzv(k))
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qv',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation qc cross section
!                          
               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=1000.*qx(k,i,jj,2)
                  end do
               end do
!               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qc',plane,'p',
!     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!c               call curve(x,hxpl,nx)
!               call frame
!
!************  y-z perturbation qr cross section
!                          
               do k=1,nz1
                  do j=1,nypl
                     jj = j+jpi-1
                     plty(j,k)=1000.*qx(k,i,jj,3)
                  end do
               end do
               call conplot(plty,ny,nypl,nz1,0.,0.,0.,'qr',plane,'p',  &
     &              time,pyl,pyr,pzl,zptop,ypll,yplr,zplb,zplt,dyp,dzp)
!               call curve(x,hxpl,nx)
               call frame
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
     &                                  /(1.+1.61*qx(k,ii,jj ,1))  &
     &                                      +t(k,ii,jp1)  &
     &                                  /(1.+1.61*qx(k,ii,jp1,1)))  &
     &                            -t0*tz(k)
                     else
                        plt(i,j) = t0*t(k,ii,jj)/(1.+1.61*qx(k,ii,jj,1))  &
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

!              call trplot(plt,nx,ny,0.,0.,0.,'u2',plane,'u2',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
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
                     plt(i,j) = 1000.*(qx(k,i,j,1)-qvzv(k))
                  end do
               end do

!               call trplot(plt,nx,ny,0.,0.,0.,'qv',plane,'h ',time,
!     &                     pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
!               call frame
!
!************  x-y perturbation qc cross section
!                          
               do j=1,ny
                  do i=1,nx
                     plt(i,j) = 1000.*qx(k,i,j,2)
                  end do
               end do

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
               
               do j=1,nypl
                  jj = j+jpi-1
                  do i=1,nxpl
                     ii = i+ipi-1
                     if(mod(ii,2).eq.0.)  then
                        jp1 =min(jj+1,ny)
                        if(jper*jj.eq.ny)  jp1 = 2
                        plt(i,j) = 1000.*.5*(qx(k,ii,jj,3)+qx(k,ii,jp1,3))
                     else
                        plt(i,j) = 1000.*qx(k,ii,jj,3)
                     end if
                  end do
               end do
               call conplot(plt,nx,nxpl,nypl,0.,0.,0.,'qr',plane,'h ',  &
     &              time,pxl,pxr,pyl,pyr,xpll,xplr,ypll,yplr,dxp,dyp)
               call frame

               end do

!  124          continue

        end if
         end if

