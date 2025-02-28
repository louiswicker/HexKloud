c
c-----------------------------------------------------------------------
c                                                                     72
      subroutine rhs_s( s,s1,fs,ww,ru1,ru2,ru3,rho,dmp,dts,dtsa,
     $           rdz,xnus,xnusz,nz1,nx,ny,iper,jper,
     &           s0,n1,n2,n3,flux1,flux2,flux3,fluxz,order)

      implicit none

      integer nz1,nx,ny,iper,jper,n1,n2,n3
      real    s    (nz1,  nx,ny),s1   (nz1  ,nx,ny),ru1  (nz1,0:nx,ny)
     &       ,fs   (nz1,  nx,ny),ww   (nz1+1,nx,ny),ru2  (nz1,nx,0:ny)
     &       ,s0   (n1 ,  n2,n3),rho  (nz1  ,nx,ny),ru3  (nz1,0:nx,ny)
     &       ,flux1(nz1,0:nx,ny),flux2(nz1,nx,0:ny),flux3(nz1,0:nx,ny)
     &       ,fluxz(0:nz1,nx,ny)
    
      real    dts,dtsa,rdz,xnus,cr,xnusz,prandtl,dmp(nz1)

      integer nx1,ny1,nz,i,j,k,ip1,im1,jp1,jm1,nz2,jpj,jpp,jpm,jmm
     &       ,nscratch,ip2,im2,jp2,jm2,ip3,jp3,ii,jj
      character*6 order

      nx1 = nx-1
      ny1 = ny-1
      nz = nz1+1
      nz2 = nz1-1

      prandtl = 3.0
c
c        right hand side of s equation
c
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
c
c  vertical flux calculation
c
			fluxz(0  ,i,j) = 0.
			fluxz(nz1,i,j) = 0.
			fluxz(1  ,i,j) = 0.5*ww(2  ,i,j)*(s(2  ,i,j)+s(1  ,i,j))
			fluxz(nz2,i,j) = 0.5*ww(nz1,i,j)*(s(nz1,i,j)+s(nz2,i,j))

		    if (order.eq.'second'.or.order.eq.'fourth') then

			   do k=2,nz1-2  
			      fluxz(k,i,j) = .5*ww(k+1,i,j)*(s(k+1,i,j)+s(k,i,j))
			   end do
		    else
			   do k=2,nz1-2  
                  fluxz(k,i,j) = (1./12.)*(ww(k+1,i,j) *(
     &                7.*(s(k+1,i,j)+s(k,i,j))-(s(k+2,i,j)+s(k-1,i,j)))
     &                                +abs(ww(k+1,i,j))*(
     &               -3.*(s(k+1,i,j)-s(k,i,j))+(s(k+2,i,j)-s(k-1,i,j))))

ccc			      fluxz(k,i,j) = .5*ww(k+1,i,j)*(s(k+1,i,j)+s(k,i,j))
			   end do
			end if
c
c  horizontal flux calculations
c
			if(order.eq.'second')  then
			   do k=1,nz1
				  flux1(k,i,j) = .5*ru1(k,i,j)*(s(k,i,j)+s(k,ip1,jpm))
				  flux2(k,i,j) = .5*ru2(k,i,j)*(s(k,i,j)+s(k,i  ,jp1))
				  flux3(k,i,j) = .5*ru3(k,i,j)*(s(k,i,j)+s(k,ip1,jpj))
			   end do
			else if(order.eq.'third '.or.order.eq.'fourth')  then
			   do k=1,nz1
				  flux1(k,i,j) =  1./12.*ru1(k,i  ,j  )
     &                                *(7.*(s(k,ip1,jpm)+s(k,i  ,j  ))
     &                                    -(s(k,ip2,jm1)+s(k,im1,jpj)))
				  flux2(k,i,j) =  1./12.*ru2(k,i  ,j  )
     &                                *(7.*(s(k,i  ,jp1)+s(k,i  ,j  ))
     &                                    -(s(k,i  ,jp2)+s(k,i  ,jm1)))
				  flux3(k,i,j) =  1./12.*ru3(k,i  ,j  )
     &                                *(7.*(s(k,ip1,jpj)+s(k,i  ,j  ))
     &                                    -(s(k,ip2,jp1)+s(k,im1,jpm)))
			   end do
			   if(order.eq.'third ')  then
			      do k=1,nz1
				     flux1(k,i,j) =  flux1(k,i,j)+1./12.*abs(ru1(k,i,j))
     &                               *(-3.*(s(k,ip1,jpm)-s(k,i  ,j  ))
     &                                    +(s(k,ip2,jm1)-s(k,im1,jpj)))
				     flux2(k,i,j) =  flux2(k,i,j)+1./12.*abs(ru2(k,i,j))
     &                               *(-3.*(s(k,i  ,jp1)-s(k,i  ,j  ))
     &                                    +(s(k,i  ,jp2)-s(k,i  ,jm1)))
				     flux3(k,i,j) =  flux3(k,i,j)+1./12.*abs(ru3(k,i,j))
     &                               *(-3.*(s(k,ip1,jpj)-s(k,i  ,j  ))
     &                                    +(s(k,ip2,jp1)-s(k,im1,jpm))) 
			      end do
			   end if
			else if(order.eq.'fifth '.or.order.eq.'sixth ')  then
			   do k=1,nz1
                  flux1(k,i,j) = (1./60.)*ru1(k,i  ,j)
     &                              *(37.*(s(k,ip1,jpm)+s(k,i  ,j  ))
     &                                -8.*(s(k,ip2,jm1)+s(k,im1,jpj)) 
     &                                   +(s(k,ip3,jmm)+s(k,im2,jp1)))
                  flux2(k,i,j) = (1./60.)*ru2(k,i  ,j)
     &                              *(37.*(s(k,i  ,jp1)+s(k,i  ,j  ))
     &                                -8.*(s(k,i  ,jp2)+s(k,i  ,jm1)) 
     &                                   +(s(k,i  ,jp3)+s(k,i  ,jm2)))
                  flux3(k,i,j) = (1./60.)*ru3(k,i  ,j)
     &                              *(37.*(s(k,ip1,jpj)+s(k,i,j))
     &                                -8.*(s(k,ip2,jp1)+s(k,im1,jpm)) 
     &                                   +(s(k,ip3,jpp)+s(k,im2,jm1)))
			   end do
			   if(order.eq.'fifth ')  then
			      do k=1,nz1
				     flux1(k,i,j) =flux1(k,i,j)-(1./60.)*abs(ru1(k,i,j))
     &                                  *((s(k,ip3,jmm)-s(k,im2,jp1))
     &                                -5.*(s(k,ip2,jm1)-s(k,im1,jpj)) 
     &                               +10.*(s(k,ip1,jpm)-s(k,i  ,j  )))
                     flux2(k,i,j) =flux2(k,i,j)-(1./60.)*abs(ru2(k,i,j))
     &                                  *((s(k,i  ,jp3)-s(k,i  ,jm2))
     &                                -5.*(s(k,i  ,jp2)-s(k,i  ,jm1)) 
     &                               +10.*(s(k,i  ,jp1)-s(k,i  ,j  )))
                     flux3(k,i,j) =flux3(k,i,j)-(1./60.)*abs(ru3(k,i,j))
     &                                  *((s(k,ip3,jpp)-s(k,im2,jm1))
     &                                -5.*(s(k,ip2,jp1)-s(k,im1,jpm)) 
     &                               +10.*(s(k,ip1,jpj)-s(k,i  ,j  )))
			      end do
			   end if
			end if
		 end do
	  end do
c
c   Flux divergence calculation
c
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
			   fs(k,i,j) = - dtsa*(flux1(k,i,j)-flux1(k,im1,jpj)
     &                           +flux2(k,i,j)-flux2(k,i  ,jm1)
     &                           +flux3(k,i,j)-flux3(k,im1,jpm))
     &                 - dts*rdz*(fluxz(k,i,j)-fluxz(k-1,i,  j))
			end do
		 end do
	  end do
c
c        horizontal mixing and Rayleigh damping terms
c
         do j = 1,ny
            jp1 = min(j+1,ny)
            if(jper*j.eq.ny )  jp1 = 2
            jp2 = min(jp1+1,ny)
            if(jper*j.eq.ny1)  jp2 = 2
            jm1 = max(j-1,1)
            if(jper*j.eq.1  )  jm1 = ny1
            jm2 = jm1-1
            if(jper*j.eq.2  )  jm2 = ny1
            do i = 1,nx
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
			   do k = 1,nz1
			      fs(k,i,j) = fs(k,i,j) + prandtl*xnus*rho(k,i,j)
     &                       *(s1(k,ip1,jpj)+s1(k,im1,jpm)
     &                        +s1(k,ip1,jpm)+s1(k,im1,jpj)
     &                        +s1(k,i  ,jp1)+s1(k,i  ,jm1)-6.*s1(k,i,j))
			   end do
			   if(n3.eq.ny)  then
                  do k=1,nz1
			         fs(k,i,j) = fs(k,i,j) - prandtl*xnus*rho(k,i,j)
     &                       *(s0(k,ip1,jpj)+s0(k,im1,jpm)
     &                        +s0(k,ip1,jpm)+s0(k,im1,jpj)
     &                        +s0(k,i  ,jp1)+s0(k,i  ,jm1)-6.*s0(k,i,j)) 
                  end do
               end if
            end do
         end do
c
c        vertical mixing terms
c
		 do j = 1,ny
		    jj = min(j,n3)
	        do i = 1,nx
	  	       ii = min(i,n2)
			   if(n1.eq.nz1)  then
                  do k = 2,nz2
                     fs(k,i,j) = fs(k,i,j)+prandtl*xnusz*rho(k,i,j)
     &                   *((s1(k+1,i ,j )-2.*s1(k,i ,j )+s1(k-1,i ,j ))
     &                    -(s0(k+1,ii,jj)-2.*s0(k,ii,jj)+s0(k-1,ii,jj)))
                  end do
			   else
                  do k=2,nz2
                     fs(k,i,j)=fs(k,i,j)+prandtl*xnusz*rho(k,i,j)
     &                   *((s1(k+1,i ,j )-2.*s1(k,i ,j )+s1(k-1,i ,j )))
               end do
			   end if
            end do
         end do
c
c        Rayleigh damping layer
c		 
		 if(n3.eq.ny)  then
            do j=1,ny
               do i=1,nx
                  do k=1,nz1
                     fs(k,i,j) = fs(k,i,j)
     &                        - dmp(k)*rho(k,i,j)*(s1(k,i,j)-s0(k,i,j)) 
                  end do
               end do
            end do
         end if


      return
      end
c
c------------------------------------------------------------------------
c
      subroutine rhs_rho (fr,ru1,ru2,ru3,ww,dts,dtsa,rdz
     &                      ,nz1,nx,ny,iper,jper      )

      implicit none

      integer nz1,nx,ny,iper,jper
      real    fr (nz1  ,nx,ny),ru1(nz1,0:nx,ny),ru2(nz1,nx,0:ny),
     &        ww (nz1+1,nx,ny),ru3(nz1,0:nx,ny)
      real    dts,dtsa,rdz

      integer i,j,k,im1,jm1,nx1,ny1,jp1,jpj,jpm


       nx1 = nx-1
       ny1 = ny-1
c
c        right hand side of rho equation
c
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
               im1 = i-1
               if(iper*i.eq.1)  im1 = nx1
               do k=1,nz1
                  fr(k,i,j) = - dtsa*(ru3(k,i,j)-ru3(k,im1,jpm)
     &                               +ru1(k,i,j)-ru1(k,im1,jpj)
     &                               +ru2(k,i,j)-ru2(k,i  ,jm1))
     &                        - dts*rdz*(ww(k+1,i,j)-ww(k,i  ,j))
               end do
            end do
         end do

         return
         end
c
c------------------------------------------------------------------------
