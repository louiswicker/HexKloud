c
c-----------------------------------------------------------------------
c                                                                     72
c
      subroutine init_sound(t,qv,u1z,u2z,u3z,zg,n)

      integer n
      real    t(n),qv(n),u1z(n),u2z(n),u3z(n),zg(n)

      pi     = 4.*atan(1.)
      nz1    = n
      ZTR    = 12000.
      ALPHA  = 4.
cc      vs     = 12.0
      vs     = 30.0
cc      vs     = 0.0
      angle  = 0.
c      angle  = 45.
      THETAR = 343.
      TTR    = 213.
      QV0M   = 0.014
c      QV0M   = 0.013
      THETAS = 300.5
cc      ZTS    = 2500.
      ZTS    = 5000.
      CC1 = 1.
      CC2 = 0.
      umax  = vs*cos(angle/180.*pi)
      vmax  = vs*sin(angle/180.*pi)
      u1max = .5*sqrt(3.)*umax-.5*vmax
      u3max = .5*sqrt(3.)*umax+.5*vmax
      u2max = vmax

      do k=1,nz1
         z = zg(k)
         ztp = cc1*z + cc2*z**2
         if(ztp .gt. ztr) then
            t(k) = thetar*exp(9.8*(ztp-ztr)/(1003.*ttr))
            qv(k) = 0.25
         else
            t(k) = 300.+43.*(ztp/ztr)**1.25
            qv(k) = (1.-0.75*(ztp/ztr)**1.25)
            if(t(k).lt.thetas) t(k)=thetas
         end if

ccc         t(k)=300.
ccc         qv(k) = 0.

         u1z(k) = u1max*min(ztp/zts,1.)
         u2z(k) = u2max*min(ztp/zts,1.)
         u3z(k) = u3max*min(ztp/zts,1.)
      end do

      write(6,*) ' initial conditions '
c
      do k=1,nz1
c         write  (6,10) k,zg(k),t(k),u1z(k),u2z(k),u3z(k),qv(k)
   10    format (1x,i3,f8.1,f7.2,3f6.2,f8.5)
      end do
      return
      end

c
