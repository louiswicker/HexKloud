      subroutine wplot(wmax,waxis,wmplt,ip,lsteps,nsteps)
      character*25 pltlab
      pltlab = 'Maximum Vertical Velocity'
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchlq(.5,.97,pltlab,.02,0,0)
      CALL SET(.05,.95,.05,.95,1.,float(nsteps+1),0.,wmplt,1)
      CALL PERIM(nsteps/ip,1,5,5)
      CALL CURVE(waxis,wmax,lsteps)
!      CALL DASHDB(255)
!      CALL CURVED()
      CALL FRAME    
      return
      end
!
!----------------------------------------------------------------------
!
      subroutine conplot(array,n1,n1p,n2,cmn,cmx,cis,field,plane,stag,  &
     &                  time,pl,pr,pb,pt,al,ar,ab,at,dx,dz)
      real array(n1,n2)
      real*4 aln,arn,abn,atn
      character*8 label(3)
      character*1 stag
      character*24 pltlab
      character*2 field
      character*6 plane
      equivalence (pltlab,label)
      COMMON /CONRE4/ ISIZEL     ,ISIZEM     ,ISIZEP     ,NREP       ,  &
     &                NCRT       ,ILAB       ,NULBLL     ,IOFFD      ,  &
     &                EXT        ,IOFFM      ,ISOLID     ,NLA        ,  &
     &                NLM        ,XLT        ,YBT        ,SIDE

!      ilab = 0
      itime = Nint(time)
      pltlab(1:2) = field
      pltlab(3:8) = plane
      label(2) = ' at t = '
      write(label(3),500)   itime
  500 format(i8)
      call setusv('LW',1000)       
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchlq(.5,.97,pltlab,.02,0,0)

      aln = al/dx 
      arn = ar/dx 
      abn = ab/dz 
      atn = at/dz 
      if(pr-pl.gt..01.and.pt-pb.gt..01)  then
         call set(pl,pr,pb,pt,aln,arn,abn,atn,1)
      else
         call set(.05,.95,.05,.95,aln,arn,abn,atn,1)
      end if
      itick1 = nint(arn-aln+.01)
      itick2 = nint(atn-abn+.01)
      call setusv('LW',1000)
      call perim(itick1,1,itick2,1)

!      call setusv('LW',2000)       
      call conrec(array,n1,n1p,n2,cmn,cmx,cis,-1,-638,-922)
!      call conrec(array,n1,n1p,n2,cmn,cmx,cis,-1,-1,0)

      return
      end
!
!----------------------------------------------------------------------
!
      subroutine trplot(array,n1,n2,cmn,cmx,cis,field,plane,stag,  &
     &                  time,pl,pr,pb,pt,al,ar,ab,at,dx,dz)
      parameter (lrwk=1500,liwk=2500)
      dimension rwrk(lrwk),iwrk(liwk),array(n1,n2)
      character*2  stag
      character*8  label(3)
      character*24 pltlab
      character*2 field
      character*6 plane
      equivalence (pltlab,label)
 
!     put header label at top of frame

      itime = nint(time)
      pltlab(1:2) = field
      pltlab(3:8) = plane
      label(2) = ' at t = '
      write(label(3),500)   itime
  500 format(i8)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.97,pltlab,.02,0.,0.)

!     set the min, max, and contour intercal 

      if(cmx-cmn.gt.0.)  then
         call cpsetr('cmn',cmn)
         call cpsetr('cmx',cmx)
      else
         call cprset
      end if
      call cpsetr('cis',cis)

!     call set and draw perimeter

      if(pr-pl.gt..01.and.pt-pb.gt..01)  then
         call set(pl,pr,pb,pt,al,ar,ab,at,1)
      else
         call set(.05,.95,.05,.95,al,ar,ab,at,1)
      end if
      itick1 = (ar-al)/dx+.01 
      itick2 = (at-ab)/dz +.01

!      call setusv('lw',2000)

      call perim(itick1,1,itick2,1)
 
!     use previous set call for contour plot

      call cpseti('set',0)

!     turn off plotting highs and lows

      call cpsetc('hlt - high/low label text',' ')

!     set special value

      call cpsetr('spv',99999.)

!     initialization for countour plot

!      call setusv('lw',2000)

      call cprect(array,n1,n1,n2,rwrk,lrwk,iwrk,liwk)

!     draw negative contours with dashed lines

      call cppkcl(array, rwrk, iwrk)
      call cpgeti ('ncl - number of contour levels',nclv)
      do  iclv=1,nclv
         call cpseti ('pai - parameter array index',iclv)
         call cpgetr( 'clv - contour level value', clv)
          if(clv .lt. 0.)then
!             call cpseti('cld - contour level dash pattern',52428)
             call cpsetc('cld - contour level dash pattern',  &
     &                ' ''''$$$''''$$$''''$$$ ')
         endif
      end do

      call setusv('LW',1000)

!     set the transformation parameter

      if(stag.eq.'h ')  then
         imap = 1
      else if(stag.eq.'u3')  then
         imap = 2
      else if(stag.eq.'u1')  then
         imap = 3
      else if(stag.eq.'u2')  then
         imap = 4
      end if
      call cpseti('map',imap)

!     draw contour lines

!      call setusv('LW',3000)

!     CALL TO TURN OFF LINE LABELS
!      CALL CPSETI('LLP',0)
      
      call cpcldr(array,rwrk,iwrk)

!     draw line labels and informational label

!      call setusv('LW',3000)

      call cplbdr(array,rwrk,iwrk)

!     check amount of workspace used

      call cpgeti('iwu',iiwu)
      call cpgeti('rwu',irwu)
      if(irwu.gt.lrwk.or.iiwu.gt.liwk)  then
         write(6,*) ' insufficient workspace for conpack'
         write(6,*) ' liwk must be at least ',iiwu
         write(6,*) ' lrwk must be at least ',irwu
      end if
      return
      end
!
!---------------------------------------------------------------------
!
      subroutine cpmpxy(imap,xinp,yinp,xotp,yotp)
!      parameter (nx=43,ny=37)
!      parameter (nx=61,ny=53)
!      parameter (nx=181,ny=157)
      parameter (nx=91,ny=79)
!      parameter (nx=47,ny=40)
!      parameter (nx=121,ny=105)
!      parameter (nx=181,ny=53)
!      parameter (nx=101,ny=5)
!      parameter (nx=5,ny=101)

      common /grid/ xh(nx,ny),xu1(nx,ny),xu2(nx,ny),xu3(nx,ny),  &
     &              yh(nx,ny),yu1(nx,ny),yu2(nx,ny),yu3(nx,ny)

      if(imap.eq.1)  then
         im  = int(xinp)
!         ip  = min(im+1,nx)
         ip  = im+1
         if(ip.gt.nx)  ip=ip-nx+1

         fri = xinp-im
         km  = int(yinp)
!         kp  = min(km+1,ny)
         kp  = km+1
         if(kp.gt.ny)  kp=kp-ny+1

         frk = yinp-km
         xm  = xh(im,km)+fri*(xh(ip,km)-xh(im,km))
         xp  = xh(im,kp)+fri*(xh(ip,kp)-xh(im,kp))
         xotp= xm+frk*(xp-xm)
         ym  = yh(im,km)+fri*(yh(ip,km)-yh(im,km))
         yp  = yh(im,kp)+fri*(yh(ip,kp)-yh(im,kp))
         yotp= ym+frk*(yp-ym)
      else if(imap.eq.2)  then
         im  = int(xinp)
!         ip  = min(im+1,nx)
         ip  = im+1
         if(ip.gt.nx)  ip=ip-nx+1

         fri = xinp-im
         km  = int(yinp)
!         kp  = min(km+1,ny)
         kp  = km+1
         if(kp.gt.ny)  kp=kp-ny+1

         frk = yinp-km
         xm  = xu3(im,km)+fri*(xu3(ip,km)-xu3(im,km))
         xp  = xu3(im,kp)+fri*(xu3(ip,kp)-xu3(im,kp))
         xotp= xm+frk*(xp-xm)
         ym  = yu3(im,km)+fri*(yu3(ip,km)-yu3(im,km))
         yp  = yu3(im,kp)+fri*(yu3(ip,kp)-yu3(im,kp))
         yotp= ym+frk*(yp-ym)
      else if(imap.eq.3)  then
         im  = int(xinp)
!         ip  = min(im+1,nx)
         ip  = im+1
         if(ip.gt.nx)  ip=ip-nx+1

         fri = xinp-im
         km  = int(yinp)
!         kp  = min(km+1,ny)
         kp  = km+1
         if(kp.gt.ny)  kp=kp-ny+1

         frk = yinp-km
         xm  = xu1(im,km)+fri*(xu1(ip,km)-xu1(im,km))
         xp  = xu1(im,kp)+fri*(xu1(ip,kp)-xu1(im,kp))
         xotp= xm+frk*(xp-xm)
         ym  = yu1(im,km)+fri*(yu1(ip,km)-yu1(im,km))
         yp  = yu1(im,kp)+fri*(yu1(ip,kp)-yu1(im,kp))
         yotp= ym+frk*(yp-ym)
      else if(imap.eq.4)  then
         im  = int(xinp)
!         ip  = min(im+1,nx)
         ip  = im+1
         if(ip.gt.nx)  ip=ip-nx+1

         fri = xinp-im
         km  = int(yinp)
!         kp  = min(km+1,ny)
         kp  = km+1
         if(kp.gt.ny)  kp=kp-ny+1

         frk = yinp-km
         xm  = xu2(im,km)+fri*(xu2(ip,km)-xu2(im,km))
         xp  = xu2(im,kp)+fri*(xu2(ip,kp)-xu2(im,kp))
         xotp= xm+frk*(xp-xm)
         ym  = yu2(im,km)+fri*(yu2(ip,km)-yu2(im,km))
         yp  = yu2(im,kp)+fri*(yu2(ip,kp)-yu2(im,kp))
         yotp= ym+frk*(yp-ym)
      end if
      return
      end
!
!***********************************************************************
!
      SUBROUTINE TRPLOTC(ARRAY,N1,N2,CMN,CMX,CIS,FIELD,PLANE,STAG,  &
     &           TIME,PL,PR,PB,PT,AL,AR,AB,AT,DX,DZ,COLOR)
!
        PARAMETER (NAMA=400000,LRWK=1000,LIWK=5000,NCRA=20000)
!
! Declare required data arrays and workspace arrays.
!
        DIMENSION ARRAY(N1,N2),RWRK(LRWK),IWRK(LIWK),IAMA(NAMA)
        DIMENSION IASF(13)
        DIMENSION XCRA(NCRA),YCRA(NCRA)
        DIMENSION IAIA(10),IGIA(10)
!
! Declare arrays to hold the list of indices and the list of labels
! required by the label-bar routine.
!
        DIMENSION LIND(14)

        CHARACTER*10 LLBS(15)
!
      CHARACTER*1  COLOR
      character*2  stag
      character*8  label(3)
      character*24 pltlab
      character*2 field
      character*6 plane
      equivalence (pltlab,label)
!
! Declare the routine which will color the areas.
!
        EXTERNAL COLRAM
!
! Define an array for GKS aspect source flags.
!
        DATA IASF / 13*1 /
!
! Define the list of indices required by the label-bar routine.
!
        DATA LIND / 2,3,4,5,6,7,8,9,10,11,12,13,14,15 /

!        Define color indices.
!
         CALL DFCLRS(1)
!
!     PUT HEADER LABEL AT TOP OF FRAME

      ITIME = TIME
      pltlab(1:2) = field
      pltlab(3:8) = plane
      LABEL(2) = ' at t = '
      WRITE(LABEL(3),500)   ITIME
  500 FORMAT(I8)
      CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL GSPLCI(2)
      CALL PLCHHQ(.5,.97,PLTLAB,.02,0.,0.)
!
!     SET LOCATION OF CONTOUR PLOT

      IF(PR-PL.GT..01.AND.PT-PB.GT..01)  THEN
         CALL SET(PL,PR,PB,PT,AL,AR,AB,AT,1)
      ELSE
         CALL SET(.05,.95,.05,.95,AL,AR,AB,AT,1)
      END IF

!     USE PREVIOUS SET CALL FOR CONTOUR PLOT

      CALL CPSETI('SET',0)

!     SET THE TRANSFORMATION PARAMETER

      IF(STAG.EQ.'P')  THEN
         IMAP = 1
      ELSE IF(STAG.EQ.'W')  THEN
         IMAP = 2
      ELSE
         IMAP = 3
      END IF
      CALL CPSETI('MAP',IMAP)
!
! Initialize the drawing of the contour plot.
!
        CALL CPRECT (ARRAY,N1,N1,N2,RWRK,LRWK,IWRK,LIWK)
!
!     SET THE MIN AND MAX

      IF(CMX-CMN.GT.0.)  THEN
         CALL CPSETR('CMN',CMN)
         CALL CPSETR('CMX',CMX)
      END IF

!     SET SPECIAL VALUE

      CALL CPSETR('SPV',99999.)

      IF(COLOR.EQ.'C')  THEN

!        Tell CONPACK to use 13 contour levels, splitting the range into
!        14 equal bands, one for each of the 14 colors available.
!
      ICLS = 13

!     ADJUST CONTOUR LEVSLS IF NOT SPECIFIED

      IF(CMX-CMN.GT.0.)  THEN
         IF(CIS.EQ.0.)  THEN
            ICLS = 9
            CIS = (CMX-CMN)/FLOAT(ICLS+1)
         ELSE
            ICLS = NINT((CMX-CMN)/CIS)-1
         END IF
      else
         if(cis.EQ.0.)  then
            icls = 6
         end if
      END IF
!      WRITE(6,*) ICLS,CIS
!
      CALL CPSETR('CIS',CIS)
!
       CALL CPSETI ('CLS - CONTOUR LEVEL SELECTOR',ICLS)
!
!        Turn off the clipping indicator.
!
         CALL GSCLIP (0)
!
!        Set all aspect source flags to "individual".
!
         CALL GSASF (IASF)
!
!        Force solid fill.
!
         CALL GSFAIS (1)
!


!        Disallow the trimming of trailing zeroes.
!
!         CALL CPSETI ('NOF - NUMERIC OMISSION FLAGS',0)
!
!        Initialize the area map and put the contour lines into it.
!
         CALL ARINAM (IAMA,NAMA)
         CALL CPCLAM (ARRAY,RWRK,IWRK,IAMA)
!
!        Color the map.
!
         CALL ARSCAM (IAMA,XCRA,YCRA,NCRA,IAIA,IGIA,10,COLRAM)
!
!        Put black contour lines over the colored map.
!
         CALL GSPLCI (0)
         CALL CPCLDR (ARRAY,RWRK,IWRK)
         CALL GSPLCI (1)
!
!        Draw a label bar for the plot, relating colors to values.
!
         call cpgeti('ncl - number of contour levels',ncl)
         if(ncl.eq.0)  go to 100
         call cpgetr('ciu - contour interval used',ciu)
         CALL CPGETR ('ZMN',ZMIN)
         CALL CPGETR ('ZMX',ZMAX)
         ISKP = 0
         NBOX = NCL+1
         IF(CMX-CMN.gt.0.) THEN
            IF(CMN.LT.ZMIN.AND.CMX.LE.ZMAX)   NBOX = NCL
            IF(CMN.GE.ZMIN.AND.CMX.GT.ZMAX)   NBOX = NCL
            IF(CMN.GE.ZMIN.AND.CMX.LE.ZMAX)   NBOX = NCL-1
            if(cmn.gt.zmin)  ISKP = 1
               DO I = 1,NCL
                  CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I)
                  CALL CPGETr ('clv - countour VALUE',clv)
                  IF(CLV.LT.ZMIN)  ISKP = ISKP+1
                  IF(CLV.LT.ZMIN.OR.CLV.GT.ZMAX)  THEN
                     NBOX = NBOX-1
                  END IF
!                  WRITE(6,*) I,CLV,ISKP,NBOX
               END DO
         END IF
         NBOX = MIN(NBOX,14)
!
!         write(6,*) CIS,ciu,NBOX,ncl,ZMIN,ZMAX

         DO I=1,NBOX+1

            IF(I.EQ.1)  THEN
               CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I)
               CALL CPGETr ('clv - countour VALUE',clv)
               CALL CPSETR ('ZDV - Z DATA VALUE',CLV+(ISKP-1)*CIU)
            ELSE IF(I.EQ.NBOX+1)  THEN
               CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I-2)
               CALL CPGETr ('clv - countour VALUE',clv)
               CALL CPSETR ('ZDV - Z DATA VALUE',CLV+(ISKP+1)*CIU)
            ELSE
               CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I-1)
               CALL CPGETr ('clv - countour VALUE',clv)
               CALL CPSETR ('ZDV - Z DATA VALUE',CLV+ISKP*CIU)
            END IF

!            CALL CPSETR ('ZDV - Z DATA VALUE',CLV)

!           if(cmx-cmn.gt.0.)   then
!              CALL CPSETR ('ZDV - Z DATA VALUE',cmn+(i-1)*ciu)
!           else
!              CALL CPSETR ('ZDV - Z DATA VALUE',clv+(i-nbox)*ciu)
!           end if

           CALL CPGETC ('ZDV - Z DATA VALUE',LLBS(I))
!           write(6,*) I,'  ',llbs(i)
         END DO
!
         CALL LBSETI ('CBL - COLOR OF BOX LINES',0)
         CALL LBSETI ('CLB - COLOR OF labels',2)
         CALL LBLBAR (0,.05,.95,.01,.08,NBOX,1.,.3,LIND(1+ISKP),0,LLBS &
     &               ,NBOX+1,1)
  100   continue

      ELSE

!        TURN OFF PLOTTING HIGHS AND LOWS

         CALL CPSETC('HLT - HIGH/LOW LABEL TEXT',' ')

!        DRAW NEGATIVE CONTOURS WITH DASHED LINES

         CALL CPPKCL(ARRAY, RWRK, IWRK)
         CALL CPGETI ('NCL - NUMBER OF CONTOUR LEVELS',NCLV)
         DO  ICLV=1,NCLV
       CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',ICLV)
       CALL CPGETR( 'CLV - CONTOUR LEVEL VALUE', CLV)
       IF(CLV .LT. 0.)THEN
          CALL CPSETC('CLD - CONTOUR LEVEL DASH PATTERN',  &
     &                ' ''''$$$''''$$$''''$$$ ')
       END IF
         END DO

!        DRAW CONTOUR LINES

!pl         CALL SETUSV('LW',3000)

         CALL CPCLDR(ARRAY,RWRK,IWRK)

!        DRAW LINE LABELS AND INFORMATIONAL LABEL

!pl         CALL SETUSV('LW',3000)

         CALL CPLBDR(ARRAY,RWRK,IWRK)
      END IF

!     DRAW PERIMETER

         ITICK1 = (AR-AL)/DX 
         ITICK2 = (AT-AB)/DZ 

!pl         CALL SETUSV('LW',2000)

         CALL GSPLCI(2)
         CALL PERIM(ITICK1,1,ITICK2,1)
 
!     CHECK AMOUNT OF WORKSPACE USED

      CALL CPGETI('IWU',IIWU)
      CALL CPGETI('RWU',IRWU)
!      WRITE(6,*) 'IWU = ',IIWU,'  RWU = ',IRWU
      IF(IRWU.GT.LRWK.OR.IIWU.GT.LIWK)  THEN
         WRITE(6,*) ' INSUFFICIENT WORKSPACE FOR CONPACK'
         WRITE(6,*) ' LIWK MUST BE AT LEAST ',IIWU
         WRITE(6,*) ' LRWK MUST BE AT LEAST ',IRWU
      END IF

      RETURN
      END

      SUBROUTINE COLRAM (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
!
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
!
! The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
! coordinates of points defining a polygon.  The area identifiers in
! the array IAIA, each with an associated group identifier in the array
! IGIA, tell us whether the polygon is to be color-filled or not.
!
!
! Assume the polygon will be filled until we find otherwise.
!
        IFLL=1
!
! If any of the area identifiers is negative, don't fill the polygon.
!
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
!
! Otherwise, fill the polygon in the color implied by its area
! identifier relative to edge group 3 (the contour-line group).
!
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
          IF (IFLL.GT.0.AND.IFLL.LT.15) THEN
            CALL GSFACI (IFLL+1)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
!
! Done.
!
        RETURN
!
      END
!
!
      SUBROUTINE CAPSAP (LABL,IAMA,LAMA)
!
        DIMENSION IAMA(*)
!
        CHARACTER*(*) LABL
!
! Compute and print the time required to draw the contour plot and how
! much space was used in the various arrays.
!
        PRINT * , 'PLOT TITLE WAS ',LABL
        CALL CPGETI ('IWU - INTEGER WORKSPACE USAGE',IIWU)
        CALL CPGETI ('RWU - REAL WORKSPACE USAGE',IRWU)
        PRINT * , 'INTEGER WORKSPACE USED ',IIWU
        PRINT * , '   REAL WORKSPACE USED ',IRWU
        IF (LAMA.NE.0) THEN
          IAMU=LAMA-(IAMA(6)-IAMA(5)-1)
          PRINT * , '   AREA MAP SPACE USED ',IAMU
        END IF
!
! Done.
!
        RETURN
!
      END

      SUBROUTINE LABTOP (LABL,SIZE)
!
        CHARACTER*(*) LABL
!
! Put a label just above the top of the plot.  The SET call is re-done
! to allow for the use of fractional coordinates, and the text extent
! capabilities of the package PLOTCHAR are used to determine the label
! position.
!
        CALL GETSET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG)
        SZFS=SIZE*(XVPR-XVPL)
        CALL    SET (0.,1.,0.,1.,0.,1.,0.,1.,1)
        CALL PCGETI ('QU - QUALITY FLAG',IQUA)
        CALL PCSETI ('QU - QUALITY FLAG',0)
        CALL PCSETI ('TE - TEXT EXTENT COMPUTATION FLAG',1)
        CALL PLCHHQ (.5,.5,LABL,SZFS,360.,0.)
        CALL PCGETR ('DB - DISTANCE TO BOTTOM OF STRING',DBOS)
        CALL PLCHHQ (.5*(XVPL+XVPR),YVPT+SZFS+DBOS,LABL,SZFS,0.,0.)
        CALL PCSETI ('QU - QUALITY FLAG',IQUA)
        CALL    SET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG)
!
! Done.
!
        RETURN
!
      END
!
      SUBROUTINE DFCLRS(IWKID)
!
! Define a set of RGB color triples for colors 0 through 15.
!
      PARAMETER (NCLRS=16)
      DIMENSION RGBV(3,NCLRS)
!
! Define the RGB color triples needed below.
!
      DATA RGBV /   &
     &     0.00 , 0.00 , 0.00 ,  &
     &     1.00 , 1.00 , 1.00 ,  &
     &     0.50 , 0.50 , 0.50 ,  &
     &     0.75 , 0.50 , 1.00 ,  &
     &     0.50 , 0.00 , 1.00 ,  &
     &     0.00 , 0.00 , 1.00 ,  &
     &     0.00 , 0.50 , 1.00 ,  &
     &     0.00 , 1.00 , 1.00 ,  &
     &     0.00 , 1.00 , 0.60 ,  &
     &     0.00 , 1.00 , 0.00 ,  &
     &     0.70 , 1.00 , 0.00 ,  &
     &     1.00 , 1.00 , 0.00 ,  &
     &     1.00 , 0.75 , 0.00 ,  &
     &     1.00 , 0.38 , 0.38 ,  &
     &     1.00 , 0.00 , 0.38 ,  &
     &     1.00 , 0.00 , 0.00 /
!
! Define 16 different color indices, for indices 0 through 15.  The
! color corresponding to index 0 is black and the color corresponding
! to index 1 is white.
!
      DO 101 I=1,NCLRS
         CALL GSCR (IWKID,I-1,RGBV(1,I),RGBV(2,I),RGBV(3,I))
 101  CONTINUE
!
! Done.
!
        RETURN
!
      END
!

