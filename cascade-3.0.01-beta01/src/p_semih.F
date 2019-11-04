      SUBROUTINE p_semih(X,WPART)
      IMPLICIT None
#include "casgki.inc"
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casshwr.inc"
#include "casmult.inc"
      COMMON/DIFFA/ PHI
      DOUBLE PRECISION ME,MP
      Double Precision XY
      COMMON/XVAR/ XY(10)
      DOUBLE PRECISION STHETA,SPHI
      DIMENSION X(20)
      Integer ndimc
      COMMON /DIMEN/ NDIMC
      Integer IMIX
      COMMON /OALPHAS/ IMIX
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      Integer NDIM
      COMMON/INTDIM/ NDIM
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      DOUBLE PRECISION QG21,QG22
      COMMON/SEMIH/ QG21,QG22
      Double Precision dcasrn
      EXTERNAL dcasrn
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision QF,W12
      DATA QF/0.0D0/,W12/0.0D0/

      Integer NDIMEN
      Integer I,J,IN,NP,NPP,KPA2,NRN,IST,NPFIN,KI,KPZ
      Integer nb2,kpf,kstart,kpao,kpao2,nacc
      Double Precision X,WPART,W02,QFT,YX,QT2
      Double Precision QG2MAX,QG2MIN,WQG22,SMIN
      Double Precision SHAT1,PHI,FGAM,FWEI,THETE,PR,CPHP,SPHP

      Double Precision STHP,CTHP
      Double precision xsi_hard,phi_k2,phi_k1,BOCHCK,ECM,WT,PT2,xp
      Double Precision wtglu,xpint,flux,wmat,XZOHRN
      Double Precision DOT,DOT1
	Double Precision XP1MIN,xp1,fwei1,wqg21
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar

	double precision xsim,ykt
	common/hzbval/xsim,ykt
	Double Precision sumkt,k1x,k1y,k2x,k2y,xmin1,ratio 

      double precision cut
      Double Precision ma,mb
      
      
      Double Precision sp,sq
      Integer KWZ,KWZP1,KWZP2
      Common/CASWZ/KWZ,KWZP1,KWZP2
      
      Integer Icoll
      Common/Colltest/Icoll
      
      Icoll = 0

      cut=pt2cut(ipro)
      cut = 0.25
cc      cut = 0.0
cc      cut = 1.
      NDIMEN = NDIM
      Do  IN=1,200
	 Do J=1,5
         K(IN,J) = 0
         P(IN,J) = 0
	 Enddo
      Enddo 
C...  GIVE BEAM  FOUR VECTORS
      DO 20 I=1,2
         DO 20 J=1,5
            K(I,J)=KBEAM(I,J)
   20 P(I,J) = PBEAM(I,J)
      If(Mult.ne.0) then
C... CALCULATE new CMS ENERGY
         DO 21 I=1,4
            CM(I) =P(1,I)+P(2,I)
   21    CONTINUE
      SSS = DOT(CM,CM)
      Endif
   
      ME =P(1,5)
      MP =P(2,5)
      NP=2
      NPP=3*NP
      N=2
      KT2_1=-99999d0
      KT2_2=-99999d0
      QF = 0.D0
      XG1SF=-9999d0
	XG2SF=-9999d0
      IF(IPRO.EQ.10) THEN
c         AM(1) = 0.140d0
         AM(1) = 0.250d0
         AM(1) = max(AM(1),0.0001d0)
         AM(2) = AM(1)
      ELSEIF(IPRO.EQ.11) THEN
         IF(IWEI.EQ.0) THEN
c            KPA = 4
	      KPA=IHFLA
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         AM(1) = PYMASS(KPA)
ccc         AM(2) = PYMASS(KPA)
         AM(2) = AM(1)
      ELSEIF(IPRO.EQ.15) THEN
         AM(1) = 0.250d0
         AM(1) = 0.0d0
c         AM(1) = max(AM(1),0.0001d0)
c         AM(2) = 2d0
      ELSEIF(IPRO.EQ.2) THEN
         AM(1) = PYMASS(Kpsi)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.3) THEN
         KPA=Kchi
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.4) THEN
         KPA=2
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.21) THEN
         KPA=Kpsi
         AM(1) = PYMASS(Kpsi)
         AM(2) = PYMASS(Kpsi)
      ELSEIF(IPRO.EQ.102) THEN
	   KPA = 25
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.500) THEN
         KPZ = 23
         KPA = KPZ
         AM(1) = PYMASS(KPA)
         AM(2) = 0.
      ELSEIF(IPRO.EQ.501) THEN
         KPZ = 23
         KPA = KPZ
         KPA2= 1
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.503) THEN
         NP=3
         NPP=3*NP
	   KPA = 1
         KPZ = 23
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.504) THEN
         NP=3
         NPP=3*NP
	   KPA = IHFLA
         KPZ = 23
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.510) THEN
      ELSEIF(IPRO.EQ.511) THEN
         KPZ = KWZ
         KPA = KPZ
         KPA2= 1
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.513) THEN
         NP=3
         NPP=3*NP
c        select charge of W
         KPZ = KWZ
	   KPA = KWZP1
	   KPA2 = KWZP2
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA2)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.514) THEN
         NP=3
         NPP=3*NP
c        select charge of W
         KPZ = KWZ
c         if(IHFLA.eq.5) then
c            write(6,*) ' p_semih: cabbibo supressed combination selected'
c            write(6,*) ' p_semih: heaviest flavor = ',IHFLA
c            write(6,*) ' p_semih: not implemented - Program stopped '
c            stop
c            endif
	   KPA = KWZP1
         KPA2 = KWZP2
c         write(6,*) ' KPZ, KPA KPA2 ',KPZ,KPA,KPA2
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA2)
         AM(3) = PYMASS(KPZ)
      ELSE
         WRITE(6,*) ' p_semih: wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
	Fwei1 = 1d0
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
      W02=(AM(1)+AM(2)+MP)**2
c      IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
      IF(AM(1).LT.1.0D0.and.pt2cut(ipro).ne.0.) then
           W02=(DSQRT(PT2CUT(IPRO)) + MP)**2
           elseif(AM(1).LT.1.0D0) then
           W02=(2*DSQRT(cut) + MP)**2
      Endif
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
c         WRITE(6,10500) YMIN,YMAX
c10500 FORMAT(' limits on y ',/,
c     +' YMIN = ',E10.5,' YMAX = ',E10.5)

C ... select particle code for light flavour production according
C ... to charges
      IF(IPRO.EQ.10) THEN
         IF(INTER.LT.2) THEN
            QF=DFLOAT(PYCHGE(1))**2 + DFLOAT(PYCHGE(2))**2 +
     +      DFLOAT(PYCHGE(3))**2
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) QF=3d0
            QF = 2.D0*QF
            KPA = -4
            QFT = - dcasrn()*QF
   30       KPA=KPA+1
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) Then
               QFT = QFT + 1d0
            else
               QFT = QFT + DFLOAT(PYCHGE(KPA))**2
            endif
            IF(QFT.LT.0.0D0.OR.KPA.EQ.0) GOTO 30
            IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!! ',
     +      KPA
            IF(KPA.EQ.0) write(6,*) 'fatal no quark !!!!!! ',
     +      KPA,QFT
            KPA2 = - KPA
         ENDIF
      ELSEIF(IPRO.EQ.11) THEN
         IF(IWEI.EQ.0) THEN

         ELSEIF(IWEI.EQ.1) THEN
            IF(dcasrn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
         AM(1)=PYMASS(KPA)

         AM(2)=AM(1)
      ELSEIF(IPRO.EQ.21) THEN
         KPA2=KPA
         AM(1)=PYMASS(KPA)
         AM(2)=AM(1)
      ELSEIF(IPRO.eq.501) THEN
      ELSEIF(IPRO.eq.503.or.IPRO.eq.504) THEN
         IF(IWEI.EQ.0) THEN

         ELSEIF(IWEI.EQ.1) THEN
            IF(dcasrn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
         AM(2)=PYMASS(KPA2)
         AM(1)=PYMASS(KPA)
      ELSEIF(IPRO.eq.513.or.IPRO.eq.514) THEN
         AM(2)=PYMASS(KPA2)
         AM(1)=PYMASS(KPA)
      ENDIF

C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... xg2 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,xg2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NRN = 0
      XMAX=0.999d0

      IF(Iabs(KE).eq.11.or.KE.eq.22) THEN
c.......
         IF(KE.NE.22) THEN
            NRN = NRN + 1

            YX = YMIN*((YMAX/YMIN)**X(NRN))
            IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 90
            Q2MIN=ME*ME*YX*YX/(1.D0-YX)
            IF(QMI.GT.Q2MIN) Q2MIN = QMI
            Q2MAX=YX*SSS - W12
            IF(QMA.LT.Q2MAX) Q2MAX = QMA
            IF(Q2MAX.LT.Q2MIN) GOTO 100
            NRN = NRN + 1
            Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
            XMIN = 0.0D0
            QG22 = 0.D0
            WQG22 = 1.D0
            IF(ISEMIH.EQ.1) THEN
               QT2 = Q2*(1.d0-yx)
               QG2MAX = SSS*YX - Q2 - (AM(1) + AM(2))**2

               QG2MIN = 0.0001d0
               QG2MIN = 0.000001d0

               QG2MAX=DMAX1(QG2MIN,QG2MAX)
               NRN = NRN + 1
               XZOHRN = X(NRN)
c 1/qt**2
               QG22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
               WQG22 = QG22*DLOG(QG2MAX/QG2MIN)
c... end semihard
            ENDIF

            SMIN = DMAX1(4*cut,2.D0*(AM(1)**2+AM(2)**2))
            XMIN=(SMIN)/(YX*SSS)

            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XMIN = DMAX1(XMIN,Q2/YX/SSS)
            ENDIF
         ELSEIF(KE.EQ.22) THEN
            Q2 = 0.D0
            YX = 1.D0
            QG22 = 0.D0
            WQG22 = 1.D0
            IF(ISEMIH.EQ.1) THEN
C... here is new for semihard
               QG2MAX = SSS - (AM(1)+AM(2))**2

               QG2MIN = 0.0001d0
               NRN = NRN + 1
               XZOHRN = X(NRN)
c 1/qt**2
               QG22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
               WQG22 = QG22*DLOG(QG2MAX/QG2MIN)

c... end new semihard
            ENDIF
c               write(6,*) ' p_semih yx,q2 ',yx,q2

            SMIN = DMAX1(4*cut,2.D0*(AM(1)**2+AM(2)**2))
            XMIN= (SMIN)/(YX*SSS)

         ELSE
            write(6,*) ' wrong KF selected; program stop '
            YX = 0.D0
         ENDIF
         XG1=YX
	   XG1SF = XG1
	   yx1 = yx
	   IF(XMIN.GE.XMAX) GOTO 110
C ... xg2 = E_parton /E_proton
         IF(IMIX.EQ.0) THEN
           NRN = NRN + 1
           xg2= XMIN*((XMAX/XMIN)**X(NRN))
         ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
           xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
         ELSE
           xg2 = 0.D0
         ENDIF
         XG2SF = XG2
         IF(IRES(1).EQ.1) THEN
	      If(Iabs(KE).EQ.11) then
               NRN = NRN + 1

               XP1MIN =DMAX1(4*cut,2.D0*(AM(1)**2+AM(2)**2))
     +           /yx/sss
               xp1 = XP1MIN*(1.D0/XP1MIN)**X(NRN)
               xg1 = xp1*yx
		   XG1SF = xp1  
               QG2MAX = sss*xg1 - Q2 -Qg22 - 4.D0 * AM(1)**2

               QG2MIN = 0.0001d0
               QG2MAX=DMAX1(QG2MIN,QG2MAX)
               NRN = NRN + 1
               QG21 = QG2MIN*((QG2MAX/QG2MIN)**X(NRN))
               WQG21 = QG21*DLOG(QG2MAX/QG2MIN)
	         Fwei1 = WQG21 * DLOG(1.D0/XP1MIN)
               IF(yx1.ne.1d0) THEN
                  xg1_bar = (-qg21)/(1.d0 - xg1)/sss
                  xg1_bar=0d0
                ELSE
                  xg1_bar=0d0
               ENDIF
		ELSEIF(KE.EQ.22) THEN
		
		ENDIF
         ENDIF
	ELSEIF(Iabs(KE).eq.2212) then
c only p-p included
         NRN = NRN + 1
         YX = YMIN*((YMAX/YMIN)**X(NRN))
         XG1=YX
	   XG1SF = XG1
         QG2MIN = 0.0001d0
         Q2MIN = QG2MIN
         Q2MAX=YX*(SSS-MP**2) - W12
                  
         IF(Q2MAX.LT.Q2MIN) GOTO 100
         NRN = NRN + 1
         Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         XMIN = 0.0D0
         XMIN = 0.0D0
         QG22 = 0.D0
         WQG22 = 1.D0
         QT2 = Q2 
c         QG21 = QT2
         QG2MAX = SSS*YX - Q2 - (AM(1)+AM(2))**2
                  
         QG2MAX=DMAX1(QG2MIN,QG2MAX)
         NRN = NRN + 1
         XZOHRN = X(NRN)
c 1/qt**2
         QG22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
         WQG22 = QG22*DLOG(QG2MAX/QG2MIN)
c        write(6,*) ' p_semih Q2,QG2: ',Q2,QG22
         IF(IPRO.EQ.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
	     Smin = AM(1)**2
           NRN = NRN + 1
           PHI_K1 = 2*PI*X(NRN)
           k1x = DSQRT(Q2)*SIN(PHI_K1)
           k1y = DSQRT(Q2)*COS(PHI_K1)
           NRN = NRN + 1
           PHI_K2 = 2*PI*X(NRN)
           k2x = DSQRT(QG22)*SIN(PHI_K2)
           k2y = DSQRT(QG22)*COS(PHI_K2)
	     
	     sumkt = (k1x+k2x)**2 + (k1y+k2y)**2
	     xmin = (SMIN+sumkt)/(YX*SSS)
           
           xg1 = yx
c           QT2 = Q2*(1.d0-yx)
           XG1_BAR = (-QT2)/(1.d0 - XG1)/SSS
	     sp = (-sss*xg1*xg1_bar - QG22 + xg1*sss + Smin +sumkt)
           sp = sp/(-xg1*sss)
           sq = sss*xg1*xg1_bar - QG22*xg1_bar - Smin -sumkt
           sq = sq/(-xg1*sss)
           xmin = -sp/2.- sqrt(sp**2/4. - sq)
           if(xmin.ne.xmin) then
            goto 120
           write(6,*) ' sp, sq ',sp,sq
           write(6,*) ' xmin1/2=',xmin,-sp/2.- sqrt(sp**2/4. - sq)
           write(6,*) ' xmin old ',xmin
         write(6,*) ' test xg1 ',xg1,xg1_bar
         write(6,*) ' test xg2 ',xmin,qg22/(1-xmin)/sss
         write(6,*) ' test qg21,qg22 ',qT2,qg22
           endif
	     xg2 = xmin 
c           write(6,*) 'masses',
c     & sqrt(SSS*(XG1*xg2+xg1*xg1_bar+xg1_bar*xg2_bar+xg2*xg2_bar)-sumkt)
c     &       ,am(1)
         Elseif(IPRO.eq.501.or.ipro.eq.511) then

           SMIN = DMAX1(4*cut,(AM(1)+AM(2))**2)
           
           XMIN=(SMIN+QT2+QG22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
           IF(IMIX.EQ.0) THEN
              NRN = NRN + 1
              xg2= XMIN*((XMAX/XMIN)**X(NRN))
           ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
              xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
           ELSE
             xg2 = 0.D0
           ENDIF
         Elseif(IPRO.eq.503.or.Ipro.eq.504.or.ipro.eq.513.or.ipro.eq.514) then

           SMIN = DMAX1(4*cut,(AM(1)+AM(2)+AM(3))**2)
           
           XMIN=(SMIN+QT2+QG22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
           IF(IMIX.EQ.0) THEN
              NRN = NRN + 1
              xg2= XMIN*((XMAX/XMIN)**X(NRN))
           ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
              xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
           ELSE
             xg2 = 0.D0
           ENDIF
	   Else

           SMIN = DMAX1(4*cut,(AM(1)+AM(2))**2)
           XMIN=(SMIN+QT2+QG22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
           IF(IMIX.EQ.0) THEN
              NRN = NRN + 1
              xg2= XMIN*((XMAX/XMIN)**X(NRN))
           ELSEIF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
              xg2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
           ELSE
             xg2 = 0.D0
           ENDIF
	   Endif
	   XG2SF = XG2
      Else
	   Write(6,*) ' IRES(1)=1 for KE = ',KE,' not implemented '
         Stop
	Endif


      IF(xg2.GE.0.9999.OR.xg2.LT.XMIN) GOTO 120
      SHAT1 = SSS*XG1*xg2
      IF(SHAT1.LT.(AM(1)+AM(2))**2) GOTO 130

      IST = 0
      YY=YX
      XG=xg2
      PHI = 99999.D0
c	If(ipro.ne.102) then 
c         NRN=NRN+1
c         PHI = 2.D0*PI*X(NRN)
c	   else
c	   phi=phi_k1
c	   endif
	If(ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then 
	   phi=phi_k1
	   else
         NRN=NRN+1
         PHI = 2.D0*PI*X(NRN)
	   endif
      CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      XP=0.D0
      IRES(2)=1
C here do some gymnastics to make cut on theta angle of electron
      IF(THEMA.LT.180.D0.OR.THEMI.GE.0.D0) THEN
C go first to ep LAB system
         CALL PYROBO(3,3,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
         THETE = PYP(3,14)
         CALL PYROBO(3,3,0.D0,0.D0,-CM(1)/CM(4),
     +                -CM(2)/CM(4),-CM(3)/CM(4))
c         CALL PYLIST(1)
         IF(THETE .GT. 180.01D0) THEN
            WRITE(6,*) ' FATAL: theta_electron > 180 deg ',THETE
         ENDIF
         IF(THETE .GT. THEMA) GOTO 140
         IF(THETE .LT. THEMI) GOTO 140
      ENDIF
C end of these gymnastics

C FINAL STATE PROTON

      NPFIN=NIA1+2+NP
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      K(NPFIN,3)=2

      PR = P(2,3)*(1.D0 - xg2)
      P(NPFIN,4) = PR
      CPHP=1.D0
      SPHP=DSQRT(1.D0 - CPHP**2)
      IF(ISEMIH.EQ.1) THEN
         STHP = DSQRT(QG22)/PR
      ENDIF
      IF(DABS(STHP).GT.1.D0) goto 190
      CTHP=DSQRT(1.D0 - STHP**2)
      P(NPFIN,1)= PR*STHP*CPHP
      P(NPFIN,2)= PR*STHP*SPHP
      P(NPFIN,3)= PR*CTHP
c     write(6,*) ' lin 372 '
c      CALL PYLIST(1)
C MOMENTA OF PARTON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KGL
      IF(IPRO.EQ.15) K(NIA1+1,2)=2
      K(NIA1+1,3)=2
      DO 40  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   40 CONTINUE
      IF(ISEMIH.EQ.0) THEN
         P(NIA1+1,4)=ABS(P(NIA1+1,3))
         P(NIA1+1,5)=0.0D0
      ELSE
         P(NIA1+1,4)=DSQRT(P(NIA1+1,1)**2+P(NIA1+1,2)**2+P(NIA1+1,3)**2
     +               - QG22)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
c here define new momenta: introduce maximum virtuality mu**2 = Q2 + 4M**2
c for the matrix element and kinematics we set Virt_max=kt**2
c     define lightcone variable xg=xg2
c         XG = xg2
c     define lightcone variable xg_bar=(-kt**2)/(1-xg)/SSS
         XG2_bar = (-QG22)/(1.d0-xg2)/SSS
c here check on virt. gamma
         QT2 = PYP(NIA1,9)
         IF(XG1.NE.1D0) THEN
            XG1_BAR = (-QT2)/(1.d0 - XG1)/SSS
         ELSE
            XG1_BAR=0D0
         ENDIF
         
         If(Icoll.eq.1) then      
            xg2_bar = 0. ! hannes hack
            xg1_bar = 0. ! hannes hack
         Endif
         
         
         KT2_1 = ABS(QT2)
                  
         xsi_hard=(xg1+xg2_bar)/(xg2+xg1_bar)

	   xsim = xsi_hard
         ykt = sqrt(abs(QG22))
         YY=YX
         if(xsi_hard.le.0) goto 200
         P(NIA1+1,4) =  xg2_bar*P(1,4) + xg2*P(2,4)
         P(NIA1+1,3) =  xg2_bar*P(1,3) + xg2*P(2,3)
	   if(ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
            else 
            NRN = NRN + 1
            PHI_K2 = 2*PI*X(NRN)
	   Endif
         P(NIA1+1,1) = DSQRT(QG22)*SIN(PHI_K2)
         P(NIA1+1,2) = DSQRT(QG22)*COS(PHI_K2)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
         KT2_2=ABS(QG22)
c	   write(6,*) KT2_2,QG22
         P(NIA1,4) =  xg1*P(1,4) + xg1_bar*P(2,4)
         P(NIA1,3) =  xg1*P(1,3) + xg1_bar*P(2,3)
c         write(6,*) ' xg1 ',xg1,xg1_bar
c         write(6,*) ' xg2 ',xg2,xg2_bar
c         write(6,*) ' qg21,qg22 ',QT2,qg22
         P(NIA1,5)=-sqrt(ABS(DOT1(NIA1,NIA1)))
         IF(IRES(1).EQ.1.AND.IABS(KE).NE.2212) Then
            NRN = NRN + 1
            PHI_K1 = 2*PI*X(NRN)
            P(NIA1,1) = DSQRT(QG21)*SIN(PHI_K1)
            P(NIA1,2) = DSQRT(QG21)*COS(PHI_K1)
		DO I=1,4
		P(NIA1-1,I) = P(4,I) - P(NIA1,I)
		ENDDO
		YY = XG1/YX
	   ENDIF
	   
         do i=1,5
            p(n+1,i) = P(NIA1,i)
            p(n+2,i) = P(NIA1+1,i)
         enddo
         P(N+2,3) =  0.5d0*(P(N+2,3)+P(N+2,4))
         P(N+2,4) =  DABS(P(N+2,3))
         P(N+1,3) =  0.5d0*(P(N+1,3)-P(N+1,4))
         P(N+1,4) =  DABS(P(N+1,3))
         K(N+1,1)=13
         K(N+2,1)=13
         N=N+2
c        write(6,*) ' p_semih xg: ',xg

         DO 50   KI=1,4
            P(NPFIN,KI)=P(2,KI)-P(NIA1+1,KI)
   50    CONTINUE
         P(NPFIN,5) = -sqrt(abs(dot1(npfin,npfin)))
      ENDIF
      NIA2 = NIA1+1
c      call PYLIST(1)
      NF1=NIA1+2
      NF2=NIA1+3
      if(IPRO.eq.503.or.Ipro.eq.504.or.Ipro.eq.513.or.IPRO.eq.514) NF2=NIA1+4
      KPF=KPA
      K(NF1,1)=2
      K(NF1,2)=KPF
      K(NF1,3)=NIA2
      P(NF1,5)=AM(1)
      K(NF2,1)=1
      K(NF2,2)=KPA2
      K(NF2,3)=NIA2
      P(NF2,5)=AM(2)
      
      if(IPRO.eq.503.or.Ipro.eq.504.or.IPRO.eq.513.or.Ipro.eq.514) then
        K(NF1+1,1)=1
        K(NF1+1,2)=KPZ
        K(NF1+1,3)=NIA2
        P(NF1+1,5)=AM(3)
        elseif(ipro.eq.15) then
        K(NF1,2)=21
        K(NF2,2)=2 
        elseif(ipro.eq.501.or.ipro.eq.511) then
        K(NF1,2)=KPA
        K(NF2,2)=KPA2 
        
      endif   
      if(kpa.eq.0) then
         write(6,*) 'p-semih kpa ',kpa
      endif
      NB2 = NIA2
C...   VECTOR OF GAMMA GLUON CMS SYSTEM
      DBCMS(1)=  P(NIA1,1) + P(NB2,1)
      DBCMS(2)=  P(NIA1,2) + P(NB2,2)
      DBCMS(3)=  P(NIA1,3) + P(NB2,3)
      DBCMS(4)=  P(NIA1,4) + P(NB2,4)
      DO 60  I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
         if(IPRO.eq.503.or.Ipro.eq.504.or.IPRO.eq.513.or.Ipro.eq.514) P(NF1+1,I)=0.0D0
   60 CONTINUE
      SHAT=DOT(DBCMS,DBCMS)
c      Write(6,*) ' SHAT ',SHAT,(DBCMS(I),I=1,4)
c      if(shat.ne.shat) stop
      IF(SHAT.LE.0.0) THEN
         GOTO 160
      ENDIF
      IF(IPRO.EQ.10) THEN
         IF(INTER.LT.2) THEN
            KSTART = -4
            QF=DFLOAT(PYCHGE(1))**2 + DFLOAT(PYCHGE(2))**2 +
     +      DFLOAT(PYCHGE(3))**2
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) QF=3d0
            QF = 2.D0*QF
            KPA = KSTART
            QFT = - dcasrn()*QF
   70       KPA=KPA+1
            IF(IABS(K(1,2)).eq.2212.and.IABS(K(2,2)).eq.2212) then
               QFT = QFT + 1.d0
            else
               QFT = QFT + DFLOAT(PYCHGE(KPA))**2
            endif
            IF(QFT.LT.0.0D0.OR.KPA.EQ.0) GOTO 70
            IF(KPA.GT.IABS(KSTART)-1) THEN
               write(6,*) 'fatal light quark = charm!!!!!! ',KPA
            ENDIF
            IF(KPA.EQ.0) write(6,*) 'fatal no quark !!!!!! ',
     +      KPA,QFT
            KPA2 = - KPA
            IF(IGENFL.EQ.0) THEN
               KPAO = KPA
               KPAO2 = KPA2
            ELSE
               KPA=KPAO
               KPA2 = KPAO2
            ENDIF
            K(NF1,2)=KPA
            K(NF2,2)=KPA2
         ENDIF
      ENDIF

      IF(IPRO.EQ.2) THEN
         K(NF1,2)=KPSI
         K(NF2,2)=21
      ENDIF
c      call PYLIST(1)
C NOW BOOST TO GAMMA GLUON
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2 +
     +(DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 150
      CALL PYROBO(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),
     +-DBCMS(2)/DBCMS(4),-DBCMS(3)/DBCMS(4))
      SPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c      write(6,*) ' after boost '
c      call PYLIST(1)
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
c      if(ECM.LE.0)  write(6,*) ECM
      ECM =DSQRT(SHAT)
c      if(ECM.ne.ECM)  write(6,*) ECM,SHAT

      If(ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500)  then 
      else
c cut only for a 2-> 2 process, if resonance, EMC = AM(1) by construction
c if cutted, only precision problems are encountered
         if(IPRO.eq.503.or.Ipro.eq.504.or.Ipro.eq.513.or.IPRO.eq.514) then
            IF(ECM.LE.(AM(1)+AM(2)+AM(3))) GOTO 170
         else
            IF(ECM.LE.(AM(1)+AM(2))) GOTO 170
         endif

	endif 
      IF(IPRO.EQ.3.or.IPRO.eq.4.or.IPRO.EQ.102.or.ipro.eq.500) Then
         K(NF2,2)=0
         K(NF2,1)=-1
         K(NF2,3)=0
         P(NF2,5)=0
         K(NF1,1)=1
         K(NF1,3)=4
        DO I=1,4
            P(NF1,I)=0
            P(NF2,I)=0
         Enddo
         P(NF1,4)=ECM
         P(NF1,5)=AM(1)
         if(ecm/am(1).le.0.98.or.ecm/am(1).gt.1.02) then
c            write(6,*) ' p_semih ecm,Am ',ecm,am(1),sqrt(sumkt)
         endif
c        WT to corrrect for it afterwards	   
	   WT = (2.D0*PI)**NPP
c	   call pylist(1)
	ELSE
c         write(6,*) ' IPRO = 13 ,NRN= ',NRN
         IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
            XY(1)=X(NDIMEN+3)
            XY(2)=X(NDIMEN+4)
            if(IPRO.eq.503.or.Ipro.eq.504.or.Ipro.eq.513.or.IPRO.eq.514) then
              XY(3)=X(NDIMEN+5)
              XY(4)=X(NDIMEN+6)
              XY(5)=X(NDIMEN+7)
           endif
         ELSE
            NRN = NRN + 1
            XY(1) = X(NRN)
            NRN = NRN + 1
            XY(2) = X(NRN)
            if(IPRO.eq.503.or.Ipro.eq.504.or.Ipro.eq.513.or.IPRO.eq.514) then
              NRN = NRN + 1
              XY(3) = X(NRN)
              NRN=NRN+1
              XY(4)=X(NRN)
              NRN=NRN+1
              XY(5)=X(NRN)
           endif
         ENDIF
c         write(6,*) am(1),am(2),ipro
c         write(6,*) 'p_semih ',XY(1),XY(2),NDIM,NDIMEN,IMIX,IWEI
         CALL PHASE(NP,ECM,AM,PCM,WT)
         IF(WT.LE.0.D0) GOTO 180
         DO 80  I=1,4
            P(NF1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,2)
            if(IPRO.eq.503.or.Ipro.eq.504.or.Ipro.eq.513.or.IPRO.eq.514) then
               P(NF1+1,I)=PCM(I,3)
            endif
   80    CONTINUE 
c         write(6,*) ' p_semih ', KPA,KPA2,K(NF1,2),K(NF2,2)
	Endif

C BOOST BACK TO OVERALL CMS 
      call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL PYROBO(0,N,0.D0,0.D0, DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +DBCMS(3)/DBCMS(4))
cc remove cuts to meoffsh routines

      XG = xg2
      if(xg2.GT.1.d0) write(6,*) 'p_semih xg2>1',xg2,ipro
c      call pylist(1)
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
      XPINT=DLOG(XMAX/XMIN)
      WTGLU = XPINT
	if(ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) wtglu=1
      WTGLU = WTGLU*WQG22
      WMAT=WT *WTGLU*FGAM*FWEI * Fwei1
c	write(6,*) FGAM,FWEI,FWEI1
      WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
      FLUX=(2.D0*PI)**4
      FLUX=FLUX/(2.D0*xg2*xg1*sss)          

c      ma = sqrt(qt2)
c      mb = sqrt(qg22)
c      mb=0
c      ma=mb
c      if((shat -(ma+mb)**2).gt.0) then          
c         FLUX=(2.D0*PI)**4/
c     &      2.D0/sqrt((shat -(ma+mb)**2)*
c     &                (shat -(ma-mb)**2))
ccc
c         write(6,*) flux,qt2,qg22
c          else 
c          flux=0.
c      endif
ccc      FLUX=(2.D0*PI)**4/(2.D0*xg2*xg1*sss)          

      WPART = WMAT *FLUX
	IF(IPRO.EQ.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
	  Wpart = WTGLU*FGAM*FWEI * Fwei1
        
        if(Ipro.eq.3) then
          Wpart=Wpart*Flux/((2.D0*PI)**3)
c 
c xsection for resonance production gg-> chi_c
c sig = \int dx1 dx2 g(x1) g(x2) sigma_hat \delta(x1x2 * sss - M2)
c     = \int dx1 dx2 g(x1) g(x2) sigma_hat 1/(x1*sss) \delta(x2 - M2/(x1*sss))
c     = \int dx1 g(x1) g(x2=M2/x1*sss) sigma_hat x2/M2 
c     = \int dx1/x1  x1g(x1) x2g(x2) sigma_hat/M2
          Wpart=Wpart/AM(1)**2
        endif
	Endif
      PT2H = PT2
      SHH = SHAT
      NDIMC = NRN
      CALL PYEDIT(14)
c      write(6,*)  'p_semih: WPART ',WPART
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSE
         KPA = KPAO
         KPA2 = KPAO2
      ENDIF
	
c      write(6,*) ' end of p_semih ',NDIMC,NRN
c      CALL PYLIST(1)
	
	
      if(wpart.ne.wpart) write(6,*) ' p_semih wpart = ',wpart,xsi_hard
c      write(6,*) ' p_semih end ,kpa,kpa1,kpa2,kpao2',kpa,kpa1,kpa2,kpao2
      RETURN
   90 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      RETURN
  100 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
         NDIMC = 9999
      ENDIF
      RETURN
  110 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: xlimit ; RETURN xmin,xmax',xmin,xmax
      ENDIF
      NDIMC = 9999
      RETURN
  120 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: xg2imit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  130 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: slimit ; RETURN ',shat1,am(1),am(2),IPRO
         write(6,*) ' sss,xg1,xg2,q2,qg22 ',sss,xg1,xg2,q2,qg22
         write(6,*) ' p_semih q2,y,sss',q2,yx,sss
         write(6,*) ' p_semih 2*pq ',2*dot1(nia1,2),sss*yx
      ENDIF
      NDIMC = 9999
      RETURN
  140 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: theta limit ; RETURN ',themi,thema,thete
      ENDIF
      NDIMC = 9999
      RETURN
  150 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: bochck limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  160 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: shat ; RETURN ',shat,shat1,kpa
c        write(6,*) ' p_semih: QG2MAX,QG22 ',QG2MAX,QG22
c        write(6,*) ' p_semih: yx,xg2,q2,s ',yx,xg2,q2,sss
c        write(6,*) ' p_semih 2*pq ',2*dot1(nia1,2),sss*yx
c        write(6,*) ' p_semih nia1,nia2 ',nia1,nia2
c        write(6,*) ' p_semih 2*ktgam*ktglu ',2*(P(Nia1,1)*P(nia2,1)+
c     &  P(Nia1,2)*P(NIA2,2))
      ENDIF
      NDIMC = 9999
      RETURN
  170 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: ECM limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  180 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: PHASE WT=0 ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  191 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: PTCUT limit ; RETURN ',PT2,KPA
      ENDIF
      NDIMC = 9999
      RETURN

  190 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: QG22 limit ; RETURN ',DSQRT(QG22),PR,STHP
      ENDIF
      NDIMC = 9999
      RETURN
  200 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_semih: xsi_hard limit ; RETURN ',xsi_hard
      ENDIF
      NDIMC = 9999
      RETURN


      END
