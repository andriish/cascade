      SUBROUTINE p_gamgam(X,WPART)
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
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar
	Integer Iel1,Iel2,Nph1,Nph2
	Common/cagaminfo/Iel1,Iel2,Nph1,Nph2
	Double Precision phi_1,phi_2 
      Double Precision QF,W12
      DATA QF/0.0D0/,W12/0.0D0/
      Integer NDIMEN
      Integer I,J,IN,NP,NPP,KPA2,NRN,IST
      Integer kpf,kstart,kpao,kpao2,nacc
      Double Precision X,WPART,W02,QFT
	Double Precision Q21
	Double Precision Q22 
	Double precision y2min
      Double Precision SMIN
      Double Precision SHAT1,PHI,THETE 
      Double Precision FGAM1,FWEI1,FGAM2,FWEI2
      Double precision xsi_hard,BOCHCK,ECM,WT,PT2
      Double Precision flux,wmat
      Double Precision DOT,DOT1
	Double Precision XP1MIN,XP1,WQG21
	Double Precision XP2MIN,XP2,WQG22
	Double Precision QG2MAX,QG2MIN
	Double Precision Alph_em
      NDIMEN = NDIM
      DO 10  IN=1,20
         K(IN,1) = 0
         K(IN,2) = 0
   10 CONTINUE

c      IDEBUG = 1

C...  GIVE BEAM  FOUR VECTORS
      DO 20 I=1,2
         DO 20 J=1,5
            K(I,J)=KBEAM(I,J)
   20 P(I,J) = PBEAM(I,J)
C... CALCULATE CMS ENERGY
      DO 21 I=1,4
         CM(I) =P(1,I)+P(2,I)
   21 CONTINUE
      SSS = DOT(CM,CM)
      ME =P(1,5)
      MP =P(2,5)
      NP=2
      NPP=3*NP
      N=2
	Nia1 = 0
	Nia2 = 0
	Nir1 = 0
	nir2 = 0
	xmax = 0.9999
	xmin = 1.d-8
      KT2_1=-99999d0
      KT2_2=-99999d0
      QF = 0.D0


      IF(IPRO.EQ.10) THEN
         AM(1) = 0.0D0
         AM(2) = 0.0D0
         AM(1) = 0.140d0
c         AM(1) = max(AM(1),0.0001d0)
         AM(1) = 0.140d0
         AM(2) = AM(1)
         IF(ISEMIH.EQ.1) THEN
            AM(1) = 0.140d0
            AM(1) = max(AM(1),0.0001d0)
            AM(2) = AM(1)
         ENDIF
      ELSEIF(IPRO.EQ.11) THEN
         IF(IWEI.EQ.0) THEN
            KPA = 4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
      ELSEIF(IPRO.EQ.2) THEN
         AM(1) = PYMASS(Kpsi)
         AM(2) = 0.d0
      ELSE
         WRITE(6,*) 'wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
      W02=(AM(1)+AM(2)+MP)**2
      IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
      IF(AM(1).LT.1.0D0.and.pt2cut(ipro).ne.0.)
     +     W02=(DSQRT(PT2CUT(IPRO)) + MP)**2
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
            KPA=4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
            IF(dcasrn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
         AM(1)=PYMASS(KPA)
         AM(2)=PYMASS(KPA2)
      ENDIF

C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... XG1 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,XG2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NRN = 0
c.......
      ALPH_EM = ALPH
      IF(KE.NE.22) THEN
         NRN = NRN + 1
         YX1 = YMIN*((YMAX/YMIN)**X(NRN))
         IF(YX1.GT.YMAX.OR.YX1.LT.YMIN) GOTO 90
         Q2MIN=ME*ME*YX1*YX1/(1.D0-YX1)
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         Q2MAX=YX1*SSS - W12
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
c	   write(6,*) Qmi,Qma,q2min,q2max
         IF(Q2MAX.LT.Q2MIN) GOTO 100
         NRN = NRN + 1
         Q21 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q21)
         SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))	
	   
         Fgam1 = (1.D0-yx1+yx1**2/2.D0)/Q21/yx1 - ME*ME*YX1/Q21**2 
         Fgam1 = alph_em * Fgam1 / pi
         Fwei1 = yx1*Q21*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
c  
         IF(yx1.ne.1d0) THEN
c            yx1_bar = (-q21)/(1.d0 - yx1)/sss
            yx1_bar = (-q21)/sss
            yx1_bar=0d0
        ELSE
            yx1_bar=0d0
         ENDIF
         NRN = NRN + 1
         PHI_1 = 2*PI*X(NRN)
c virtual photon 
	   K(N+2,1) = 21
	   K(N+2,2) = KPH
	   K(N+2,3) = N+1
	   P(N+2,1) = dsqrt(q21*(1.d0 - yx1))*sin(phi_1)
	   P(N+2,2) = dsqrt(q21*(1.d0 - yx1))*cos(phi_1)
	   P(N+2,3) = yx1*P(1,3) + yx1_bar*P(2,3)
	   P(N+2,4) = yx1*dabs(P(1,4)) + yx1_bar*dabs(P(2,4))	   	   
         P(N+2,5)=-sqrt(dabs(dot1(N+2,N+2)))
	   Nph1 = N+2
c	   write(6,*) ' q2 ',P(N+2,5)**2,Q21 
C scattered electron	   
	   K(N+1,1) = 1
	   K(N+1,2) = KE
	   K(N+1,3) = 1
	   P(N+1,1) = P(1,1) - P(N+2,1)
	   P(N+1,2) = P(1,2) - P(N+2,2)
	   P(N+1,3) = P(1,3) - P(N+2,3)
	   P(N+1,4) = P(1,4) - P(N+2,4)
         P(N+1,5)= me
	   Iel1 = N+1
c	   write(6,*) ' electron ',P(N+1,5)**2  
	   
         N = N + 2
		
      ELSEIF(KE.EQ.22) THEN
         Q21 = 0.D0
         YX1 = 1.D0
C interacting photon  
	   K(N+1,1) = 21
	   K(N+1,2) = KPH
	   K(N+1,3) = N+1
	   P(N+1,1) = P(1,1)
	   P(N+1,2) = P(1,2)
	   P(N+1,3) = P(1,3)
	   P(N+1,4) = P(1,4)	   	   
         P(N+1,5) = P(1,5)
         N = N + 1
	   Nph1 = N 
	   Nia1 = N 
	   Fgam1 = 1
	   Fwei1 = 1
      ELSE
         write(6,*) ' wrong KF selected; program stop '
         YX1 = 0.D0
      ENDIF
      XG1=YX1
	NIA1 = N
c	write(6,*) ' 1st NIA1 ',nia1
      IF(IRES(1).EQ.1) THEN
	   If(KE.EQ.22) Then
c  gam gam  
            NRN = NRN + 1
            XG1 = YMIN*((YMAX/YMIN)**X(NRN))
            QG2MAX = sss*xg1 - 4.D0 * AM(1)**2
            IF(ICCFM.EQ.0) QG2MAX = DMIN1(4D0*AM(1)**2,QG2MAX)
            QG2MIN = 0.0001d0
            QG2MAX=DMAX1(QG2MIN,QG2MAX)
            NRN = NRN + 1
            QG21 = QG2MIN*((QG2MAX/QG2MIN)**X(NRN))
            WQG21 = QG21*DLOG(QG2MAX/QG2MIN)
c  factor of xg1 is already included since we use xg(x)
		Fwei1 = DLOG(YMAX/YMIN) * WQG21
            NRN = NRN + 1
            PHI_1 = 2*PI*X(NRN)
            xg1_bar=0d0
	   Elseif(Iabs(KE).EQ.11) Then
            NRN = NRN + 1
            XP1MIN =DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
     +           /yx1/sss
            xp1 = XP1MIN*(1.D0/XP1MIN)**X(NRN)
            xg1 = xp1*yx1         
		Fwei1 = Fwei1 * DLOG(1.D0/XP1MIN)
            QG2MAX = sss*xg1 - Q21 - 4.D0 * AM(1)**2
            IF(ICCFM.EQ.0) QG2MAX = DMIN1(Q21+4D0*AM(1)**2,QG2MAX)
            QG2MIN = 0.0001d0
            QG2MAX=DMAX1(QG2MIN,QG2MAX)
            NRN = NRN + 1
            QG21 = QG2MIN*((QG2MAX/QG2MIN)**X(NRN))
            WQG21 = QG21*DLOG(QG2MAX/QG2MIN)
	      Fwei1 = Fwei1 * WQG21
            IF(yx1.ne.1d0) THEN
c               xg1_bar = (-qg21)/(1.d0 - xg1)/sss
               xg1_bar = (-qg21)/sss
               xg1_bar=0d0
            ELSE
               xg1_bar=0d0
            ENDIF
	   Endif
	   
	   
	   N = N + 1 
	   K(N,1) = 21 
	   K(N,2) = 21
	   K(N,3) = N-1
	   P(N,1) = dsqrt(qg21*(1.d0 - xg1))*sin(phi_1)
	   P(N,2) = dsqrt(qg21*(1.d0 - xg1))*cos(phi_1)
	   P(N,3) = xg1*P(1,3) + xg1_bar*P(2,3)
	   P(N,4) = xg1*P(1,4) + xg1_bar*P(2,4)	   	   
         P(N,5)=-sqrt(ABS(DOT1(N,N)))
	   NIA1 = N	  
c	   write(6,*) ' kt1 : ',q21*(1.d0 - xg1),qg21,qg2min,qg2max
      ENDIF
	
c	write(6,*) 'after 1st step '
c	call pylist(1)
	
c      IF(XMIN.GE.XMAX) GOTO 110
C ... XG2 = E_parton /E_proton
      IF(KP.NE.22) THEN
         NRN = NRN + 1
	   Y2min = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
     &	   /xg1/sss
         YX2 = Y2MIN*((YMAX/Y2MIN)**X(NRN))
         IF(YX2.GT.YMAX.OR.YX2.LT.Y2MIN) GOTO 90
         Q2MIN=ME*ME*YX2*YX2/(1.D0-YX2)
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         Q2MAX=YX2*SSS - W12
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
         IF(Q2MAX.LT.Q2MIN) GOTO 100
         NRN = NRN + 1
         Q22 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q22)
         IF(yx2.ne.1d0) THEN
c            yx2_bar = (-q22)/(1.d0 - yx2)/sss
           yx2_bar = (-q22)/sss
            yx2_bar=0d0
         ELSE
            yx2_bar=0d0
         ENDIF
         Fgam2 = (1.D0-yx2+yx2**2/2.D0)/Q22/yx2 - ME*ME*YX2/Q22**2 
         Fgam2 = alph_em * Fgam2 / pi
         Fwei2 = yx2*Q22*DLOG(YMAX/Y2MIN)*DLOG(Q2MAX/Q2MIN)
c  
         NRN = NRN + 1
         PHI_2 = 2*PI*X(NRN)
c virtual photon 
	   K(N+2,1) = 21
	   K(N+2,2) = KPH
	   K(N+2,3) = N+1
	   P(N+2,1) = dsqrt(q22*(1.d0 - yx2))*sin(phi_2)
	   P(N+2,2) = dsqrt(q22*(1.d0 - yx2))*cos(phi_2)
	   P(N+2,3) = yx2*P(2,3) + yx2_bar*P(1,3)
	   P(N+2,4) = yx2*dabs(P(2,4)) + yx2_bar*dabs(P(1,4))	   	   
         P(N+2,5)=-sqrt(dabs(dot1(n+2,n+2)))
         Nph2 = N+2
c      P(N+2,5)=-sqrt(q22)
c      write(6,*) ' q2 ',P(N+2,5) ,Q22 
c      write(6,*) ' q2 ',P(N+2,5)**2,Q22 
c      write(6,*) ' 1st  ',yx1,yx1_bar,q21
c      write(6,*) ' 2nd  ',yx2,yx2_bar,q22
	   
C scattered electron	   
	   K(N+1,1) = 1
	   K(N+1,2) = KP
	   K(N+1,3) = 1
	   P(N+1,1) = P(2,1) - P(N+2,1)
	   P(N+1,2) = P(2,2) - P(N+2,2)
	   P(N+1,3) = P(2,3) - P(N+2,3)
	   P(N+1,4) = P(2,4) - P(N+2,4)
         P(N+1,5)= me
c	  write(6,*) ' electron ',P(N+1,5)**2  
	   Iel2 = N + 1
         N = N + 2
         NIA2 = N	  

      ELSEIF(KP.EQ.22) THEN
         Q22 = 0.D0
         YX2 = 1.D0
C interacting photon  
	   K(N+1,1) = 21
	   K(N+1,2) = KPH
	   K(N+1,3) = N+1
	   P(N+1,1) = P(2,1)
	   P(N+1,2) = P(2,2)
	   P(N+1,3) = P(2,3)
	   P(N+1,4) = P(2,4)	   	   
         P(N+1,5) = P(2,5)
	   N = N + 1
	   Nph2 = N 
         NIA2 = N	  
	   Fgam2 = 1
	   Fwei2 = 1
      ELSE
         write(6,*) ' wrong KF selected; program stop '
         YX2 = 0.D0
      ENDIF
      XG2=YX2
	NIA2=N
c	write(6,*) ' 2nd NIA2 ', NIA2
      IF(IRES(2).EQ.1) THEN
	   If(KP.EQ.22) Then
c  gam gam  
            XP2MIN =DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
     &           /xg1/sss
            NRN = NRN + 1
            XG2 = XP2MIN*((YMAX/XP2MIN)**X(NRN))
            QG2MAX = sss - 4.D0 * AM(1)**2
            IF(ICCFM.EQ.0) QG2MAX = DMIN1(4D0*AM(1)**2,QG2MAX)
            QG2MIN = 0.0001d0
            QG2MAX=DMAX1(QG2MIN,QG2MAX)
            NRN = NRN + 1
            QG22 = QG2MIN*((QG2MAX/QG2MIN)**X(NRN))
            WQG22 = QG22*DLOG(QG2MAX/QG2MIN)
c  factor of xg2 is already included since we use xg(x)
		Fwei2 = DLOG(YMAX/XP2MIN) * WQG22
            xg2_bar=0d0
            NRN = NRN + 1
            PHI_2 = 2*PI*X(NRN)
	   Elseif(Iabs(KP).EQ.11) Then
            NRN = NRN + 1
            XP2MIN =DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
     &           /xg1/yx2/sss
            XP2 = XP2MIN*(1.D0/XP2MIN)**X(NRN)
            xg2 = xp2*yx2
		Fwei2 = Fwei2 * DLOG(1.D0/XP2MIN)
            QG2MAX = sss*xg1*xg2 - Q21 -Q22 - 4.D0 * AM(1)**2
            IF(ICCFM.EQ.0) QG2MAX = DMIN1(Q22+4D0*AM(1)**2,QG2MAX)
            QG2MIN = 0.0001d0
            QG2MAX=DMAX1(QG2MIN,QG2MAX)
            NRN = NRN + 1
            QG22 = QG2MIN*((QG2MAX/QG2MIN)**X(NRN))
            WQG22 = QG22*DLOG(QG2MAX/QG2MIN)
	      Fwei2 = Fwei2 * WQG22
            IF(yx2.ne.1d0) THEN
c               xg2_bar = (-qg22)/(1.d0 - xg2)/sss
               xg2_bar = (-qg22)/sss
               xg2_bar=0d0
            ELSE
               xg2_bar=0d0
            ENDIF
	   Endif
	   N = N + 1 
	   K(N,1) = 21 
	   K(N,2) = 21
	   K(N,3) = N-1
	   P(N,1) = dsqrt(qg22*(1.d0 - xg2))*sin(phi_2)
	   P(N,2) = dsqrt(qg22*(1.d0 - xg2))*cos(phi_2)
	   P(N,3) = xg2*P(2,3) + xg2_bar*P(1,3)
	   P(N,4) = xg2*P(2,4) + xg2_bar*P(1,4)	   	   
         P(N,5)=-sqrt(ABS(DOT1(N,N)))
	   NIA2 = N	  
c	   write(6,*) ' kt2 : ',qg22*(1.d0 - xg2),qg22,qg2min,qg2max
      ENDIF
c      write(6,*) ' ires(1) ',ires(1),XG1,YX1,Q21
c      write(6,*) ' ires(2) ',ires(2),XG2,YX2,Q22
c      call pylist(1)
	
      IF(XG2.GE.1.0001d0.OR.XG2.LT.XMIN) GOTO 120
      SHAT1= SSS*XG1*XG2 - QG21 - QG22
      IF(ISEMIH.EQ.1)  then 
	   if(xg1.eq.1.or.xg2.eq.1) then
	   else
	     SHAT1=SSS*XG1*XG2 - Q21/(1-xg1) - QG22/(1-xg2)
	   endif
	endif
      IF(SHAT1.LT.(AM(1)+AM(2))**2) GOTO 130

      IST = 0
      YY=XG1
	XG1SF = XG1
	IF(IRES(1).EQ.1) XG1SF=XG1/YX1
      XG=XG2
	XG2SF = XG1
	IF(IRES(2).EQ.1) XG2SF=XG2/YX2
      KT2_1 = P(NIA1,1)**2 + P(NIA1,2)**2
      KT2_2 = P(NIA2,1)**2 + P(NIA2,2)**2
	
c      write(6,*) ' nia1,nia2 ',nia1,nia2,N
c      write(6,*) ' nir1,nir2 ',nir1,nir2
	
	Do I = 1,5
	   K(N + 1,I) = K(NIA1,I)
	   P(N + 1,I) = P(NIA1,I)
	   K(N + 2,I) = K(NIA2,I)
	   P(N + 2,I) = P(NIA2,I)
	Enddo
c	K(Nia1,1)=13
c	K(Nia2,1)=13	
	Nia1 = N+1
	Nia2 = N+2
	N = N + 3
c	write(6,*) ' bef N= ',N,' Nia1, Nia2 ',Nia1, Nia2
	IF(K(NIA1,2).EQ.21) THEN
	  Nir1 = N + 1
	  do I = 1,5
	  P(Nir1,I) = P(nph1,I) - P(NIA1,I) 
	  Enddo
	  K(NIR1,1)=2	  
	  K(NIR1,2)=21	  
	  K(NIR1,3)=nph1
	  N = N + 1	  
	ELSE
	  Nir1 = 0
	ENDIF
c	write(6,*) ' 1st N= ',N,' Nia1, Nia2 ',Nia1, Nia2
c	write(6,*) ' 1st N= ',N,' Nir1, Nir2 ',Nir1, Nir2
	
	IF(K(NIA2,2).EQ.21) THEN
	  Nir2 = N + 1
	  do I = 1,5
	  P(Nir2,I) = P(nph2,I) -  P(NIA2,I)
	  Enddo
	  K(NIR2,1)=2
	  K(NIR2,2)=21
	  K(NIR2,3)=nph2
	  N = N + 1
	ELSE
	  Nir2 = 0
	ENDIF
	N = N + 1
c	write(6,*) ' 2nd N= ',N,' Nia1, Nia2 ',Nia1, Nia2
c	write(6,*) ' 2nd N= ',N,' Nir1, Nir2 ',Nir1, Nir2
c      write(6,*) '1st nia2 ',(P(NIA2,i),i=1,4)
      NF1=NIA2+1
      NF2=N
      K(NF1,1)=2
c      call PYLIST(1)
c      pause
      K(NF2,2)=KPA2
      KPF=KPA
      if(kpa.eq.0) then
         write(6,*) 'p-semih kpa ',kpa
      endif
      K(NF1,2)=KPF
      K(NF1,3)=NIA1
      K(NF2,1)=1
      K(NF2,3)=NIA1
      P(NF1,5)=AM(1)
      P(NF2,5)=AM(2)
C...   VECTOR OF GAMMA GLUON CMS SYSTEM

      DBCMS(1)=  P(NIA1,1) + P(NIA2,1)
      DBCMS(2)=  P(NIA1,2) + P(NIA2,2)
      DBCMS(3)=  P(NIA1,3) + P(NIA2,3)
      DBCMS(4)=  P(NIA1,4) + P(NIA2,4)
      DO 60  I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
   60 CONTINUE
c      write(6,*) ' dbcms ', (DBCMS(i),i=1,4)
c      write(6,*) ' nia1 ',nia1, (P(nia1,i),i=1,4)
c      write(6,*) ' nia2 ',nia2, (P(nia2,i),i=1,4)
      SHAT=DOT(DBCMS,DBCMS)
c	write(6,*) ' SHAT ',SHAT
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
      CALL PYROBO(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +-DBCMS(3)/DBCMS(4))
      SPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c      write(6,*) ' after boost '
c      call PYLIST(1)
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
c         write(6,*) ECM
      ECM =DSQRT(SHAT)
C      IF(ECM.LE.(AM(1)+AM(2))) WRITE(6,*) ' ECM LE MASS',ECM
      IF(ECM.LE.(AM(1)+AM(2))) GOTO 170
      IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
         XY(1)=X(NDIMEN+3)
         XY(2)=X(NDIMEN+4)
      ELSE
         NRN = NRN + 1
         XY(1) = X(NRN)
         NRN = NRN + 1
         XY(2) = X(NRN)
      ENDIF
c         write(6,*) 'p_gamgam ',XY(1),XY(2),NDIM,NDIMEN,IMIX,IWEI
c      Write(6,*) ' NP,ECM,AM ', NP,ECM,AM(1),AM(1),SHAT
      CALL PHASE(NP,ECM,AM,PCM,WT)
      IF(WT.LE.0.D0) GOTO 180
      DO 80  I=1,4
         P(NF1,I)=PCM(I,1)
         P(NF2,I)=PCM(I,2)
   80 CONTINUE
c      call pylist(1)
c         write(6,*) ' p_gamgam ', KPA,KPA2,K(NF1,2),K(NF2,2)
      PT2 = PYP(NF1,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) GOTO 191

C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
      call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)

      if(xg2.GT.1.d0) write(6,*) 'p_gamgam xg2>1',xg2,ipro


c      call pylist(1)
c	pause
      WMAT=WT  * FGAM1 * FWEI1 * FGAM2 * FWEI2
      WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
      FLUX=(2.D0*PI)**4
      FLUX=FLUX/(2.D0*xg2*xg1*sss)
      WPART = WMAT *FLUX

c      write(6,*) ' Fwei ',FGAM1,FWEI1,FGAM2,FWEI2
c      write(6,*) ' wmat ',wmat,flux

      PT2H = PT2
      SHH = SHAT


      NDIMC = NRN
      CALL PYEDIT(14)
c      write(6,*) ' end of p_gamgam ',NDIMC,NRN
c      CALL PYLIST(1)
c      write(6,*)  'p_gamgam: WPART ',WPART
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSE
         KPA = KPAO
         KPA2 = KPAO2
      ENDIF
	
      CALL PYROBO(0,N,0.D0,0.D0, DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +DBCMS(3)/DBCMS(4))
	
	
      if(wpart.ne.wpart) write(6,*) ' p_gamgam wpart = ',wpart,xsi_hard
c      write(6,*) ' p_gamgam end ,kpa,kpa1,kpa2,kpao2',kpa,kpa1,kpa2,kpao2
      RETURN
   90 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: ylimit ; RETURN ',yx1,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      RETURN
  100 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
         NDIMC = 9999
      ENDIF
      RETURN
  110 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: xlimit ; RETURN xmin,xmax',xmin,xmax
      ENDIF
      NDIMC = 9999
      RETURN
  120 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: xg2 limit ; RETURN ',xg2
      ENDIF
      NDIMC = 9999
      RETURN
  130 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: slimit ; RETURN ',shat1,am(1),am(2),IPRO
         write(6,*) ' sss,xg1,xg2,qg21,qg22 ',sss,xg1,xg2,qg21,qg22
         write(6,*) ' p_gamgam q21,y1,sss',q21,yx1,sss
         write(6,*) ' p_gamgam q22,y2,sss',q22,yx2,sss
         write(6,*) ' p_gamgam 2*pq ',2*dot1(nia1,2),sss*yx1
      ENDIF
      NDIMC = 9999
      RETURN
  140 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: theta limit ; RETURN ',themi,thema,thete
      ENDIF
      NDIMC = 9999
      RETURN
  150 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: bochck limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  160 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: shat ; RETURN ',shat,shat1,kpa
c        write(6,*) ' p_gamgam: QG2MAX,QG2 ',QG2MAX,QG21
c        write(6,*) ' p_gamgam: yx,xg2,q2,s ',yx,xg2,q2,sss
c        write(6,*) ' p_gamgam 2*pq ',2*dot1(nia1,2),sss*yx
c        write(6,*) ' p_gamgam nia1,nia2 ',nia1,nia2
c        write(6,*) ' p_gamgam 2*ktgam*ktglu ',2*(P(Nia1,1)*P(nia2,1)+
c     &  P(Nia1,2)*P(NIA2,2))
      ENDIF
      NDIMC = 9999
      RETURN
  170 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: ECM limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  180 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: PHASE WT=0 ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  191 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: PTCUT limit ; RETURN ',PT2,KPA
      ENDIF
      NDIMC = 9999
      RETURN

  190 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: QG21,QG22 limit ; RETURN ',
     &	   DSQRT(QG21),DSQRT(QG21)
      ENDIF
      NDIMC = 9999
      RETURN
  200 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_gamgam: xsi_hard limit ; RETURN ',xsi_hard
      ENDIF
      NDIMC = 9999
      RETURN


      END
