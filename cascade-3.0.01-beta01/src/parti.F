      SUBROUTINE PARTI(IP1,YX,YVAL,WEIGHT,IFL,IST)
      IMPLICIT NONE
      Integer IP1,IFL,IST
      Double Precision YX,YVAL,WEIGHT
#include "caslujets.inc"
#include "casluco.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "caspart.inc"
#include "casgki.inc"
#include "casvalues.inc"
#include "caspara5.inc"
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar
      Double Precision PHI,DBGAM
      COMMON/DIFFA/ PHI
      DOUBLE PRECISION ME
      DIMENSION DBGAM(4)
      Double Precision ALPH_EM,FGAM,X11,COST2,SPHE,CPHE
      Double Precision PEP,PEZ,PDLE,XPP,BOCHCK,PHIE,XGP
      Double Precision SPHI,STHETA,PT,DOT1
      Double Precision yg,yg_bar
      Integer NIPH
C IST = 0, ONLY PARTICLE MOMENTA, BUT SCALE Q2Q IN STRUCTURE FUNCTION
C NOT YET DEFINED.
C IST = 1, CALCULATE ONLY STRUCTURE FUNCTION SINCE Q2Q NOW DEFINED
C
C IF PARTICLE IP1 = 11 (ELECTRON) USE EQUIVALENT PHOTON APPROXIMATION
C
c      write(6,*) 'parti  Q2 = ',Q2,' yx = ',yx
c      write(6,*) 'parti IRES(1) = ',IRES(1)
      YVAL =  -999999.D0
      WEIGHT = -999999.D0
      ME = PYMASS(IP1)
      ALPH_EM = ALPH
	if(IFL.EQ.1) Then
	  XGP=XG1
	  ELSE
	  XGP=XG2
	  ENDIF
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)

C IF 1ST INCOMING PARTICLE IS A ELECTRON AND LOOKING FOR DIRECT PHOTON
      IF((IABS(IP1).EQ.11.OR.IABS(IP1).EQ.13)) THEN
C
C...          WEIZAECKER WILLIAMS APPROXIMATION FOR BREMSTRAHLUNG (EPA)
C...          WITH Q2 DEPENDENCE
         FGAM = (1.D0 - YX +YX*YX/2.D0)/Q2/YX - ME*ME*YX/Q2/Q2

         YVAL = ALPH_EM * FGAM / PI
         WEIGHT = YX*Q2*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
C...

C...      CALCULATE COS OF ELECTRON VIA FORMULA GIVEN IN CAMPBELL ET AL
C...      CAN.J.PHYS. 39 1981 (1742)
C...      VIA X Q2 Y
C...      X11= Q2/(2P.Q)
         X11= Q2/YX/SSS
         COST2= ((1.D0 - YX) - X11*YX)/((1.D0 - YX) + X11*YX)
         IF(COST2.GT.1.) WRITE(6,*) 'YX Q2 COST2 ',YX,Q2,COST2
C...      COST2 IS SCATTERING ANGLE OF ELECTRON
         PHIE = PHI
         IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
         SPHE=DSIN(PHIE)
         CPHE=DCOS(PHIE)
C...      PEL1 = 4 VECTOR OF SCATTERED ELECTRON IN EP CMS SYSTEM
         N=N+1
         K(N,1)=1
         IF(INTER.LT.2) THEN
            K(N,2)=K(1,2)
         ELSEIF(INTER.EQ.2) THEN
            K(N,2)=ISIGN(1,K(1,2))*12
         ENDIF
         K(N,3)=IFL
c. 2*P.l_e
c new gymnastics for numerical stability
         PEP = DSQRT(P(2,3)**2 + P(2,5)**2)
         PEZ = P(2,3)
         PDLE = 2.D0*DOT1(1,2)
         X11 = Q2/PDLE/YX
         P(N,5)=PYMASS(K(N,2))
         P(N,4) = P(1,4) + (PEZ*(Q2+P(N,5)**2)/ (2.D0*P(1,4)) - YX*
     +   PDLE/2.D0)/(PEP+PEZ)
         P(N,3) = P(1,4) - (PEP*(Q2+P(N,5)**2)/ (2.D0*P(1,4)) + YX*
     +   PDLE/2.0D0)/(PEP+PEZ)
         P(N,3)=-P(N,3)
         PT = Q2 - (Q2 + P(N,5)**2)/P(1,4)/(PEP+PEZ)* (YX*PDLE/2.D0 +
     +   P(2,5)**2*(Q2+P(N,5)**2)/4.D0/ P(1,4)/(PEP+PEZ))
         P(N,1)=DSQRT(DMAX1(0.D0,PT))*CPHE
         P(N,2)=DSQRT(DMAX1(0.D0,PT))*SPHE

C...      PPH = 4 VECTOR OF (VIRTUAL) PHOTON IN EP CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KEB
         K(N,3)=IFL
         P(N,1)= P(1,1) - P(3,1)
         P(N,2)= P(1,2) - P(3,2)
         P(N,3)= P(1,3) - P(3,3)
         P(N,4)= P(1,4) - P(3,4)
         P(N,5)= -SQRT(ABS(DOT1(N,N)))
c            write(6,*) 'parti: ',P(N,5)**2,q2
         NIPH = N
         IF(IFL.EQ.1) THEN
            NIA1 = N
         ELSE
	      NIA2 = N
	   ENDIF
         IF(IRES(IFL).EQ.0) THEN
	      
	      IF(IFL.EQ.1) THEN
               NIR1 = -99999
		   XGP = yx
		ELSE
		   NIR2 = -99999
		   XGP = yx
		ENDIF
         ELSEIF(IRES(IFL).EQ.1) THEN
	      IF(IFL.EQ.1) THEN
		   XGP = xg1
		ELSE
		   XGP = xg2
		ENDIF
C         XGP = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XPP = E_GLUON/E_PHOTON
            XPP=XGP/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
            N = N + 1
            K(N+1,1)=21
            K(N+1,2)=21
            K(N+1,3)=N-1
            P(N+1,1)= P(N-1,1)*XPP
            P(N+1,2)= P(N-1,2)*XPP
            P(N+1,3)= P(N-1,3)*XPP
            P(N+1,4)= DABS(P(N-1,4))*XPP
            P(N+1,5)= PYMASS(21)
		IF(IFL.EQ.1) THEN
               NIA1 = N+1
		ELSE
		   NIA2 = N+1
		ENDIF
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
            K(N,1)=1
            K(N,2)=21
            K(N,3)=N-1
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,5)= PYMASS(21)
            P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
            P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
		IF(IFL.EQ.1) THEN
               NIR1 = N 
		ELSE
		   NIR2 = N 
		ENDIF

c            write(6,*) ' parti yx,xgp,xpp',yx,xgp,xpp
c            call PYLIST(1)
C NOW BOOST TO GAMMA PROTON FRAME
            DBGAM(1) = P(NIPH,1) + P(2,1)
            DBGAM(2) = P(NIPH,2) + P(2,2)
            DBGAM(3) = P(NIPH,3) + P(2,3)
            DBGAM(4) = P(NIPH,4) + P(2,4)
            BOCHCK = (DBGAM(1)/DBGAM(4))**2 + (DBGAM(2)/DBGAM(4))**2 +
     +      (DBGAM(3)/DBGAM(4))**2
            BOCHCK = DSQRT(BOCHCK)
            IF(BOCHCK.GT.0.99999999D0) goto 10
            CALL PYROBO(0,N,0.D0,0.D0,-DBGAM(1)/DBGAM(4),-DBGAM(2)/
     +      DBGAM(4), -DBGAM(3)/DBGAM(4))
            SPHI = PYANGL(P(NIPH,1),P(NIPH,2))
            call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
            STHETA = PYANGL(P(NIPH,3),P(NIPH,1))
            call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
            P(N+1,1) = 0.D0
            P(N+1,2) = 0.D0
            P(N+1,3) = XPP*P(NIPH,3)
            P(N+1,4) = XPP*(P(NIPH,4)*P(2,4)-P(NIPH,3)*P(2,3))/
     +                 (P(2,4)-P(2,3))
            P(N+1,3) = P(N+1,4)
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,4)= P(1,4) - P(N+1,4)- P(3,4)
            P(N,5)= DSQRT(DABS(P(N,1)**2+P(N,2)**2+P(N,3)**2-P(N,4)**2))
            N = N + 1
c            call PYLIST(1)
            call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
            call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
            CALL PYROBO(0,N,0.D0,0.D0,DBGAM(1)/DBGAM(4),DBGAM(2)/
     +      DBGAM(4), DBGAM(3)/DBGAM(4))
c            call PYLIST(1)
         ENDIF
C IF 1ST INCOMING PARTICLE IS A PHOTON
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(IFL).EQ.0) THEN
C...      4 VECTOR OF (REAL) PHOTON IN CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KPH
         K(N,3)=IFL
         P(N,1)= P(1,1)
         P(N,2)= P(1,2)
         P(N,3)= P(1,3)
         P(N,4)= P(1,4)
         P(N,5)=PYMASS(KPH)
	   IF(IFL.EQ.1) THEN 
            NIA1 = N
            NIR1 = -99999
	   ELSE
            NIA2 = N
            NIR2 = -99999
	   ENDIF
         YVAL = 1.D0
         WEIGHT = 1.D0
c	   write(6,*) ' in parti ',ires(1),(p(Nia1,i),i=1,5)
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(1).EQ.1) THEN
         WEIGHT = 1.0D0
         YVAL = 1.0D0

C         XGp = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XPp = E_GLUON/E_PHOTON
         XPP=XGP/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
         N = N + 1
         K(N+1,1)=21
         K(N+1,2)=21
         K(N+1,3)=IFL
         P(N+1,1)= P(N-1,1)*XPP
         P(N+1,2)= P(N-1,2)*XPP
         P(N+1,3)= P(N-1,3)*XPP
         P(N+1,4)= DABS(P(N-1,4))*XPP
         P(N+1,5)= PYMASS(21)
	   IF(IFL.EQ.1) THEN 
            NIA1 = N+1
	   ELSE
            NIA2 = N+1
	   ENDIF
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
         K(N,1)=1
         K(N,2)=21
         K(N,3)=IFL
         P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
         P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
         P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
         P(N,5)= PYMASS(21)
         P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
         P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
	   IF(IFL.EQ.1) THEN 
            NIR1 = N
	   ELSE
            NIR2 = N 
	   ENDIF
         N = N + 1
c            write(6,*) ' parti yx,xgp,xpp',yx,xgp,xpp
c            call PYLIST(1)
      ELSEIF(IABS(IP1).EQ.2212.AND.IRES(IFL).EQ.1) THEN
c no factor of pi here, this is accounted for in the structure function routine
         WEIGHT = YX*Q2*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
         YVAL = 1.D0/YX
         PHIE = PHI
c	   write(6,*) yx,q2,phi
         IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
         N=N+1
         K(N,1)=1
         K(N,2)=K(1,2)
         K(N,3)=IFL
C...      PPH = 4 VECTOR OF gluon IN PP CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KEB
         K(N,3)=IFL
         YG = YX
         YG_BAR = -Q2/(1.d0 - YG)/SSS

         P(N,4) =  yg*P(1,4) + yg_bar*P(2,4)
         P(N,3) =  yg*P(1,3) + yg_bar*P(2,3)
         PHIE = PHI
         P(N,1) = DSQRT(Q2)*SIN(PHIE)
         P(N,2) = DSQRT(Q2)*COS(PHIE)
         P(N,5)= -SQRT(ABS(DOT1(N,N)))
c            write(6,*) 'parti: ',P(N,5)**2,q2
c         write(6,*) ' parti ',KEB,K(1,2),IFL
	   IF(IFL.EQ.1) THEN 
            NIA1 = N
	   ELSE
            NIA2 = N
	   ENDIF
         P(N-1,1)= P(1,1) - P(N,1)
         P(N-1,2)= P(1,2) - P(N,2)
         P(N-1,3)= P(1,3) - P(N,3)
         P(N-1,4)= P(1,4) - P(N,4)
         P(N-1,5)= -SQRT(ABS(DOT1(N-1,N-1)))
         NIPH = N
	   IF(IFL.EQ.1) THEN 
         NIR1 = -99999
	   ELSE
         NIR2 = -99999
	   ENDIF  

      ENDIF
      RETURN
   10 write(6,*) ' PARTI boost error '
      RETURN
      END
