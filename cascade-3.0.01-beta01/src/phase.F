      SUBROUTINE PHASE(NP,ET,XM,PCM,WT)
c to generate n-body phase space
c taken from: Collider Physics
c             V.D. Barger and R.J.N Phillips
c             Adison-Wesley Publishing Company, Inc.
c             1987
	Implicit none
	Integer NP
	Double Precision ET,XM,PCM,WT
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
      DIMENSION XM(18),PCM(4,18)
	Integer I,J,Jmax
      N=NP
      PIE=4.D0*DATAN(1.D0)
      WPS=0.D0
      DO 10 I=1,N
         M(I)=0.D0
         MT(I)=0.D0
         MS(I)=0.D0
         DO 10 J=1,4
            K(I,J)=0.D0
            P(I,J)=0.D0
   10 CONTINUE
      DO 20 I=1,N
         M(I)=XM(I)
   20 CONTINUE
      K(1,1)=ET
      MT(1) =ET
      DO 30 I=1,N
   30 MS(1) = MS(1) + M(I)
      CALL ORDER
      JMAX = N - 1
      WPS = JMAX*((MT(1)-MS(1))/(4.D0*PIE**2))**(N-2)/
     +      (4.D0*PIE*MT(1))
      DO 40 J =1,JMAX
         WPS = WPS/J
         MS(J+1) = MS(J) - M(J)
         CALL GEN(J)
   40 CALL BOO(J)
      DO 50 J=1,N
         PCM(1,J) = P(J,2)
         PCM(2,J) = P(J,3)
         PCM(3,J) = P(J,4)
         PCM(4,J) = P(J,1)
   50 CONTINUE
      if(p(1,2).ne.p(1,2)) write(6,*) ' error in phase'
      WT=WPS*(2.D0*PIE)**(3*N-4)
c      write(6,*) ' phase ',wt
      RETURN
      END
      SUBROUTINE GEN(I)
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
	Double Precision X,AALAM
      COMMON/XVAR/ X(10)
#include "casparas.inc"
      Double Precision PSQ,EI,PI,COST,COSTM,SINT,PHI
      Integer I
      IF(I.EQ.N-1) THEN
         MT(I+1) = M(N)
      ELSE
         MT(I+1) = MS(I+1) + (MT(1)-MS(1))*X(3*I)
      ENDIF
      PSQ = AALAM(MT(I)**2,M(I)**2,MT(I+1)**2)/(4.D0*MT(I)**2)
      IF(PSQ.LT.0.0D0) WRITE(6,*) 'PSQ.LT.0 ',PSQ
      IF(PSQ.LT.0.0D0) PSQ=0.D0
      EI = DSQRT(PSQ+M(I)**2)
      PI = DSQRT(PSQ)
      COST = 1.D0-2.D0*X(3*I-1)
      SINT = DSQRT(1.D0-COST**2)
      PHI = 2.D0*PIE*X(3*I-2)
      P(I,1) = EI
      P(I,2) = PI*SINT*DSIN(PHI)
      P(I,3) = PI*SINT*DCOS(PHI)
c      write(6,*) ' GEN: PT**2 ',P(I,2)**2+P(I,3)**2,I
      P(I,4) = PI*COST
      if(p(I,2).ne.p(i,2)) then
         write(6,*) '  error in gen '
         write(6,*) ' pi,sint,phi ',pi,sint,phi,X(3*I-2),X(3*I-1)
      endif
      WPS = WPS*PI
      RETURN
      END
      SUBROUTINE ORDER
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K	
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
      IF(N.LE.3) RETURN
      WRITE(6,*) ' N .GT. 3; N = ',N
      STOP
      END
      SUBROUTINE BOO(I)
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K	
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
	Integer I,J
	Double Precision BG1,G1,BG2,G2,BG3,G3,E,X,Z,Y
      BG1 = K(I,2)/MT(I)
      G1 = DSQRT(1.D0 + BG1**2)
      BG2 = K(I,3)/(MT(I)*G1)
      G2 = DSQRT(1.D0 + BG2**2)
      BG3 = K(I,4)/(MT(I)*G1*G2)
      G3 = DSQRT(1.D0 + BG3**2)
      E = G1*G2*G3*P(I,1) + BG1*G2*G3*P(I,2)
     +    + BG2*G3*P(I,3) + BG3*P(I,4)
      X = BG1*P(I,1) +G1*P(I,2)
      Y = BG2*G1*P(I,1) + BG2*BG1*P(I,2) +G2*P(I,3)
      Z = BG3*G2*G1*P(I,1) + BG3*G2*BG1*P(I,2)
     +    + BG3*BG2*P(I,3) +G3*P(I,4)
      P(I,1) = E
      P(I,2) = X
      P(I,3) = Y
      P(I,4) = Z
      IF(I.EQ.N-1) GOTO 20
      DO 10 J=1,4
   10 K(I+1,J) = K(I,J) - P(I,J)
C     WRITE(6,*) 'BOO ',(P(I,IK),IK=1,4)
      RETURN
   20 DO 30 J =1,4
   30 P(N,J) = K(I,J) - P(I,J)
      RETURN
      END
