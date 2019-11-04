      SUBROUTINE ELEMNZ(WT1)
C
C    PHOTON GLUON ----> JPSI GLUON
C
C         P1(G)-----//////--------Q1(JPSI)
C                  ////////
C         P2(PH)-----//////--------Q2(G)
C
C  ( calculation of N.Zotov and V.Saleev
C
C
      IMPLICIT None
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casgki.inc"
#include "casshwr.inc"
      Double Precision MJ,MJ2
      Double Precision SS(4),UU(4),TT(4)
      Double Precision RRR,DOT,alphas,wpsi,k12,k22,yph,yph_bar
      Double Precision xsi_hard,wtg,PT2Q
      Double Precision WT1,SUM2,S,U,T,Z,AX,AJ,PSI2,QG2,ALPHS1,ALPHS2
      Double Precision XL,XK,XQ,EPSM,G
      Double Precision DOT1,dcasrn
      Double Precision dscale,alpha_s
      Double Precision FQ12,QG1,A1,A2,B1,B2
      Integer IKPSI,I
      Double Precision XPQ(-6:6),SFSCAL1,SFSCAL2,STMAX
      Double Precision CUT
c      DATA CUT /1.d-10/
c      DATA CUT /1.d-4/
       DATA CUT /1.d-3/
      Double Precision Scale
      Integer KF1,KF2,KIA1,KIA2,Kin1,Kin2
      Double Precision DSH(4),DSHN(4),SH_NEW,PCM_FIN,COS_FIN,PCM_NEW
      Double Precision cphi,cphio,e1,e2
      Double Precision QG0
      Integer Ikincut
      COMMON /GLUDAT/QG0,ikincut
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Data KPSI/443/
      EXTERNAL RRR,dcasrn
      Logical first
      Data first/.true./
      data wpsi/4.7D-6/
      Double Precision Wfact,Bochck
      SAVE cphio   
      
      IKPSI = -99
      if(first) then
c now the new value from 2000
         WPSI=5.26D-6
c       alphs = 0.3
c        WFACT = 1./(1.-5.33*ALPHS/PI)
c        WRITE(6,*) ' W_psi= ',WPSI,' correction factor = ',WFACT
c        WPSI =WPSI * WFACT
         first=.false.
      endif
c boost to ep cms
c      CALL PYROBO(0,N,0.D0,0.D0,DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
c     +   DBCMS(3)/DBCMS(4))
      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))

      DO 10  I=1,N
   10 IF(K(I,2).EQ.KPSI) IKPSI=I
      IF(IKPSI.LT.0) THEN
         WRITE(6,*) ' NO PSI FOUND IN ELEMDI; FATAL; STOP;'
         STOP
      ENDIF
      QG1 = DABS(DOT1(NIA1,NIA1))
      QG2 = DABS(DOT1(NIA2,NIA2))


      KIA1=N+1
      KIA2=N+2
      KF1=N+3
      KF2=N+4
      Kin1=N+5
      Kin2=N+6
      N=N+6
      DO I=1,4
         P(KIA1,I) = P(NIA1,I)
         P(KIA2,I) = P(NIA2,I)
         P(KF1,I) = P(NF1,I)
         P(KF2,I) = P(NF2,I)
         DSH(I)= P(NIA1,I) + P(NIA2,I)
      ENDDO
      K(KIA1,1)=13
      K(KIA2,1)=13
      K(KF1,1)=13
      K(KF2,1)=13

      K(Kin1,1)=13
      K(Kin2,1)=13
c define incoming particles as massless, to ensure proper cancellation in
c matrix-element
      P(Kin1,4)=Dabs(P(1,3))
      P(Kin1,5)=0
      P(Kin2,4)=Dabs(P(2,3))
      P(Kin2,5)=0
      Do I=1,3
         P(Kin1,I)=P(1,I)
         P(Kin2,I)=P(2,I)
      enddo


c     Goto 123
c
      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
      P(KIA1,4) = DABS(P(KIA1,3))
      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
      P(KIA2,4) = DABS(P(KIA2,3))

      IF(ICCFM.EQ.0) Then
         P(KIA1,1) = 0D0
         P(KIA1,2) = 0D0
         P(KIA2,1) = 0D0
         P(KIA2,2) = 0D0
      ENDIF

      DO I=1,4
         DSHN(I)= P(KIA1,I) + P(KIA2,I)
      ENDDO
      SH_NEW = DOT(DSHN,DSHN)
      IF(SH_NEW.LE.(AM(1)+AM(2))**2) GOTO 30

      CALL PYROBO(0,0,0.D0,0.D0,-DSH(1)/DSH(4), -DSH(2)/DSH(4), -DSH(3)
     +/DSH(4))
      PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
      COS_FIN=P(NF1,3)/PCM_FIN
      CALL PYROBO(0,0,0.D0,0.D0,DSH(1)/DSH(4) ,DSH(2)/DSH(4), DSH(3)/
     +DSH(4))
      CALL PYROBO(0,0,0.D0,0.D0,-DSHN(1)/DSHN(4), -DSHN(2)/DSHN(4),
     +-DSHN(3)/DSHN(4))

      P(KF1,4) =(SH_NEW + AM(1)**2)/2./DSQRT(SH_NEW)
      PCM_NEW=SQRT(P(KF1,4)**2 - AM(1)**2)

      P(KF1,3) =PCM_NEW*COS_FIN
      IF(IMEPHI.EQ.0) THEN
         cphi=cos(2.d0*PI*dcasrn())
         cphio = cphi
      ELSEIF(IMEPHI.EQ.1) THEN
         cphi=cphio
      ENDIF
      P(KF1,1) =DSQRT(PCM_NEW**2-P(KF1,3)**2)
      P(KF1,2) = 0.d0
      P(KF1,1) =DSQRT(PCM_NEW**2-P(KF1,3)**2)*cphi
      P(KF1,2) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*dsqrt(1.d0-cphi**2)
      DO I=1,4
         P(KF2,I) = P(KIA1,I) + P(KIA2,I) - P(KF1,I)
      ENDDO

      CALL PYROBO(0,0,0.D0,0.D0,DSHN(1)/DSHN(4) ,DSHN(2)/DSHN(4),
     +DSHN(3)/DSHN(4))

  123 Continue


      PT2Q = P(IKPSI,1)**2 + P(IKPSI,2)**2
      MJ = P(IKPSI,5)
      MJ2 = MJ*MJ
      DO 20 I=1,4
         SS(I)=P(KF1,I)+P(KF2,I)
         TT(I)=P(KIA1,I)-P(KF1,I)
         UU(I)=P(KIA2,I)-P(KF1,I)
   20 CONTINUE

      S  =DOT(SS,SS)
      U  =DOT(UU,UU)
      T  =DOT(TT,TT)
      K12 = -DABS(DOT1(KIA1,KIA1))
      K22 = -DABS(DOT1(KIA2,KIA2))

      E1 = -K12/AM(1)**2
      E2 = -K22/AM(1)**2
c	if(e2.lt.cut) then
c	  k22= -1.0
c	  E2 = -K22/AM(1)**2
c	  endif

      SCALE = PT2Q+MJ2
      IF(IQ2.EQ.1) THEN
         SCALE = MJ2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = S
      ELSEIF(IQ2.EQ.3) THEN
         SCALE = MJ2 + PT2Q
      ELSEIF(IQ2.EQ.4) THEN
         SCALE = Q2
      ELSEIF(IQ2.EQ.5) THEN
         SCALE = Q2 + PT2Q + MJ2
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K22)
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE.',IQ2,' STOP '
         STOP
      ENDIF
      Dscale = scale



      alpha_s = ALPHAS(sqrt(dscale))




      WT1=0.D0
      Qg2 = DABS(K22)
	
      if(E2.lt.cut) then
         XL = 1.D0 - S/MJ2
         XK = 1.D0 - T/MJ2
         XQ = 1.D0 - U/MJ2
         EPSM = 0.D0
c         EPSM = 0.16D0
         PSI2 = G(XL,XK,XQ,EPSM)+G(XK,XQ,XL,EPSM)+G(XQ,XL,XK,EPSM)
         PSI2=PSI2/MJ2/MJ2
         IF(IRUNA.EQ.1) ALPHS=ALPHA_S
C.....MARTIX ELEMENT GAMMA GLU --> JPSI GLU
         SUM2 = PSI2*128.D0*PI*PI*ALPHS*ALPHS*(MJ**3) * WPSI
         SUM2 = SUM2 / (3.D0*ALPH)
c	   write(6,*) ' should not enter '
c      elseif(E1.LT.CUT.AND.E2.GT.CUT) THEN
      elseif(E2.GT.CUT) THEN

         Z = DOT1(KF1,2)/DOT1(KIA1,2)
         AX = 1.D0 - Z
         AJ = Z
         AJ = DOT1(KF1,2)/DOT1(KIA1,2)
         AX = DOT1(KF2,2)/DOT1(KIA1,2)
C           WRITE(6,*) 'S+T+U',S+T+U
C           WRITE(6,*) 'M_PSI**2 - QG2=',MJ2 - QG2
C.....MARTIX ELEMENT
         PSI2 = RRR(S,T,U,AJ,AX,Dabs(k22),MJ)
         IF(PSI2.LT.0.0d0) THEN
C          WRITE(6,*) 'PSI2 = ',PSI2
C          WRITE(6,*) 'S+T+U  ',S+T+U
C          WRITE(6,*) 'M_PSI**2 - QG2=',MJ2 - QG2
C          WRITE(6,*) ' S ',S,' U ',U,' T ',T
C          WRITE(6,*) ' AX ',AX,' AJ ',AJ,' QG2 ',QG2
C          WRITE(6,*) ' Z = ',AJ
            PSI2=0.0d0
         ENDIF
         SCALE = PT2Q+MJ2
         alphs1=alpha_s
         alphs2 =alphs1
         PSI2 = PSI2 * ALPHS1 * ALPHS2
         SUM2 = PSI2 * 16.D0 * PI * PI * WPSI*MJ / 3.D0 / ALPH
C ... x*s_tot = s_hat + q_t**2
         SUM2 = SUM2 * (S + Dabs(k22))**2 / Dabs(k22)
c apply additional factor 1/pi according to N. Zotov (emails march 2000)
         sum2 = sum2/pi
c        sum2 = 0
      endif
      xg_bar = k22/(1.d0-xg)/SSS
      yph = YY
      yph_bar = k12/(1.d0 - yph)/SSS
      SFSCAL2 = DABS(K22)
      if(E2.gt.cut) then
         xg_bar =  k22/(1.d0-xg)/SSS
      else
         xg_bar=0
         SFSCAL2 = 0.01D0
         IF(ICCFM.EQ.0) SFSCAL2=QG2
c        write(6,*) ' E2 < cut ',k22
      endif
c      yph = dble(YY)
c      if(E1.gt.cut) then
c         yph_bar =  k12/(1.d0 - yph)/SSS
c      else
c         yph_bar = 0.d0
c         SFSCAL1 = 0.01D0
c         IF(ICCFM.EQ.0) SFSCAL1=QG1
c      endif

      xsi_hard = (yph+xg_bar)/(xg+yph_bar)
      if(xsi_hard.le.0) then
c            write(6,*) 'meoffsh xsi ',xsi_hard
c            write(6,*) 'k12 ',k12,' k22 ',K22
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
         stmax = 0.54D0
      else
         stmax = xg*sqrt(xsi_hard*SSS)
      endif

      if(stmax.ge.sqrt(sss)) then
c           write(6,*) ' meoffsh stmax ',stmax,' sqrt(s) ',sqrt(sss)
c            write(6,*) ' meoffsh xsi ',xsi_hard
c            write(6,*) ' k12 ',k12,' k22 ',k22
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
         XPQ(0) = 0.d0
      else
c            write(6,*) ' x,kt2,p ',xg,sfscal,stmax
         call cauniglu(K(2,2),XG,SFSCAL2,STMAX,XPQ)
      endif
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4), -CM(2)/CM(4),-CM(3)/CM(4)
     +)
c      CALL PYROBO(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),-DBCMS(2)/ DBCMS(4),
c     +-DBCMS(3)/DBCMS(4))
      WTG= dble(XPQ(0))
      WT1=SUM2*WTG
c     if(wt1.le.0) write(6,*) sum2,wtg,xg,xg_bar,k22,sfscal2,qg2
      CALL PYEDIT(13)
      RETURN
   30 continue
      WT1=0d0
      CALL PYEDIT(13)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RRR(S1,T1,U1,AJ,AX,G2,MJ)
C THIS IS FROM SALEEV END OF DECEMBER 1993
C     FUNCTION RRR(S1,T1,U1,AJ,AX,G2,MJ)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION MJ
      D1=ABS(S1-MJ**2)
      D2=ABS(T1-MJ**2)
      D3=ABS(U1-G2-MJ**2)
C
      R11=(G2*(AJ**2*MJ**2-2d0*AJ*AX*MJ**2+2d0*AJ*AX*T1-4d0*AX
     + **2*MJ**2+2d0*AX**2*T1)+(-3d0*AJ**2*MJ**4+2d0*AJ**2*MJ
     + **2*T1-10d0*AJ*AX*MJ**4+2d0*AJ*AX*MJ**2*T1-4d0*AJ*AX*MJ**2*
     + U1+2d0*AJ*AX*T1**2+2d0*AJ*AX*T1*U1-8d0*AX**2*MJ**4-2d0*AX**2*
     + MJ**2*T1-6d0*AX**2*MJ**2*U1+2d0*AX**2*T1**2+2d0*AX**2*T1*U1)
     + )/(D1**2*D2**2)
      R22=(G2**2*AJ*(-AJ-2d0*AX)+2d0*G2*(-AJ**2*T1-AJ*AX
     + *MJ**2-AJ*AX*T1+2d0*AX**2*MJ**2)+AJ*(-6d0*AJ*MJ**4+2d0*
     + AJ*MJ**2*T1+6d0*AJ*MJ**2*U1-AJ*T1**2+AJ*U1**2-4d0*AX*MJ**4-
     + 2d0*AX*MJ**2*T1+2d0*AX*MJ**2*U1+2d0*AX*T1*U1+2d0*AX*U1**2))/(
     + 2d0*D1**2*D3**2)
      R33=(2d0*G2*(2d0*AJ**2*MJ**2-AJ**2*T1+5d0*AJ*AX*MJ**2-AJ
     + *AX*T1+2d0*AX**2*MJ**2)+AJ*(AJ*MJ**4-2d0*AJ*MJ**2*T1-
     + AJ*T1**2+6d0*AX*MJ**4-2d0*AX*MJ**2*T1-6d0*AX*MJ**2*U1+2d0*AX*
     + T1*U1))/(2d0*D2**2*D3**2)
      R12=(-2d0*G2**2*AX**2+G2*(-AJ**2*MJ**2+4d0*AJ*AX*
     + MJ**2-2d0*AJ*AX*U1-2d0*AX**2*MJ**2-2d0*AX**2*T1-4d0*AX**2*U1)
     + +(6d0*AJ**2*MJ**4-3d0*AJ**2*MJ**2*T1-3d0*AJ**2*MJ**2*U1
     + +12d0*AJ*AX*MJ**4-4d0*AJ*AX*MJ**2*U1-2d0*AJ*AX*T1*U1-2d0*AJ*
     + AX*U1**2+4d0*AX**2*MJ**4+2d0*AX**2*MJ**2*T1-2d0*AX**2*MJ**2*
     + U1-2d0*AX**2*T1*U1-2d0*AX**2*U1**2))/(D1**2*D2*D3)
      R13=(2d0*G2*(-4d0*AJ**2*MJ**2+AJ**2*T1-8d0*AJ*AX*MJ**2+
     + 2d0*AJ*AX*T1-3d0*AX**2*MJ**2+AX**2*T1)+(-3d0*AJ**2*MJ**4
     + +AJ**2*MJ**2*T1-4d0*AJ*AX*MJ**2*T1-6d0*AJ*AX*MJ**2*U1+2d0*AJ
     + *AX*T1*U1+6d0*AX**2*MJ**4-2d0*AX**2*MJ**2*T1-6d0*AX**2*MJ**2
     + *U1+2d0*AX**2*T1*U1))/(D1*D2**2*D3)
      R23=(G2**2*(AJ**2-3d0*AJ*AX-2d0*AX**2)+G2*(9d0*AJ**
     + 2*MJ**2-2d0*AJ**2*T1-3d0*AJ**2*U1-2d0*AJ*AX*T1-4d0*AJ*AX*U1-
     + 4d0*AX**2*U1)+(3d0*AJ**2*MJ**4-4d0*AJ**2*MJ**2*U1-AJ**2*
     + T1**2-AJ*AX*MJ**4-2d0*AJ*AX*MJ**2*T1+2d0*AJ*AX*MJ**2*U1+2d0*
     + AJ*AX*T1*U1-AJ*AX*U1**2-2d0*AX**2*MJ**4+4d0*AX**2*MJ**2*U1-
     + 2d0*AX**2*U1**2))/(D1*D2*D3**2)
      RRR=R11+R22+R33+R12+R13+R23
      RETURN
      END
      DOUBLE PRECISION FUNCTION G(XL,XK,XQ,EDM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C     MATRIX ELEMENT A LA KEUNG MUZINICH PHYS REV D 27 (1518) 1983
C
C     INCLUDING RELATIVISTIC CORRECTIONS
C
C     EPSM = EPS/M
C     M_JPSI = 2*M + EPS
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      G=0.D0
      G1=
     +  (1+5.d0/3.d0*EDM)*((XL-1)/XQ/XK)**2
      G2= - EDM/12.d0 *
     + ( 8*(XQ**4+XK**4)-4*(XQ**3+XK**3)*(10+6*XQ*XK)
     +      +(XQ**2+XK**2)*(88+82*XQ*XK)
     +        -(XQ+XK)*(88+198*XQ*XK+21*XQ**2*XK**2)
     +          +32+228*XQ*XK+126*XQ*XQ*XK*XK)
      XR=XQ**3 *XK**3 *XL
      IF(XR.NE.0.) G2=G2/XR
      G=(G1 + G2)

      return
      END
      Double Precision function fq12(s0,qt1,qt2,a1,a2,b1,b2,m)
      Implicit None
      Double Precision s0,qt1,qt2,a1,a2,b1,b2,m
      Double Precision ans,ans1,ans2
      ans2=a2*b2*(-3.0*m**8+10.0*m**6*qt1+10.0*m**6*qt2-12.0*m**4*qt1
     + **2-20.0*m**4*qt1*qt2-12.0*m**4*qt2**2+6.0*m**2*qt1**3+10.0*m
     + **2*qt1**2*qt2+10.0*m**2*qt1*qt2**2+6.0*m**2*qt2**3-qt1**4+2.0
     + *qt1**2*qt2**2-qt2**4)
      ans1=2.0*s0*(-2.0*a1**2*b2**2*m**2*qt1*qt2+2.0*a1**2*b2**2*qt1
     + **2*qt2+3.0*a1*a2*b1*b2*m**6-7.0*a1*a2*b1*b2*m**4*qt1-7.0*a1*
     + a2*b1*b2*m**4*qt2+5.0*a1*a2*b1*b2*m**2*qt1**2+8.0*a1*a2*b1*b2*
     + m**2*qt1*qt2+5.0*a1*a2*b1*b2*m**2*qt2**2-a1*a2*b1*b2*qt1**3-a1
     + *a2*b1*b2*qt1**2*qt2-a1*a2*b1*b2*qt1*qt2**2-a1*a2*b1*b2*qt2**3
     + -2.0*a1*a2*b2**2*m**2*qt1*qt2+2.0*a1*a2*b2**2*qt1**2*qt2-2.0*
     + a1*a2*b2**2*qt1*qt2**2-2.0*a2**2*b1**2*m**2*qt1*qt2+2.0*a2**2*
     + b1**2*qt1*qt2**2-2.0*a2**2*b1*b2*m**2*qt1*qt2-2.0*a2**2*b1*b2*
     + qt1**2*qt2+2.0*a2**2*b1*b2*qt1*qt2**2-4.0*a2**2*b2**2*m**6+9.0
     + *a2**2*b2**2*m**4*qt1+9.0*a2**2*b2**2*m**4*qt2-6.0*a2**2*b2**2
     + *m**2*qt1**2-10.0*a2**2*b2**2*m**2*qt1*qt2-6.0*a2**2*b2**2*m**
     + 2*qt2**2+a2**2*b2**2*qt1**3-a2**2*b2**2*qt1**2*qt2-a2**2*b2**2
     + *qt1*qt2**2+a2**2*b2**2*qt2**3)+ans2
      ans=2.0*s0**3*a2**2*b2**2*m**2*(-a1**2*b1**2+2.0*a1*a2*b1*b2-a2
     + **2*b2**2)+s0**2*a2*b2*(-3.0*a1**2*b1**2*m**4+4.0*a1**2*b1**2*
     + m**2*qt1+4.0*a1**2*b1**2*m**2*qt2-a1**2*b1**2*qt1**2-2.0*a1**2
     + *b1**2*qt1*qt2-a1**2*b1**2*qt2**2-4.0*a1**2*b1*b2*qt1*qt2-4.0*
     + a1*a2*b1**2*qt1*qt2+10.0*a1*a2*b1*b2*m**4-12.0*a1*a2*b1*b2*m**
     + 2*qt1-12.0*a1*a2*b1*b2*m**2*qt2+2.0*a1*a2*b1*b2*qt1**2-4.0*a1*
     + a2*b1*b2*qt1*qt2+2.0*a1*a2*b1*b2*qt2**2+4.0*a1*a2*b2**2*qt1*
     + qt2+4.0*a2**2*b1*b2*qt1*qt2-7.0*a2**2*b2**2*m**4+8.0*a2**2*b2
     + **2*m**2*qt1+8.0*a2**2*b2**2*m**2*qt2-a2**2*b2**2*qt1**2+6.0*
     + a2**2*b2**2*qt1*qt2-a2**2*b2**2*qt2**2)+ans1
      fq12=-ans*16.*s0
      return
      end
