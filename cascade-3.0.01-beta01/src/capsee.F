      SUBROUTINE capsee
      IMPLICIT NONE
#include "caslujets.inc"
#include "casgki.inc"
#include "caspara1.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "caspara5.inc"
#include "casherup.inc"
#include "caspar.inc"
#include "casscal.inc"
	Integer Iel1,Iel2,Nph1,Nph2
	Common/cagaminfo/Iel1,Iel2,Nph1,Nph2
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar
	Integer NVMAX
	Parameter (NVMAX=5000)
      Double Precision X,Y,Q2LP
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /CASTOU/PARL(30),X,Y,Q2LP,LLST(40)
      DOUBLE PRECISION DETOT,DBETA,STHETA,SPHI
      COMMON /MYBOOST/ DETOT,DBETA(2,3),STHETA(2),SPHI(2)

      double Precision XPY,SH
      Integer KFL
      COMMON /MYPROC/KFL(3,2),XPY(2),SH
      Double Precision PS(nvmax,5),PSAVE(2,5)
      Double Precision QMAX,DOT1
      Double Precision XPRO
      Integer ks,ns
      COMMON /COLR/ KS(nvmax,5),NS
      double precision psc(nvmax,5)
      Integer ksc(nvmax,5)
      logical first
      data    first /.TRUE./

      Integer I,J,it,ipu1,ipu2,Ncall,Nprin,Ncheck,niph,nb2,nso
      Integer ipu1o,ipu2o,mstjo,irfl,irad,nia1t,niel
      Double Precision SHRAP,xgtest,chec,p1,p2,p3,p4
      Double Precision pc1,pc2,pc3,pc4,umas2,xr,dboost
      Double Precision MT1,MT2
	Integer MSTJO14,MSTJO41

      DATA NCALL/0/
c      DATA NPRIN/100/
      DATA NPRIN/0/
      DATA NCHECK/0/
	integer nev
	data nev/0/
      If(first) Then
         first = .FALSE.
      EndIf
	nev=nev+1
      DO I=1,40
         DO J=1,5
            KS(I,J) = 0
            PS(I,J) = 0.d0
         ENDDO
      ENDDO

c..hju
      NCALL = NCALL + 1
      IPY(48)=0
      LST(21)=0
C ILEPTO=1 parton shower a la LEPTO
C ILEPTO=0 parton shower a la PYTHIA
      SHRAP=SHH
      Q2LP = Q2
      NIPH = Nph1
      NB2 = Nph2
c change to e in + z direction
      DO 10  I=1,N
         P(I,3) = - P(I,3)
   10 continue
      UMAS2 = PYMASS(K(2,2))**2
c.hju parl(21) = 2P.k = invariant mass s
      If(Ires(1).eq.0) then 
        PARL(21) = 2.D0*DOT1(1,Nph2)
	Else
        PARL(21) = 2.D0*DOT1(Nph1,Nph2)+dot1(Nph1,Nph1)+dot1(nph2,nph2)
	Endif
      PARL(22) = 2.D0*DOT1(NB2,NIPH)  

      XPRO = XG
      PYVAR(31) = 1.
      IF(IRES(1).EQ.1) PYVAR(31) = XG1SF
	XPRO = XG1SF
      PYVAR(32) = XPRO
      XR = XPRO
      IF(NCALL.LE.NPRIN) then
         write(6,*) ' capee: 1st '
	   CALL PYLIST(1)
      ENDIF
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' capee: 1st '
         CALL PYLIST(1)
      ENDIF

C...Transform to  cms, boost parameters in double precision.
         DETOT=P(NIPH,4)+P(NB2,4)
         DBETA(2,1)=(P(NIPH,1)+P(NB2,1))/DETOT
         DBETA(2,2)=(P(NIPH,2)+P(NB2,2))/DETOT
         DBETA(2,3)=(P(NIPH,3)+P(NB2,3))/DETOT
         CALL PYROBO(0,0,0.D0,0.D0,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
         SPHI(2)=PYANGL(P(NIPH,1),P(NIPH,2))
         CALL PYROBO(0,0,0.D0,-SPHI(2),0.D0,0.D0,0.D0)
         STHETA(2)=PYANGL(P(NIPH,3),P(NIPH,1))
         CALL PYROBO(0,0,-STHETA(2),0.D0,0.D0,0.D0,0.D0)
         DO I=1,5
            PSAVE(1,I)=P(1,I)
            PSAVE(2,I)=P(2,I)
            psc(7,i) = p(Iel1,i)
            psc(9,i) = p(Iel2,i)
         ENDDO

       IF(IRES(1).EQ.1) Then
         DO I=1,5
            PSAVE(1,I)=P(1,I)
            PSAVE(2,I)=P(2,I)
         ENDDO
         IF(NCALL.LE.NPRIN) Then
            write(6,*) ' capee: before  pyvar boost '
            CALL PYLIST(1)
         Endif
         dboost = (PYVAR(31)-PYVAR(32)) /(PYVAR(31)+PYVAR(32))
C...Transform to gamma g cms, boost parameters in double precision.
         CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
      ENDIF


      IF(NCALL.LE.NPRIN) Then
         write(6,*) ' capee: after boost '
         CALL PYLIST(1)
      Endif
C...Save event record in  cms
      DO 20 I=1,N
         DO 20 J=1,5
            KS(I,J)=K(I,J)
   20 PS(I,J)=P(I,J)
C...Rearrange event record to PYSSPA standard
      DO 30 J=1,5
         K(3,J)=0
         P(3,J)=0.D0
         K(4,J)=0
         P(4,J)=0.D0
         K(5,J)=KS(NIA1,J)
         P(5,J)=PS(NIA1,J)
         K(6,J)=KS(NIA2,J)
         P(6,J)=PS(NIA2,J)
         IF(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13) THEN
            K(7,J)=KS(Iel1,J)
            P(7,J)=PS(Iel1,J)
            K(9,J)=KS(Iel2,J)
            P(9,J)=PS(Iel2,J)
         ENDIF
         K(8,J)=KS(NF1,J)
         P(8,J)=PS(NF1,J)

         K(11,J) = 0
         K(12,J) = 0
         K(13,J) = 0
         K(14,J) = 0

         K(10,J)=KS(NF2,J)
   30 P(10,J)=PS(NF2,J)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(10,3)=6
      DO 40 I=5,10
   40 K(I,1)=21
      N=15
C...Partons with colour information in hadronic cms frame.
      DO 60  I=16,MSTU(4)
         DO 60 J=1,5
            K(I,J)=0
            P(I,J)=0.D0
   60 V(I,J)=0.d0
      NS=20

      DO 70  J=1,5
         N=NS+7
         K(NS+1,J)=K(5,J)
         P(NS+1,J)=P(5,J)
         K(NS+3,J)=K(6,J)
         P(NS+3,J)=P(6,J)
         K(NS+5,J)=K(8,J)
         P(NS+5,J)=P(8,J)
         K(NS+6,J)=K(10,J)
         P(NS+6,J)=P(10,J)
         IF(IRES(1).NE.0) THEN
            K(NS+1,1) = 13
		Else
            K(NS+1,1) = 21		
         ENDIF
         IF(IRES(2).NE.0) THEN
            K(NS+3,1)=13
            K(NS+3,3)=2
            K(NS+3,4)=27
            K(NS+3,5)=27
		Else
            K(NS+3,1) = 21		
         ENDIF
   70 CONTINUE
C...Old standard continuation lines
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      P(NS+4,3)=27
      P(NS+4,4)=27
C...Colour info for two outgoing partons
      K(NS+5,1)=3
      K(NS+6,1)=3
      K(NS+5,4)=(NS+7)*MSTU(5)
      K(NS+5,5)=(NS+7)*MSTU(5)
      K(NS+6,4)=(NS+7)*MSTU(5)
      K(NS+6,5)=(NS+7)*MSTU(5)
C...Effective outgoing parton = sum of both outgoing partons
      K(NS+7,1)=14
      K(NS+7,3)=3
C... qqbar event
c      write(6,*) ' ires(1),ires(2) ',ires(1),ires(2)
      IF(IRES(1).EQ.0) THEN
         K(NS+7,2)=21
         IF(K(NS+5,2).GT.0) THEN
            K(NS+7,4)=(NS+3)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSE
            K(NS+7,4)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ENDIF
         IF(IRES(2).EQ.0) THEN
            K(NS+7,1)=13
            K(NS+7,2)=22
            K(NS+7,3)=22
            K(NS+7,4)=0
            K(NS+7,5)=0
            K(NS+5,4)=(NS+6)*MSTU(5)
            K(NS+5,5)=(NS+6)*MSTU(5)
            K(NS+6,4)=(NS+5)*MSTU(5)
            K(NS+6,5)=(NS+5)*MSTU(5)
	   ENDIF
      ELSE
c gg --> qqbar
         K(NS+7,2)=21
         K(NS+3,4)=21
         K(NS+3,5)=27
         IF(K(NS+5,2).GT.0) THEN
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=0
            K(NS+6,4)=0
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSE
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=0
            K(NS+5,4)=0
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ENDIF
      ENDIF
      DO 80  J=1,4
         P(NS+7,J)=P(8,J)+P(10,J)
   80 CONTINUE
      P(NS+7,5)=DSQRT(DMAX1(0.D0, P(NS+7,4)**2-P(NS+7,1)**2-P(NS+
     +7,2)**2- P(NS+7,3)**2))
      N=NS+7
c this is for jpsi production
      IF(IPRO.EQ.2) THEN
         K(NS+5,1)=1
         K(NS+7,4)=(NS+3)*MSTU(5)+NS+6
         K(NS+7,5)=(NS+3)*MSTU(5)+NS+6
      ENDIF

CPHI
C...Scale for bremsstrahlung etc.

      SH = SHRAP
      IPY(40)=10
      IPY(47)=N
C...Save quantities for later use.
      IPY(41)=K(1,2)
      IPY(42)=KINT(2,2)
      XPY(1)=1.
      XPY(2)=PYVAR(32)
      IF(IRES(1).EQ.1) THEN
         XPY(1)=PYVAR(31)
         IPY(41)=KINT(1,2)

      ENDIF
c.hju CALL PYSTFU(K(2,2),XR,Q2,XPQ)
      IF(XPY(2).GT.1.) write(6,*) 'capee : XPY(2) > 1. ',XPY(2)

      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' IPY(42),XPY(2)',
     +    IPY(42),XPY(2)
      ENDIF
      IPY(11)=1
      KFL(2,1)=KS(NIA1,2)
      KFL(2,2)=K(6,2)
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      KFL(3,1)=K(1,2)
      KFL(3,2)=K(27,2)
      PYVAR(2)=PARL(21)
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(5)=PYVAR(1)/2.


C...Generate timelike parton shower (if required)
      IF(IPY(13).EQ.1) THEN
         IF(ILHA.GE.1) then
            MSTJO41=MSTJ(41)
            MSTJ(41)=1
	   ENDIF
         IF(IFINAL.EQ.1) THEN
            MT1 = P(25,1)**2 + P(25,2)**2 + P(25,5)**2
            MT2 = P(26,1)**2 + P(26,2)**2 + P(26,5)**2
            QMAX = (MT1 + MT2)/2.d0
            QMAX = 4d0*QMAX
         ELSEIF(IFINAL.EQ.2) THEN
            QMAX=P(27,5)**2
         ELSEIF(IFINAL.EQ.3) THEN
            MT1 = P(25,5)**2
            MT2 = P(26,5)**2
            QMAX = (MT1 + MT2)/2.d0
            QMAX = 4d0*QMAX
         ENDIF
         QMAX = QMAX * SCALFAF
         write(6,*) 'caps final ',ifinal,scalfaf
         QMAX=DSQRT(QMAX)
         CALL PYSHOW(25,26,QMAX)
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      NS=N
	if(N.GT.NVMAX) Then
	  write(6,*) ' capee: N > NVMAX; N=',N,' NVMAX = ',NVMAX
	  Stop
	  Endif
      do I = 20,N
         do j=1,5
            psc(i,j) = p(i,j)
            ksc(i,j) = k(i,j)
         enddo
      enddo

C...Generate spacelike parton shower (if required)
      IPU1=0
changed
      IF(IRES(1).EQ.1) IPU1=21
      IPU2=23

      XGTEST = (P(21,5)**2 + P(IT,5)**2+PYPAR(22))/PARL(22)


c     call pylist(1)
      IF(Ires(1).Eq.0.and.Ires(2).EQ.0) then
	ELSE
      IF(IPY(14).GE.1.AND.XGTEST.LT.0.999) THEN
         NSO=NS
         IPU1O=IPU1
         IPU2O=IPU2
         CALL cascps(IPU1,IPU2)

         IF(LST(21).NE.0) THEN
c           write(6,*) ' capee after cascps: LST(21) = ',LST(21)
            IF(LST(21).EQ.55.OR.LST(21).EQ.58) THEN
            ELSE
               RETURN
            ENDIF
         ENDIF
      ENDIF
	ENDIF
      IF(.NOT.(IPY(14).GE.1.AND.XGTEST.LT.0.999)
     +   .OR.LST(21).EQ.55.OR.LST(21).EQ.58.OR.
     +   (Ires(1).Eq.0.and.Ires(2).EQ.0)) THEN
c      IF(.NOT.(IPY(14).GE.1.AND.XGTEST.LT.0.999)
c     +   .OR.LST(21).EQ.55.OR.LST(21).EQ.58) THEN
         IF(LST(21).EQ.55) LST(21)=0
         IF(LST(21).EQ.58) THEN
            PYVAR(32) = XPRO
            PYVAR(31) = XG1SF
            NS = NSO
            IPU1=IPU1O
            IPU2=IPU2O
            LST(21)=0
            do i=20,27
               do j=1,5
                  p(i,j) = psc(i,j)
                  k(i,j) = ksc(i,j)
               enddo
            enddo
            DO I=NS+1,N
               DO J=1,5
                  K(I,J)=0
                  P(I,J)=0.D0
                  V(I,J)=0.D0
               ENDDO
            ENDDO
            N=NS
         ENDIF

         DO 90  I=NS+1,NS+4
            DO 90  J=1,5
               K(I,J)=0
               P(I,J)=0.D0
   90    V(I,J)=0.d0
Cadded
         IF(IRES(1).EQ.0) THEN
            K(NS+1,1)=11
            K(NS+1,3)=21
            IF(IRES(2).EQ.0) THEN
               K(NS+1,1)=13
               K(NS+1,3)=21
               K(NS+3,1)=13
               K(NS+3,3)=21
		ENDIF
         ELSE
            K(NS+2,1)=-1
            K(NS+2,3)=NS+1
            P(NS+2,3)=21
            P(NS+2,4)=21
            P(22,1)=NS+1
            P(22,2)=NS+1
            K(NS+1,1)=13
            K(NS+1,3)=21
            K(NS+1,4)=21
            K(NS+1,5)=21
            K(21,4)=K(21,4) + (NS+1)*MSTU(5)
            K(21,5)=K(21,5) + (NS+1)*MSTU(5)
         ENDIF
         DO 100 J=1,5
  100    P(NS+1,J)=P(21,J)
         IF(IRES(1).EQ.0.AND.IRES(2).EQ.0) THEN
            K(NS+1,2)=KFL(2,1)
            K(NS+2,1)=-1
            K(NS+3,2)=KFL(2,2)
            K(NS+3,4)=0
            K(NS+3,5)=0
	   else
            K(NS+1,2)=KFL(2,1)
            K(NS+2,1)=-1
            K(NS+2,3)=NS+1
            K(NS+3,1)=13
            K(NS+3,2)=KFL(2,2)
            K(NS+3,3)=23
            K(NS+3,4)=23
            K(NS+3,5)=23
	   ENDIF
         P(NS+3,1)= 0.D0
         P(NS+3,2)= 0.D0
         P(NS+3,5)= 0.D0
         P(NS+3,1) = PS(NIA2,1)
         P(NS+3,2) = PS(NIA2,2)
         P(NS+3,5) = PS(NIA2,5)
         P(NS+3,3) = P(23,3)
         P(NS+3,4) = P(23,4)
         IF(NCALL.LT.NPRIN) THEN
            write(6,*) ' capee NS+3,23,it ',P(NS+3,3),P(23,3),P(IT,5)
            write(6,*) ' capee Q2 ',Q2,P(21,5)**2
         ENDIF

         K(NS+4,1)=-1
         K(NS+4,3)=NS+3
         P(NS+4,3)=23
         P(NS+4,4)=23
         P(24,1)=NS+3
         P(24,2)=NS+3
	   IF(Ires(1).EQ.0.and.Ires(2).eq.0) then
	   Else
            K(23,4)=K(23,4)+(NS+3)*MSTU(5)
            K(23,5)=K(23,5)+(NS+3)*MSTU(5)
	   Endif
         IPU1=0
         IF(IRES(1).EQ.1) IPU1=NS+1
         IPU2=NS+3
         N=N+4
      ENDIF
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' before CAREMN '
         CALL PYLIST(1)
      ENDIF
C...Hadron remnant and primordial kt
      IPY(47)=N

      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' capee before CAREMN: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF

      DO I=1,5
         P(1,I)=PSAVE(1,I)
         P(2,I)=PSAVE(2,I)
      ENDDO
	If(Ires(1).eq.0.and.Ires(2).eq.0) then
	Else
        CALL CAREMN(IPU1,IPU2)
	Endif

      IF(NCALL.LE.NPRIN) THEN
         write(6,*) 'after CAREMN IPY(48)=',IPY(48)
         CALL PYLIST(1)
      ENDIF
      IF(IPY(48).GE.1) THEN
         LST(21)=47+IPY(48)
         RETURN
      ENDIF
C...Rearrange partons along strings
      MSTU(24)=0
      MSTU(28)=0
      IF(IFPS.EQ.10) THEN
         MSTJO = MSTJ(105)
         MSTJ(105) = -1
      ENDIF
      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' capee before PYPREP: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
      CALL PYPREP(0)
      IF(IFPS.EQ.10) THEN
         MSTJ(105) = MSTJO
      ENDIF
      IF(MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN
c         WRITE(6,*) ' capee: LUPREP error MSTU(24)= ',MSTU(24),MSTU(28)
c     &   ,' Mass: ',P(27,5)
c         call PYLIST(2)
         LST(21)=50
         RETURN
      ENDIF

C...Clean up event record -> order:
C...1=inc. lepton; 2=inc. nucleon; 3=exch boson; 4=scat. lepton;
C... (4+1) = rad. gamma if any
C...5+1=inc. parton before initial shower; 6+1=inc. parton at hard scattering
C...after shower; 7+1,8+1=first,second parton from hard scattering
C...before final shower
c      IF(NCALL.LT.NPRIN) CALL PYLIST(1)
      DO 110 J=1,5
         K(N+2,J)=K(3,J)
         P(N+2,J)=P(3,J)
         K(N+1,J)=K(4,J)
         P(N+1,J)=P(4,J)
         K(10,J)=0
         K(N+3,J)=KS(32,J)
         P(N+3,J)=PS(32,J)
         K(12,J)=KS(32,J)
         P(12,J)=PS(32,J)

         K(N+4,J)=0
         P(N+4,J)=0.d0
         K(N+5,J)=0
         P(N+5,J)=0.d0
         K(N+4,J)=KS(NPH1,J)
         P(N+4,J)=P(1,J)-P(7,J)
         K(N+5,J)=KS(NPH2,J)
         P(N+5,J)=P(2,J)-P(9,J)
  110 CONTINUE
      P(N+4,5)=-DSQRT(DABS(DOT1(N+4,N+4)))
      P(N+5,5)=-DSQRT(DABS(DOT1(N+5,N+5)))
	K(N+4,3) = 1
	K(N+5,3) = 2
      DO 120 I=1,20
      DO 120 J=1,5
            PS(I,J) = P(I,J)
  120 KS(I,J) = K(I,J)
      IRFL = 0
      DO 130 J=1,5
         IRAD = 0
         IF(ILEPTO.EQ.1) THEN
            K(3,J)=K(N+4,J)
            P(3,J)=P(N+4,J)
            K(4,J)=K(N+5,J)
            P(4,J)=P(N+5,J)
            NIA1T = 3
            NIEL=5
            K(NIEL,J)=KS(7,J)
		K(NIEL,3) = 1
            P(NIEL,J)=psc(7,J)
            K(NIEL+1,J)=KS(9,J)
		K(NIEL+1,3) = 2
            P(NIEL+1,J)=psc(9,J)
            NIEL = 6
		NIA1 = NIA1T
            IF(IRES(1).EQ.1) THEN
               IRFL = 2
               K(NIEL+1,J) = K(N+2,J)
               K(NIEL+2,J) = K(21,J)
               K(NIEL+1,1) = 21
               K(NIEL+1,3) = NIA1T
               K(NIEL+2,1) = 21
               K(NIEL+2,3) = NIEL+1
               P(NIEL+1,J) = P(N+2,J)
               P(NIEL+2,J) = P(21,J)
               NIA1 = NIEL + 2 
            ENDIF
            K(Niel-1,1)=1
            K(Niel,1)=1

         ENDIF
         K(NIEL+1+IRFL,J)=K(N+1,J)
	   K(NIEL+1+IRFL,3)=NIA1T+1
         P(NIEL+1+IRFL,J)=P(N+1,J)
         K(NIEL+2+IRFL,J)=K(23,J)
         P(NIEL+2+IRFL,J)=P(23,J)
         NIA2 = NIEL + 2 + IRFL
         IF(Ires(1).EQ.0.and.Ires(2).EQ.0) then
            K(NIEL+1+IRFL,J)=K(21,J)
            K(NIEL+1+IRFL,3) = NIA1T
            P(NIEL+1+IRFL,J)=P(21,J)
            NIA1 = NIEL + 1 + IRFL
            K(NIEL+2+IRFL,3) = NIA1T+1
	   ELSE
            K(NIEL+2+IRFL,1) = 21	   
            K(NIEL+2+IRFL,3) = NIEL+1+IRFL 	   
	   Endif
         K(NIEL+3+IRFL,J)=K(25,J)
         P(NIEL+3+IRFL,J)=P(25,J)
         NF1 = NIEL + 3 + IRFL
         K(NIEL+4+IRFL,J)=K(26,J)
         P(NIEL+4+IRFL,J)=P(26,J)
         NF2 = NIEL + 4 + IRFL
  130 CONTINUE
      K(NIEL+2+IRFL,4)=0
      K(NIEL+2+IRFL,5)=0
	K(25,3)=NIA1
      K(NIEL+3+IRFL,4)=0
      K(NIEL+3+IRFL,5)=0
 	K(26,3)=NIA2
      K(NIEL+4+IRFL,4)=0
      K(NIEL+4+IRFL,5)=0
      IF(IRES(1).EQ.1) THEN
         P(NIEL+1,5)=-DSQRT(dabs(P(NIEL+1,1)**2 + P(NIEL+1,2)**2
     +   + P(NIEL+1,3)**2 - P(NIEL+1,4)**2))
         P(NIEL+2,5)=-DSQRT(dabs(P(NIEL+2,1)**2 + P(NIEL+2,2)**2
     +   + P(NIEL+2,3)**2 - P(NIEL+2,4)**2))
      ENDIF
      IF(iabs(k(1,2)).eq.2212.and.iabs(k(2,2)).eq.2212) then
         Nia1=4
         Nia2=6
         NF1=7
         NF2=8
      endif
C...Deactivate obsolete lines 9, 10, 21, NS+1 (extra lines with boson)
      K(21,1)=0
	IF(IRES(2).EQ.0) K(23,1)=0
      IF(IRFL.EQ.0) THEN
         K(NIEL-IRAD+6,1)=0
      ELSE

      ENDIF
      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
      K(13,1)=0
      K(14,1)=0
      K(15,1)=0
C...Zero irrelevant lines with K(I,1)<0
      DO 150 I=1,N
         IF(P(I,4).eq.0) THEN
            DO J=1,5
               K(I,J)=0
            ENDDO
         ENDIF
         IF(K(I,1).LT.0) THEN
            DO 140 J=1,5
               K(I,J)=0
  140       P(I,J)=0.D0
         ENDIF
  150 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' capee: before PYEDIT '
         CALL PYLIST(1)
      ENDIF
c         CALL PYLIST(1)
C...Delete internal parton lines, i.e. with K(I,1)=13,14
C...Delete lines from propagator gluons in initial state ps
      IF(IPST.EQ.0)  CALL pyedit(14)

      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' capee:after LUEDIT(14) '
         CALL PYLIST(1)
      ENDIF
C...Delete empty lines
      CALL PYEDIT(12)
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' capee:after LUEDIT(12) '
         CALL PYLIST(1)
      ENDIF
c check for enegy momentum conservation
      DO 170 I=1,N
         IF(P(I,5).LT.0.0D0.OR.K(I,1).GT.2) GOTO 170
         DO 160 J=1,4
            IF(ABS(P(I,J)).LE.1.E-6) THEN
               P(I,J)=0.0D0
            ENDIF
  160    CONTINUE
  170 CONTINUE
      CHEC = 5.D-3
      CHEC = 1.D-2
c      CHEC = 5.D-1
c      CHEC = 5.D-2
cccccc      CHEC = 5.D2
      P1 =ABS(P(1,1)+P(2,1))
      P2 =ABS(P(1,2)+P(2,2))
      P3 =ABS(P(1,3)+P(2,3))
      P4 =ABS(P(1,4)+P(2,4))
      PC1 = (ABS(PYP(0,1)) - P1)/P4
      PC2 = (ABS(PYP(0,2)) - P2)/P4
      PC3 = (ABS(PYP(0,3)) - P3)/P4
      PC4 = (ABS(PYP(0,4)) - P4)/P4
c      write(6,*) ' before check '
c      CALL PYLIST(1)
      IF(DABS(PC1).GT.CHEC.
     +   OR.DABS(PC2).GT.CHEC.
     +   OR.DABS(PC3).GT.CHEC.
     +   OR.DABS(PC4).GT.CHEC) THEN
         write(6,*) ' capee: energy of final particles not correct'
     +   ,chec
         write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
         write(6,*) ' Q2 = ',Q2,' Q2_calc = ',DOT1(3,3)
         write(6,*) ' yy = ',yy,' y_calc = ',DOT1(2,3)/DOT1(2,1)
         write(6,*) ' x_bj = ',Q2/PARL(22),' W = ',SQRT(-Q2+PARL(22))
         write(6,*) ' IPRO = ',IPRO
         CALL PYROBO(0,0,STHETA(2),SPHI(2),0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
         do  i=1,n
            p(i,3) = - p(i,3)
         enddo
         CALL PYROBO(0,0,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
         call PYLIST(1)
         LST(21) = 100
         return
      ENDIF
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' in cm system '
         CALL PYLIST(1)
      ENDIF

         CALL PYROBO(0,0,STHETA(2),0d0,0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0d0,SPHI(2),0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
c... change z coordinate  to CASCADE standard
      do 180 i=1,n
  180 p(i,3) = - p(i,3)
      DO 190 I=1,4
         IF(ABS(P(NIEL,I)).LE.1.E-6) THEN
            P(NIEL,I)=0.0D0
         ENDIF
         P(1,I)=PBEAM(1,I)
         P(2,I)=PBEAM(2,I)
  190 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' after boost dbeta in capee'
         CALL PYLIST(1)
      ENDIF
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' capee final '
         CALL PYLIST(1)
      ENDIF
c      write(6,*) ' capee final '
c      CALL PYLIST(1)
      RETURN
      END
