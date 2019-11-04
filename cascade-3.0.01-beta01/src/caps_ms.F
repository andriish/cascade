      SUBROUTINE caps_ms
      IMPLICIT NONE
#include "caslujets.inc"
#include "casgki.inc"
#include "caspara1.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "caspara5.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspar.inc"
#include "casscal.inc"

	Integer NVMAX
	Parameter (NVMAX=5000)
      Double Precision X,Y,Q2LP
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /CASTOU/PARL(30),X,Y,Q2LP,LLST(40)
      DOUBLE PRECISION DETOT,DBETA,STHETA,SPHI
      COMMON /MYBOOST/ DETOT,DBETA(2,3),STHETA(2),SPHI(2)

	Double Precision pyvar1_1,pyvar2_1,pyvar31_1,pyvar32_1
      Double Precision XPY,SH
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

      Integer I,J,it,ipu1,ipu2,ipu11,ipu21
	Integer Ncall,Nprin,Ncheck,niph,nb2,nso
      Integer ipu10,ipu20,ipu1ol,ipu2ol
	Integer mstjo,irfl,irad,nia1t,niel
      Double Precision SHRAP,xgtest,chec,p1,p2,p3,p4
      Double Precision pc1,pc2,pc3,pc4,umas2,xr,dboost
      Double Precision MT1,MT2
	Integer NSTB
	Double Precision DETOTR,DBETAR(2,3),STHETAR(2),SPHIR(2)
	Double Precision px,py
	Common/carespt/px(2),py(2)

	Double Precision multscatt,dummy,xmul_max,xmul_test
	Integer MSTJO14,MSTJO41
	Double Precision weig_ms
	Integer K_ms(100,5),In,Inc,Nmult
	Integer itest,itest2,k41,k42,k51,k52,k4,k5
	Integer it41,it42,it51,it52
	Double Precision P_ms(100,5)
	Double Precision xnbar
	Integer nma,ierr
	
	Integer In1
	
	
      DATA NCALL/0/
      DATA NPRIN/0/
      DATA NCHECK/0/
	integer nev
	data nev/0/
	ns_mult=0
	pyvar31_1=0.
	pyvar32_1=0.
      PYVAR(32) = 0
	PYVAR(31) = 0
	ipu10 = 0
	ipu11 = 0
	ipu20 = 0
	ipu21 = 0
	If(ires(1).eq.0) then
	  write(6,*) ' no multiple interaction for IRES(1)=0 '
	  stop
	endif
	nev=nev+1
      first = .TRUE.

c	write(6,*) ' new  event in caps-ms ' 
666   continue
cc	write(6,*) ' event caps-ms ',ns_mult
c	if(nev.ge.1000) then
c	  ncheck = 1
c	  ncall = 0
c	  nprin =100
c	  endif
cGB
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
      NIPH = NIA1
      NB2 = 2
      NB2  = NIA2
	IF(IABS(K(1,2)).EQ.11) THEN
         DO 11 I=1,N
         IF(K(I,1).EQ.21.AND.K(I,2).EQ.22.AND.K(I,3).EQ.1) NIPH = I
  11     CONTINUE
         NB2 = 2
	ENDIF
c change to e in + z direction
      DO 10  I=1,N
         P(I,3) = - P(I,3)
   10 continue
      UMAS2 = PYMASS(K(2,2))**2
c.hju parl(21) = 2P.k = invariant mass s
c.hjutest      PARL(21) = 2.D0*DOT1(1,NB2)
      PARL(21) = 2.D0*DOT1(1,2)
c.hju parl(22) = 2P.q P = proton q = photon
      PARL(22) = 2.D0*DOT1(NB2,NIPH)
c	write(6,*) ' caps-ms  nb2,niph ',NB2,NIPH
      PARL(21) =PARL(21)+PYMASS(K(1,2))**2 +UMAS2

	XPRO = XG
      PYVAR(32) = XPRO*(1.-PYVAR(32))
      XR = XPRO
	PYVAR(31) = XG1SF*(1.-PYVAR(31))
c      write(6,*) ' caps_ms ',pyvar(31),pyvar(32)
      IF(NCALL.LE.NPRIN) CALL PYLIST(1)
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' caps-ms: 1st '
         CALL PYLIST(1)
      ENDIF
c      write(6,*) ' caps-ms: 1st '
c      CALL PYLIST(1)

C...Transform to gamma p cms, boost parameters in double precision.
      IF(IRES(1).EQ.1.and.Iabs(K(1,2)).eq.11)THEN
         DETOT=P(NIPH,4)+P(NB2,4)
         DBETA(2,1)=(P(NIPH,1)+P(NB2,1))/DETOT
         DBETA(2,2)=(P(NIPH,2)+P(NB2,2))/DETOT
         DBETA(2,3)=(P(NIPH,3)+P(NB2,3))/DETOT
         CALL PYROBO(0,0,0.D0,0.D0,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
         SPHI(2)=PYANGL(P(NIPH,1),P(NIPH,2))
         CALL PYROBO(0,0,0.D0,-SPHI(2),0.D0,0.D0,0.D0)
         STHETA(2)=PYANGL(P(NIPH,3),P(NIPH,1))
         CALL PYROBO(0,0,-STHETA(2),0.D0,0.D0,0.D0,0.D0)
         NSTB = 5
         IF(NCALL.LE.NPRIN) then
            write(6,*) ' in gamma p  '
            CALL PYLIST(1)
         Endif

         DETOTR=P(NIA1,4)+P(NIA2,4)
         DBETAR(2,1)=(P(NIA1,1)+P(NIA2,1))/DETOTR
         DBETAR(2,2)=(P(NIA1,2)+P(NIA2,2))/DETOTR
         DBETAR(2,3)=(P(NIA1,3)+P(NIA2,3))/DETOTR


         dboost = (PYVAR(31)-PYVAR(32)) /(PYVAR(31)+PYVAR(32))

         CALL PYROBO(nstb,0,0d0,0d0,0.0d0,0.0d0,-dboost)


         DO I=1,5
            PSAVE(1,I)=P(1,I)
            PSAVE(2,I)=P(2,I)
         ENDDO
      ELSE
         DO I=1,5
            PSAVE(1,I)=P(1,I)
            PSAVE(2,I)=P(2,I)
         ENDDO
         DBETA(2,1)=0d0
         DBETA(2,2)=0d0
         SPHI(2)=0d0
         STHETA(2)=0d0
         dboost = (PYVAR(31)-PYVAR(32)) /(PYVAR(31)+PYVAR(32))
         DBETA(2,3)=dboost
C...Transform to gamma g cms, boost parameters in double precision.
         CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
      ENDIF
c         write(6,*) ' caps-ms: after boost '
c         CALL PYLIST(1)


      IF(NCALL.LE.NPRIN) Then
         write(6,*) ' caps-ms: after boost '
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
         IF(IABS(K(1,2)).EQ.11.OR.
     &	   IABS(K(1,2)).EQ.13) THEN
            K(7,J)=KS(3,J)
            P(7,J)=PS(3,J)
            K(9,J)=KS(3,J)
            P(9,J)=PS(3,J)
         ENDIF
         K(8,J)=KS(NF1,J)
         P(8,J)=PS(NF1,J)

         K(11,J) = 0
         K(12,J) = 0

         K(10,J)=KS(NF2,J)
   30 P(10,J)=PS(NF2,J)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(9,3)=5
      K(10,3)=6
      DO 40 I=5,10
   40 K(I,1)=21
      K(9,1)=0
C...Incoming parton
      DO 50 J=1,4
         P(6,J)=P(8,J)+P(10,J)-P(5,J)
c this redinition of incoming and outgoing parton is needed
c because in CAREMN there is a boost to correct for the long. frame in case
c there are different masses for outgoing parton and remnant.
   50 CONTINUE
c end of redefintion
      P(6,5)=PS(NIA2,5)
      K(6,2)=KS(NIA2,2)
      N=15
c      IF(NCALL.LE.NPRIN) THEN
c         CALL PYLIST(1)
c      ENDIF
C...Partons with colour information in hadronic cms frame.
      IF(FIRST) then
         DO 60  I=16,MSTU(4)
         DO 60 J=1,5
            K(I,J)=0
            P(I,J)=0.D0
   60    V(I,J)=0.d0
      Endif
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
   70 CONTINUE
C...Old standard continuation lines
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      P(NS+4,3)=27
      P(NS+4,4)=27
C...Origin and Colour info for incoming parton
      K(NS+3,1)=13
      K(NS+3,3)=2
      K(NS+3,4)=27
      K(NS+3,5)=27
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
c gg --> qqbar
      K(NS+7,2)=21
      K(NS+1,4)=27
      K(NS+1,5)=23
      K(NS+3,4)=21
      K(NS+3,5)=27
      IF(K(NS+5,2).GT.0) THEN
         K(NS+5,4)=(NS+7)*MSTU(5)
         K(NS+5,5)=0
         K(NS+6,4)=0
         K(NS+6,5)=(NS+7)*MSTU(5)
         K(NS+7,4)=(NS+1)*MSTU(5)+25
         K(NS+7,5)=(NS+3)*MSTU(5)+26
      ELSE
         K(NS+6,4)=(NS+7)*MSTU(5)
         K(NS+6,5)=0
         K(NS+5,4)=0
         K(NS+5,5)=(NS+7)*MSTU(5)
         K(NS+7,4)=(NS+1)*MSTU(5)+26
         K(NS+7,5)=(NS+3)*MSTU(5)+25
      ENDIF
      DO 80  J=1,4
         P(NS+7,J)=P(NS+1,J)+P(NS+3,J)
   80 CONTINUE
      P(NS+7,5)=DSQRT(DMAX1(0.D0, P(NS+7,4)**2-P(NS+7,1)**2-P(NS+
     +7,2)**2- P(NS+7,3)**2))
      N=NS+7

CPHI
C...Scale for bremsstrahlung etc.

      shrap = dot1(NS+7,NS+7)
      SH = SHRAP
      IPY(40)=10
      IPY(47)=N
C...Save quantities for later use.
      IPY(41)=K(1,2)
      IPY(42)=KINT(2,2)
      XPY(1)=1.
      XPY(2)=XR
      PYVAR(2)=PARL(21)
 	IF(IABS(K(1,2)).EQ.11) THEN
         PYVAR(2) = -Q2LP +PARL(22)
	   IPY(41)=22
	ENDIF	   
      XPY(1)=PYVAR(31)
      IF(XPY(2).GT.1.) write(6,*) 'caps-ms : XPY(2) > 1. ',XPY(2)

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

         QMAX=DSQRT(QMAX)
         CALL PYSHOW(25,26,QMAX)
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      NS=N
	if(N.GT.NVMAX) Then
	  write(6,*) ' caps-ms: N > NVMAX; N=',N,' NVMAX = ',NVMAX
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
      IPU1=21
      IPU2=23

      XGTEST = (P(21,5)**2+P(IT,5)**2+PYPAR(22))/PARL(22)


c     call pylist(1)
      IF(IPY(14).GE.1.AND.XGTEST.LT.0.999) THEN
         NSO=NS
         IPU1Ol=IPU1
         IPU2Ol=IPU2
         CALL cascps(IPU1,IPU2)
c         write(6,*) ' after cascps ',LST(21),IPU1,IPU2
c     call pylist(1)

         IF(LST(21).NE.0) THEN
c           write(6,*) ' caps-ms after cascps: LST(21) = ',LST(21)
            IF(LST(21).EQ.55.OR.LST(21).EQ.58) THEN
            ELSE
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF(.NOT.(IPY(14).GE.1.AND.XGTEST.LT.0.999)
     +   .OR.LST(21).EQ.55.OR.LST(21).EQ.58) THEN
         IF(LST(21).EQ.55) LST(21)=0
         IF(LST(21).EQ.58) THEN
            PYVAR(32) = XPRO
            NS = NSO
            IPU1=IPU1Ol
            IPU2=IPU2Ol
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
         K(NS+1,2)=KFL(2,1)
         DO 100 J=1,5
  100    P(NS+1,J)=P(21,J)
         K(NS+2,1)=-1
         K(NS+2,3)=NS+1
         K(NS+3,1)=13
         K(NS+3,2)=KFL(2,2)
         K(NS+3,3)=23
         K(NS+3,4)=23
         K(NS+3,5)=23
         P(NS+3,1)= 0.D0
         P(NS+3,2)= 0.D0
         P(NS+3,5)= 0.D0
         P(NS+3,1) = PS(NIA2,1)
         P(NS+3,2) = PS(NIA2,2)
         P(NS+3,5) = PS(NIA2,5)
         P(NS+3,3) = P(23,3)
         P(NS+3,4) = P(23,4)
         IF(NCALL.LT.NPRIN) THEN
            write(6,*) ' caps-ms NS+3,23,it ',
     &		P(NS+3,3),P(23,3),P(IT,5)
            write(6,*) ' caps-ms Q2 ',Q2,P(21,5)**2
         ENDIF

         K(NS+4,1)=-1
         K(NS+4,3)=NS+3
         P(NS+4,3)=23
         P(NS+4,4)=23
         P(24,1)=NS+3
         P(24,2)=NS+3
         K(23,4)=K(23,4)+(NS+3)*MSTU(5)
         K(23,5)=K(23,5)+(NS+3)*MSTU(5)

         IPU1=NS+1
         IPU2=NS+3
         N=N+4
c         write(6,*) ' before lulist '
c        call PYLIST(2)
      ENDIF
      IF(NCALL.LE.NPRIN) THEN
c         write(6,*) ' before CAREMN '
c         CALL PYLIST(1)
      ENDIF
c     write(6,*) ' caps-ms bef. remn PYVAR(31)= ',PYVAR(31)
c     write(6,*) ' caps-ms bef. remn PYVAR(32)= ',PYVAR(32)
C...Hadron remnant and primordial kt
      IPY(47)=N

      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' caps-ms before CAREMN: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
	
c      write(6,*) 'before CAREMN '
c      CALL PYLIST(1)

      DO I=1,5
         P(1,I)=PSAVE(1,I)
         P(2,I)=PSAVE(2,I)
      ENDDO
	
c multiple scattering ?
      If(Mult.eq.1.and.ns_mult.lt.90) then 
         first = .FALSE.
	   ipu10 = ipu1
	   ipu20 = ipu2
	   pyvar31_1 = pyvar(31)
	   pyvar32_1 = pyvar(32)
	   pyvar1_1 = pyvar(1)
	   pyvar2_1 = pyvar(2)
c check on number of interactions ...
         xnbar=1.	 
	   CALL RNPSSN(xnbar,Nma,IERR)
c	   write(6,*) ' nr of mult. interactions N = ',Nma
	   if(Nma.le.1) goto 777
c	   write(6,*) ' before mul-scatt:' 
c boost back to cms	   
         dboost = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
         CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,dboost)
c         write(6,*) ' 1st part ',ipu1,ipu2
c	   call pylist(2)
	   call pyedit(21)
	   Xmul_max = 10000.
	   ns_mult = 100
	   Xmul_test = multscatt(weig_ms)
c	   write(6,*) ' caps_ms dummy,weig_ms ',dummy,weig_ms 
	   if(weig_ms.ne.0) then 
c	      write(6,*) ' weight mul-scatt:',Xmul_test
c	      call pylist(1)
c	      call pyedit(22)
c	      write(6,*) ' after pyedit(22) '
c	      call pylist(1)
            goto 666
	   else 	   
	      Do I=1,N
	      Do J=1,5
	         P(I,J) = 0
	         K(I,J) = 0
	      enddo
	      enddo
	      call pyedit(22)
c boost back to parton parton cms	   
	      pyvar(1) = pyvar1_1
	      pyvar(2) = pyvar2_1 
	      pyvar(31) = pyvar31_1 
	      pyvar(32) = pyvar32_1 
            dboost = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
	      ipu1 = ipu10
	      ipu2 = ipu20 
            CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
c		write(6,*) ' old event record restored '
	   endif
777	   continue
	elseif(Mult.eq.1.and.ns_mult.gt.90) then
c	   write(6,*) ' after mul-scatt:',N 
c	   call pylist(2)
c boost back to cms	   
         dboost = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
         CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,dboost)
	   Nmult=N
	   ipu11 = ipu1
	   ipu21 = ipu2
c         write(6,*) ' 2nd part ',ipu1,ipu2
c	   call pylist(2)
	   do I =20,Nmult
	      do J=1,5
		K_ms(I,J) = K(I,J)
		P_ms(I,J) = P(I,J)
		enddo
         enddo
	   Do I=1,N
	   Do J=1,5
	      P(I,J) = 0
	      K(I,J) = 0
	   enddo
	   enddo
	   call pyedit(22)
	   If(N.le.50) then 
	      Inc = 50
		elseif(N.gt.50.and.N.le.100) then
		Inc = 100
		else
		Inc = 200
	   endif
cc	   K(21,5) = 0
cc	   K(23,4) = 0
	   In = Inc
	   Inc = Inc -20
c	   N=In
c	   In1 = In
c	   write(6,*) ' before rearrangement '		
c	   do I =21,Nmult
c	      In1 = In1 + 1
c	      do J=1,5
c		   K(In1,J) = K_ms(I,J) 
c	         P(In1,J) = P_ms(I,J)  
c            enddo
c         enddo
c	   N=In1

c	   call pylist(2)
	   do I =21,Nmult
	      In = In + 1
	      do J=1,5
		If(J.lt.3) then
		   K(In,J) = K_ms(I,J) 
c		   if(I.eq.25.or.I.eq.26) then
c		     K(In,2) = -K(In,2)
c		   Endif
		   elseif(j.eq.3) then 
		   K(In,J) = K_ms(I,J)+Inc
		   elseif(j.eq.4) then
		   itest=MOD(K_ms(I,4),MSTU(5))
		   itest2=MOD(K_ms(I,4)/10000,10000)
		   k41=0
		   k42=0
		   if(itest.ne.0) k41=itest+Inc
		   if(itest2.ne.0) k42=itest2+inc
		   K(In,4)=k41 + MSTU(5)*k42
c		   write(6,*) ' test K_ms(I,4) ',I,k_ms(I,4),itest,itest2
		   elseif(j.eq.5) then
		   itest=MOD(K_ms(I,5),MSTU(5))
		   itest2=MOD(K_ms(I,5)/10000,10000)
		   k51=0
		   k52=0
		   if(itest.ne.0) k51=itest+Inc
		   if(itest2.ne.0)  k52=itest2+Inc
		   K(In,5)=k51 + MSTU(5)*k52
		   k4=K(In,4)
		   k5=K(In,5)
		   if(I.eq.23) then
		     
		     K(In,4) = k5 
		     K(In,5) = k4 
ccc		     K(In,5) = 0
		   endif
		   if(I.eq.21) then
		     K(In,5) = k4 
 		     K(In,4) = k5
ccc		     K(In,4) = 0
		   endif
		endif 
		if(I.eq.27) then
		   it41=MOD(K_ms(I,4),MSTU(5))
		   it42=MOD(K_ms(I,4)/10000,10000)
		   k41=0
		   k42=0
		   if(it41.ne.0) k41=it41+Inc
		   if(it42.ne.0)  k42=it42+Inc
		   it51=MOD(K_ms(I,5),MSTU(5))
		   it52=MOD(K_ms(I,5)/10000,10000)
		   k51=0
		   k52=0
		   if(it51.ne.0) k51=it51+Inc
		   if(it52.ne.0) k52=it52+Inc
		   K(In,4)=k41 + MSTU(5)*k52
		   K(In,5)=k51 + MSTU(5)*k42
		   endif
		P(In,J) = P_ms(I,J)  
		enddo
         enddo
c         write(6,*) ' caps IPU11,21 ',IPU11,ipu21
         IPU11 = Inc + ipu11
         IPU21 = Inc + ipu21
c         K(Inc+27,4)=IPU11*MSTU(5)+ Inc + 26
c         K(Inc+27,5)=IPU21*MSTU(5)+ Inc + 25
	   N=In
	   
c	   write(6,*) ' after copy back',N 
c	   call pylist(2)
	   do J=1,5
	   P(N+1,J) = P(ipu10,J) + P(ipu11,J)
	   P(N+3,J) = P(ipu20,J) + P(ipu21,J)
	   Enddo
	   pyvar(31) = pyvar31_1 + pyvar(31)
	   pyvar(32) = pyvar32_1 + pyvar(32)
	   pyvar(1) = pyvar1_1
	   pyvar(2) = pyvar2_1 
	   inc = inc+20
	   K(N+1,1) = 13
	   K(N+1,2) = 21
	   K(N+1,3) = ipu10 
	   K(N+1,4) = ipu10  
	   K(N+1,5) = ipu11
	   ipu1 = N+1 
         K(N+2,1)=-1
         K(N+2,3)=N+1
         P(N+2,3)=21
         P(N+2,4)=21
 	   K(ipu10,1) =  14
 	   K(ipu10,4) =  (N+1)*mstu(5) + K(ipu10,4)
 	   K(ipu10,5) =  ipu11*mstu(5) + K(ipu10,5)
cc  	   K(ipu10,5) = 0
	   K(ipu11,1) =  14
	   K(ipu11,4) =  ipu10*mstu(5) + K(ipu11,4) 
	   K(ipu11,5) =  (N+1)*mstu(5) + K(ipu11,5)
cc 	   K(ipu11,5) =  0 
	   K(N+3,1) = 13
	   K(N+3,2) = 21
	   K(N+3,3) = ipu20 
	   K(N+3,4) = ipu21  
	   K(N+3,5) = ipu20  
	   ipu2 = N+3
         K(N+4,1)=-1
         K(N+4,3)=N+3
         P(N+4,3)=21
         P(N+4,4)=21
 	   K(ipu20,1) =  14
	   K(ipu20,5) =  (N+3)*mstu(5) + K(ipu20,5)
 	   K(ipu20,4) =  ipu21*mstu(5) + K(ipu20,4)
cc  	   K(ipu20,4) =  0
 	   K(ipu21,1) =  14
	   K(ipu21,5) =  ipu20*mstu(5) + K(ipu21,5) 
cc 	   K(ipu21,5) =  0
	   K(ipu21,4) =  (N+3)*mstu(5) + K(ipu21,4)
c	   write(6,*) ' ipus ',ipu10,ipu20,ipu11,ipu21
c	   write(6,*) ' inc  ',inc
	   N=N+4
c boost again to parton parton cms	   
         dboost = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
         CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
c	   write(6,*) ' after mult. scatt. color conn.',ipu1,ipu2
c	   write(6,*) ' ipu10,ipu20',ipu10,ipu20
c	   write(6,*) ' ipu11,ipu21',ipu11,ipu21
c	   call pylist(1)
	Endif
c	write(6,*) ' before caremn '
c	call pylist(2)
	
      CALL CAREMN(IPU1,IPU2)

c      write(6,*) 'after CAREMN '	
c      CALL PYLIST(2)
	
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) 'after CAREMN IPY(48)=',IPY(48)
         CALL PYLIST(1)
      ENDIF
      IF(IPY(48).GE.1) THEN
         LST(21)=47+IPY(48)
c        write(6,*) ' caps-ms after CAREMN: LST(21) = ',LST(21)
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
         WRITE(6,*) ' caps-ms before PYPREP: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
      CALL PYPREP(0)
      IF(IFPS.EQ.10) THEN
         MSTJ(105) = MSTJO
      ENDIF
      IF(MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN
         WRITE(6,*) ' caps-ms: PYPREP error MSTU(24)= ',MSTU(24)
	   write(6,*) ' ipu10,ipu20',ipu10,ipu20
	   write(6,*) ' ipu11,ipu21',ipu11,ipu21
         call PYLIST(2)
         LST(21)=50
         RETURN
	   else
c         WRITE(6,*) ' caps-ms: PYPREP ok  ',MSTU(24)
c	   write(6,*) ' ipu10,ipu20',ipu10,ipu20
c	   write(6,*) ' ipu11,ipu21',ipu11,ipu21
c         call PYLIST(2)
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

         IF(IRES(1).EQ.1.AND.IABS(K(1,2)).EQ.11) THEN
            K(N+4,J)=KS(NIPH,J)
            P(N+4,J)=PS(NIPH,J)
         ELSE
            K(N+4,J)=0
            P(N+4,J)=0.d0
         ENDIF
  110 CONTINUE
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
            NIA1T = 3
            NIEL=4
            K(NIEL,J)=KS(9,J)
            P(NIEL,J)=PS(9,J)
            IF(KS(13,1).NE.0) THEN
               K(NIEL+1,J)=KS(13,J)
               P(NIEL+1,J)=PS(13,J)
               IRAD = 1
               NIEL = 5
            ENDIF
            IF(IRES(1).EQ.1) THEN
               IRFL = 2
               K(NIEL+1,J) = K(N+2,J)

               K(NIEL+2,J) = K(21,J)
               K(NIEL+1,1) = 21
               K(NIEL+1,3) = NIEL-1
               K(NIEL+2,1) = 21
               K(NIEL+2,3) = NIEL+1
               P(NIEL+1,J) = P(N+2,J)

               P(NIEL+2,J) = P(21,J)
            ENDIF
            IF(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13) K(4,1)=1

         ENDIF
         K(NIEL+1+IRFL,J)=K(N+1,J)
         P(NIEL+1+IRFL,J)=P(N+1,J)

         K(NIEL+2+IRFL,J)=K(23,J)
         K(NIEL+2+IRFL,3) = NIEL+1+IRFL

         P(NIEL+2+IRFL,J)=P(23,J)
         NIA2 = NIEL + 2 + IRFL
         K(NIEL+3+IRFL,J)=K(25,J)
         P(NIEL+3+IRFL,J)=P(25,J)
         NF1 = NIEL + 3 + IRFL
         K(NIEL+4+IRFL,J)=K(26,J)
         P(NIEL+4+IRFL,J)=P(26,J)
         NF2 = NIEL + 4 + IRFL
  130 CONTINUE
c      write(6,*) ' NIEL,IRFL,NS ',NIEL,IRFL,NS
      K(3,3)=1
      K(4,3)=1
      K(NIEL+2+IRFL,1)=21
      K(NIEL+2+IRFL,3)=5+IRFL
      K(NIEL+2+IRFL,4)=0
      K(NIEL+2+IRFL,5)=0
      K(NIEL+3+IRFL,1)=21
      K(NIEL+3+IRFL,3)=NIA1T
      K(NIEL+3+IRFL,4)=0
      K(NIEL+3+IRFL,5)=0
      K(NIEL+4+IRFL,1)=21
      K(NIEL+4+IRFL,3)=NIA1T
      K(NIEL+4+IRFL,4)=0
      K(NIEL+4+IRFL,5)=0
      K(25,3) = NIEL + 3 + IRFL
      K(26,3) = NIEL + 4 + IRFL
      IF(IRES(1).EQ.1.AND.IABS(K(1,2)).EQ.2212) THEN
         P(NIEL+1,5)=-DSQRT(dabs(P(NIEL+1,1)**2 + P(NIEL+1,2)**2
     +   + P(NIEL+1,3)**2 - P(NIEL+1,4)**2))
         P(NIEL+2,5)=-DSQRT(dabs(P(NIEL+2,1)**2 + P(NIEL+2,2)**2
     +   + P(NIEL+2,3)**2 - P(NIEL+2,4)**2))
      ENDIF
      NIA1 = NIA1T
      IF(iabs(k(1,2)).eq.2212.and.iabs(k(2,2)).eq.2212) then
         Nia1=4
         Nia2=6
         NF1=7
         NF2=8
      endif
C...Deactivate obsolete lines 9, 10, 21, NS+1 (extra lines with boson)
      K(21,1)=0
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
         IF(K(I,1).LT.0) THEN
            DO 140 J=1,5
               K(I,J)=0
  140       P(I,J)=0.D0
         ENDIF
  150 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' caps-ms: before PYEDIT '
         CALL PYLIST(1)
      ENDIF
c         write(6,*) ' caps-ms: before PYEDIT '
c         CALL PYLIST(1)
C...Delete internal parton lines, i.e. with K(I,1)=13,14
C...Delete lines from propagator gluons in initial state ps
      IF(IPST.EQ.0)  CALL pyedit(14)


c         write(6,*) ' caps-ms: after LUEDIT(14) '
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' caps-ms:after LUEDIT(14) '
         CALL PYLIST(1)
      ENDIF
C...Delete empty lines
      CALL PYEDIT(12)
c         write(6,*) ' caps:after LUEDIT(12) '
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' caps:after LUEDIT(12) '
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
      CHEC = 5.D-1
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
         write(6,*) ' caps-ms: energy of final particles not correct'
     +   ,chec
         write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
         write(6,*) ' IPRO = ',IPRO	      
	   write(6,*) ' pyvar(31) = ',pyvar(31),pyvar31_1 
	   write(6,*) ' pyvar(32) = ',pyvar(32),pyvar32_1 

         IF(IRES(1).EQ.0) THEN
            CALL PYROBO(0,0,STHETA(2),SPHI(2),0.D0,0.D0,0.D0)
            CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
         ENDIF
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

c      write(6,*) ' caps-ms: 2nd '
c      CALL PYLIST(1)
cc      IF(IRES(1).EQ.0) THEN
      IF(IRES(1).EQ.0.OR.(IRES(1).EQ.1.and.IABS(K(1,2)).EQ.11)) THEN
	   IF(IRES(1).EQ.1.and.IABS(K(1,2)).EQ.11) then
c         CALL pyrobo(nstb,0,STHETAR(2),0.D0,0.D0,0.D0,0.D0)
c         CALL pyrobo(nstb,0,0.D0,SPHIR(2),0.D0,0.D0,0.D0)
	   endif
c         write(6,*) ' after rot,bef boost '
c         CALL PYLIST(1)
         CALL PYROBO(0,0,STHETA(2),0d0,0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0d0,SPHI(2),0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
      ENDIF
c      write(6,*) ' after boost '
c      CALL PYLIST(1)
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
         write(6,*) ' after boost dbeta in caps-ms'
         CALL PYLIST(1)
      ENDIF
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' caps-ms final '
         CALL PYLIST(1)
      ENDIF
c      write(6,*) ' caps-ms final '
c      CALL PYLIST(1)
      RETURN
      END
