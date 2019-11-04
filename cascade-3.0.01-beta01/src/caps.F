      SUBROUTINE caps
      IMPLICIT NONE
#include "caslujets.inc"
#include "casgki.inc"
#include "caspara1.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "casludat2.inc"
#include "caspara5.inc"
#include "casherup.inc"
#include "caspsipol.inc"
#include "caspar.inc"
#include "casscal.inc"
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

      Integer Iremn
      common/casflags/Iremn	
      
      Integer NVMAX
	Parameter (NVMAX=5000)
      Double Precision X,Y,Q2LP
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /CASTOU/PARL(30),X,Y,Q2LP,LLST(40)
      DOUBLE PRECISION DETOT,DBETA,STHETA,SPHI
      COMMON /MYBOOST/ DETOT,DBETA(2,3),STHETA(2),SPHI(2)
      Double Precision XPY,SH
      Integer KFL
      COMMON /MYPROC/KFL(3,2),XPY(2),SH
      Double Precision PS(nvmax,5),PSAVE(2,5)
      Double Precision QMAX,DOT1,DOT
      Double Precision XPRO
      Integer ks,ns
      COMMON /COLR/ KS(nvmax,5),NS
      double precision psc(nvmax,5)
      Integer ksc(nvmax,5)
      logical first
      data    first /.TRUE./

      Integer I,J,it,ipu1,ipu2,Ncall,Nprin,Ncheck,niph,nb2,nso
      Integer ipu1o,ipu2o,mstjo,irfl,irad,nia1t,niel
      Double Precision SHRAP,chec,p1,p2,p3,p4
      Double Precision pc1,pc2,pc3,pc4,umas2,xr,dboost
      Double Precision MT1,MT2
	Integer NSTB
	Double Precision DETOTR,DBETAR(2,3),STHETAR(2),SPHIR(2)
	Double Precision px,py
	Common/carespt/px(2),py(2)

	Double Precision multscatt,dummy,xmul_max,xmul_test
      Double Precision PT1,PT2
	Integer MSTJO14,MSTJO41,MSTJO45
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar
      
      double precision cmspart(4),E1,E2,qm1,qm2,sh_test
      Integer MINT
      Double Precision VINT
      COMMON/PYINT1/MINT(400),VINT(400)
      Double Precision PARP,PARI
      INTEGER MSTP,MSTI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)	
      
      Integer Nps1,Nps2,Iveto,Iggveto
      Integer Nacc
        
      Data iggveto/1/
      DATA NCALL/0/
      DATA NPRIN/0/
      DATA NCHECK/0/
	integer nev
	data nev/0/
      If(first) Then
         first = .FALSE.
C...Calculates full and effective widths of gauge bosons, stores
C...masses and widths, rescales coefficients to be used for
C...resonance production generation.
         call pyinre
      end if
	nev=nev+1
      IPY(48)=0
      LST(21)=0
      NCALL = NCALL + 1
      if(ILHA.eq.-2) then
         call cas_lheread
         If(nup.le.0) then
            MSTU(24)=10
            Return
         endif 
         call cas_lheprocess
         return
      endif
      DO I=1,40
         DO J=1,5
            KS(I,J) = 0
            PS(I,J) = 0.d0
         ENDDO
      ENDDO

C ILEPTO=1 parton shower a la LEPTO
C ILEPTO=0 parton shower a la PYTHIA
      SHRAP=SHH
      Q2LP = Q2
      NIPH = NIA1
      NB2 = 2
      IF(IRES(1).EQ.1) THEN
         NB2  = NIA2
	   IF(IABS(K(1,2)).EQ.11) THEN
            DO 11 I=1,N
            IF(K(I,1).EQ.21.AND.K(I,2).EQ.22.AND.K(I,3).EQ.1) NIPH = I
  11        CONTINUE
		NB2 = 2
	   ENDIF
      ENDIF
c change to e in + z direction
      DO 10  I=1,N
         P(I,3) = - P(I,3)
   10 continue
      UMAS2 = PYMASS(K(2,2))**2
c.hju parl(21) = 2P.k = invariant mass s
      PARL(21) = 2.D0*DOT1(1,2)
c.hju parl(22) = 2P.q P = proton q = photon
      PARL(22) = 2.D0*DOT1(NB2,NIPH)
      IF(IRES(1).EQ.1) PARL(21) =PARL(21)+PYMASS(K(1,2))**2 +UMAS2

      XPRO = XG
      XR = XPRO
c redefinition of PYVAR(31),PYVAR(32) needed to ensure energy mom conservation in
c case of large virtualities
      IF(IRES(1).EQ.1) Then
         PYVAR(31) = xg1+xg2_bar
	ENDIF
      PYVAR(32) = xg2+xg1_bar
      
      IF(NCALL.LE.NPRIN) CALL PYLIST(1)
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' caps: 1st nev=',nev
         CALL PYLIST(1)
      ENDIF

C...Transform to gamma p cms, boost parameters in double precision.
      IF(IRES(1).EQ.0)THEN
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
         ENDDO
      ELSEIF(IRES(1).EQ.1.and.Iabs(K(1,2)).eq.11)THEN
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
         dboost = (PYVAR(31)-PYVAR(32))/
     &            (PYVAR(31)+PYVAR(32))
         DBETA(2,3)=dboost
C...Transform to g g cms, boost parameters in double precision.
         CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
      ENDIF


      IF(NCALL.LE.NPRIN) Then
         write(6,*) ' caps: after boost '
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
         IF(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13) THEN
            K(7,J)=KS(3,J)
            P(7,J)=PS(3,J)
            K(9,J)=KS(3,J)
            P(9,J)=PS(3,J)
         ENDIF
         K(8,J)=KS(NF1,J)
         P(8,J)=PS(NF1,J)

         if(NF1+1.EQ.NF2) then
         K(11,J) = 0
         K(12,J) = 0
         else
         K(11,J) = KS(NF1+1,J)
         P(11,J) = PS(NF1+1,J)
         endif

C check for radiative gamma
         IF(KS(4,2).EQ.22.AND.KS(4,1).EQ.1) THEN
            K(13,J) = KS(4,J)
            P(13,J) = PS(4,J)
         ENDIF
         K(10,J)=KS(NF2,J)
   30 P(10,J)=PS(NF2,J)
c      write(6,*) ipro,K(1,2)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(9,3)=5
      K(10,3)=6
      DO 40 I=5,10
   40 K(I,1)=21
      K(9,1)=0
      if(ipro.eq.4.and.(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13)) then
        K(10,1) = 0
        K(10,2) = 0
      endif
C...Incoming parton
      DO 50 J=1,4
         P(6,J)=P(8,J) + P(10,J) - P(5,J)
         IF(IRES(1).NE.0) P(6,J)=PS(NIA2,J)
c this redinition of incoming and outgoing parton is needed
c because in CAREMN there is a boost to correct for the long. frame in case
c there are different masses for outgoing parton and remnant.

   50 CONTINUE
      IF(NCALL.LE.NPRIN) Then
         write(6,*) ' caps: check '
         CALL PYLIST(1)
      Endif
c end of redefintion
      P(6,5)=PS(NIA2,5)
      K(6,2)=KS(NIA2,2)
      N=15
      if(ipsiel1.ge.1) then
        do j=1,5
          p(14,j) = PS(ipsiel1,j)
          k(14,j) = KS(ipsiel1,j)
          p(15,j) = PS(ipsiel2,j)
          k(15,j) = KS(ipsiel2,j)
        end do
      endif
C...Partons with colour information in hadronic cms frame.
      DO 60  I=16,MSTU(4)
         DO 60 J=1,5
            K(I,J)=0
            P(I,J)=0.D0
   60 V(I,J)=0.d0
      NS=20
         IF(NCALL.LE.NPRIN) then
            write(6,*) ' after rearragen  ',nia1,nia2,nf1,nf2
            CALL PYLIST(1)
         Endif

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
         ENDIF
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
	If(ipro.eq.3.or.ipro.eq.102.or.ipro.eq.500) K(NS+5,1)=1
      if(ipro.eq.2.and.ipsiel1.ge.1) K(NS+5,1)=11
      K(NS+5,4)=(NS+7)*MSTU(5)
      K(NS+5,5)=(NS+7)*MSTU(5)
      K(NS+6,4)=(NS+7)*MSTU(5)
      K(NS+6,5)=(NS+7)*MSTU(5)
      IF(IRESPRO.EQ.3) THEN
c qg event
         IF(K(NS+5,2).GT.0) THEN
            K(NS+5,4)=(NS+6)*MSTU(5)
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+5)*MSTU(5)
         ELSE
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=(NS+6)*MSTU(5)
            K(NS+6,4)=(NS+5)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
         ENDIF
      ENDIF

C...Effective outgoing parton = sum of both outgoing partons
      K(NS+7,1)=14
      K(NS+7,3)=3
C... qqbar event
      IF(IRES(1).EQ.0) THEN
         K(NS+7,2)=21
         IF(K(NS+5,2).GT.0) THEN
            K(NS+7,4)=(NS+3)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSE
            K(NS+7,4)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ENDIF
	   IF(Ipro.eq.4) then
            K(NS+6,1)=0
            K(NS+6,2)=0
            K(NS+6,3)=0
            K(NS+6,4)=0
            K(NS+6,5)=0
         endif
      ELSE
	   IF(Ipro.eq.3.or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then
            K(NS+3,4) = 21
            K(NS+3,5) = 21
            K(NS+5,4)=0
            K(NS+5,5)=0
            K(NS+6,1)=0
            K(NS+6,2)=0
            K(NS+6,3)=0
            K(NS+6,4)=0
            K(NS+6,5)=0
            K(NS+1,4) = 23
            K(NS+1,5) = 23
            K(NS+7,4)=0
            K(NS+7,5)=0
	   else
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
	   IF(IRESPRO.EQ.3) THEN
c... qg event
            K(NS+7,2)=K(NS+5,2)
            IF(K(NS+7,2).GT.0) THEN
               K(NS+7,4)=(NS+3)*MSTU(5)+26
               K(NS+7,5)=(NS+3)*MSTU(5)+25
            ELSE
               K(NS+7,4)=(NS+3)*MSTU(5)+25
               K(NS+7,5)=(NS+3)*MSTU(5)+26
            ENDIF
            Call Colorflow
	   ELSEIF(IRESPRO.EQ.2) THEN
            Call Colorflow
         ENDIF
	   Endif
      ENDIF
      DO 80  J=1,4
         P(NS+7,J)=P(8,J)+P(10,J)
         IF(IRES(1).NE.0) P(NS+7,J)=P(NS+1,J)+P(NS+3,J)
   80 CONTINUE
      P(NS+7,5)=DSQRT(DMAX1(0.D0, P(NS+7,4)**2-P(NS+7,1)**2-P(NS+
     +7,2)**2- P(NS+7,3)**2))
      N=NS+7
c this is for jpsi production
      IF(IPRO.EQ.2.or.ipro.eq.4) THEN
         call colorflow
      ENDIF
      IF(IPRO.eq.501.or.ipro.eq.511) THEN
         call colorflow
      ENDIF
678   continue           
c      write(6,*) ' in caps '
c      call pylist(1)
      
c write out LHA event record
      if(ILHA.le.-1) then
         call caupevnt
         if(ILHA.eq.-1) return
      endif
      
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
      IF(IRES(1).EQ.1) then
 	   IF(IABS(K(1,2)).EQ.11) THEN
            PYVAR(2) = -Q2LP +PARL(22)
	      IPY(41)=22
	   ENDIF	   
         XPY(1)=PYVAR(31)
	Endif

      IF(XPY(2).GT.1.) write(6,*) 'caps : XPY(2) > 1. ',XPY(2)

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
	
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(5)=PYVAR(1)/2.

      IF(ILHA.GE.1) then
         MSTJO41=MSTJ(41)
 	ENDIF

C...Generate timelike parton shower (if required)
      Nps1 =N
      IF(IPY(13).EQ.1) THEN
         IF(ILHA.GE.1) then
            MSTJO41=MSTJ(41)
            MSTJ(41)=1
	   ENDIF
         PT1 = P(25,1)**2 + P(25,2)**2
         PT2 = P(26,1)**2 + P(26,2)**2
         IF(IFINAL.EQ.1) THEN
            MT1 = PT1 + P(25,5)**2
            MT2 = PT2 + P(26,5)**2
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
c check for heavy flavor production in PS
         MSTJO45 = MSTJ(45)
cc         MSTJ(45) = 3
cc         if((PT1.gt.4.*PMAS(4,1)**2).or.(PT2.gt.4.*PMAS(4,1)**2)) MSTJ(45)=4
cc         if((PT1.gt.4.*PMAS(5,1)**2).or.(PT2.gt.4.*PMAS(5,1)**2)) MSTJ(45)=5        
c                
         QMAX = QMAX * SCALFAF
         QMAX=DSQRT(QMAX)
cc         write(6,*) ' before final state PS ', N, Qmax
ccc         write(6,*) ' caps ',QMAX
         if(ipro.eq.4) then
           CALL PYSHOW(25,0,QMAX)
         
         else
           CALL PYSHOW(25,26,QMAX)
         endif 
         MSTJ(45) = MSTJO45
      ENDIF
c      call pylist(1)
      Nps2 = N
cc      write(6,*) ' after final state PS ', N
CC start gg veto check      
      If(iggveto.eq.1) then
        Iveto = 0
        If(IPRO.EQ.10.and.IRESPRO.EQ.2) then
           Do I=Nps1,Nps2
           if(I.gt.0) then 
              if(iabs(K(I,2)).le.6.and.K(I,1).eq.3) then
                 PT2 = PYP(I,9)
                 CALL CUTG(PT2,NACC)
                 if(NACC.EQ.1) then
C veto those events where 2 quarks are above threshold, as they are 
C treated in gg->qqbar            
                    Iveto = Iveto + 1
                 Endif
              endif
           endif
           end do
        endif
      endif
      If(Iveto.ge.2) Iveto = 1
CC end gg veto check      
      MINT(51) = 0
      MSTP(71) = IPY(13)
c check if there are resonances to decay
      if((K(25,2).gt.22.and.K(25,2).le.25).or.iabs(K(25,2)).eq.6) Then
        CALL PYRESD(25)
      Endif
      if((K(26,2).gt.22.and.K(26,2).le.25).or.iabs(K(25,2)).eq.6) then 
        CALL PYRESD(26)  
      Endif
      IF(MINT(51).EQ.1) then
             write(6,*) ' caps: problem with res decay '
      endif
      IT=25
      IF(N.GE.27) IT=27
      NS=N
	if(N.GT.NVMAX) Then
	  write(6,*) ' caps: N > NVMAX; N=',N,' NVMAX = ',NVMAX
	  Stop
	  Endif
      do I = 20,N
         do j=1,5
            psc(i,j) = p(i,j)
            ksc(i,j) = k(i,j)
         enddo
      enddo
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' before PS '
         CALL PYLIST(1)
      ENDIF
c      write(6,*) ' aftter FPS '
c      CALL PYLIST(1)
      
C...Generate spacelike parton shower (if required)
      IPU1=0
changed
      IF(IRES(1).EQ.1) IPU1=21
      IPU2=23

      IF(IPY(14).GE.1) THEN
         NSO=NS
         IPU1O=IPU1
         IPU2O=IPU2
         CALL cascps(IPU1,IPU2)

         IF(LST(21).NE.0) THEN
            IF(LST(21).EQ.55.OR.LST(21).EQ.58) THEN
            ELSE
               RETURN
            ENDIF
         ENDIF
      ENDIF

      IF(.NOT.(IPY(14).GE.1)
     +   .OR.LST(21).EQ.55.OR.LST(21).EQ.58) THEN
         IF(LST(21).EQ.58.OR.LST(21).EQ.55) THEN
         IF(LST(21).EQ.55) LST(21)=0
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

         IF(IRES(1).EQ.0) THEN
            K(NS+1,1)=11
            K(NS+1,3)=21
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
            write(6,*) ' caps NS+3,23,it ',P(NS+3,3),P(23,3),P(IT,5)
            write(6,*) ' caps Q2 ',Q2,P(21,5)**2
         ENDIF

         K(NS+4,1)=-1
         K(NS+4,3)=NS+3
         P(NS+4,3)=23
         P(NS+4,4)=23
         P(24,1)=NS+3
         P(24,2)=NS+3
         K(23,4)=K(23,4)+(NS+3)*MSTU(5)
         K(23,5)=K(23,5)+(NS+3)*MSTU(5)
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
         WRITE(6,*) ' caps before CAREMN: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
	
      
      
c      write(6,*) 'before CAREMN '
c      CALL PYLIST(2)

      DO I=1,5
         P(1,I)=PSAVE(1,I)
         P(2,I)=PSAVE(2,I)
      ENDDO
      
c
      
      CALL CAREMN(IPU1,IPU2)
	
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
         WRITE(6,*) ' caps before PYPREP: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
      IF(ILHA.GE.1) then
	   MSTJO14=MSTJ(14)
	   MSTJ(14)=-1   ! to avoid pyprep producing clusters
	ENDIF
      if (Iremn.eq.1) CALL PYPREP(0)
      IF(ILHA.GE.1) Then 
        MSTJ(14) = MSTJO14
        MSTJ(41) = MSTJO41
      ENDIF
      IF(IFPS.EQ.10) THEN
         MSTJ(105) = MSTJO
      ENDIF
      IF(MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN
         WRITE(6,*) ' caps: PYPREP error MSTU(24)= ',MSTU(24),
     &     ' MSTU(28) = ',MSTU(28)
         call PYLIST(2)
         LST(21)=50
         RETURN
      ENDIF

C...Clean up event record -> order:
      NIEL = 0
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

         IF(IRES(1).EQ.0) THEN
            K(N+4,J)=KS(NIPH,J)
            P(N+4,J)=PS(NIPH,J)
         ELSEIF(IRES(1).EQ.1.AND.IABS(K(1,2)).EQ.11) THEN
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
      IF(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13) then
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
            K(4,1)=1
         ENDIF
         K(NIEL+1+IRFL,J)=K(N+1,J)
         P(NIEL+1+IRFL,J)=P(N+1,J)

         K(NIEL+2+IRFL,J)=K(23,J)
         K(NIEL+2+IRFL,3) = NIEL+1+IRFL

         P(NIEL+2+IRFL,J)=P(23,J)
         NIA2 = NIEL + IRFL
         K(NIEL+3+IRFL,J)=K(25,J)
         P(NIEL+3+IRFL,J)=P(25,J)
         NF1 = NIEL + 2 + IRFL
         K(NIEL+4+IRFL,J)=K(26,J)
         P(NIEL+4+IRFL,J)=P(26,J)
         NF2 = NIEL + 3 + IRFL
  130    CONTINUE
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
         NIA1 = NIA1T
         IF(IRES(1).EQ.1.AND.IABS(K(1,2)).EQ.2212) THEN
            P(NIEL+1,5)=-DSQRT(dabs(P(NIEL+1,1)**2 + P(NIEL+1,2)**2
     +      + P(NIEL+1,3)**2 - P(NIEL+1,4)**2))
            P(NIEL+2,5)=-DSQRT(dabs(P(NIEL+2,1)**2 + P(NIEL+2,2)**2
     +      + P(NIEL+2,3)**2 - P(NIEL+2,4)**2))
         ENDIF
         IF(IRFL.EQ.0) THEN
            K(NIEL-IRAD+6,1)=0
         ENDIF
	   IF(Ipro.eq.4) then 
            k(nf2,1)=-1
            nf2=-999
	      K(25,1)=0
	      K(26,1)=0
            K(NIEL+4,1)=0
         endif
      ELSEIF(IABS(K(1,2)).EQ.2212.and.IABS(K(1,2)).EQ.2212) then
         Nia1=3
         Nia2=5
         NF1=7
         NF2=8
	   K(25,3) = Nia1+1
	   K(26,3) = Nia2+1
	   do J=1,5
            K(3,J)=K(N+2,J)
            P(3,J)=P(N+2,J)
            K(Nia1+1,J)=K(21,J)
            P(Nia1+1,J)=P(21,J)
            K(5,J)=K(N+1,J)
            P(5,J)=P(N+1,J)
            K(Nia2+1,J)=K(23,J)
            P(Nia2+1,J)=P(23,J)
            K(NF1,J)=K(25,J)
            P(NF1,J)=P(25,J)
            K(NF2,J)=K(26,J)
            P(NF2,J)=P(26,J)		
	   enddo
	   K(3,1)=21
	   K(3,3)=1
	   K(4,1)=21
	   K(4,3)=3
	   K(5,1)=21
	   K(5,3)=2
	   K(6,1)=21
	   K(6,3)=5
         if(ipro.ne.3.and.IPRO.ne.102) then
           k(nf1,1)=21  !!!!!
           k(nf2,1)=21  !!!!!
         endif

	   IF(Ipro.eq.3..or.ipro.eq.4.or.ipro.eq.102.or.ipro.eq.500) then 
            k(nf2,1)=-1
            nf2=-999
	      K(25,1)=0
	      K(26,1)=0
	   Endif
      ENDIF

C...Deactivate obsolete lines 9, 10, 21, NS+1 (extra lines with boson)
      K(21,1)=0
      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
      K(13,1)=0
c do not delete lines, if we have jpsi decayed into leptons already
      if(ipsiel1.ge.1) then
        else
        K(14,1)=0
        K(15,1)=0
      endif
C...Zero irrelevant lines with K(I,1)<0
      DO 150 I=1,N

         IF(K(I,1).LE.0) THEN
            DO 140 J=1,5
               K(I,J)=0
  140       P(I,J)=0.D0
         ENDIF
	   
  150 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' caps: before PYEDIT '
         CALL PYLIST(1)
      ENDIF
C...Delete internal parton lines, i.e. with K(I,1)=13,14
C...Delete lines from propagator gluons in initial state ps
      IF(IPST.EQ.0)  CALL pyedit(14)

      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' caps:after LUEDIT(14) '
         CALL PYLIST(1)
      ENDIF
C...Delete empty lines
      CALL PYEDIT(12)
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' caps:after LUEDIT(12) '
         CALL PYLIST(1)
      ENDIF
	if(ilha.le.100.and.Iremn.eq.1) then 
c check for enegy momentum conservation
         DO 170 I=1,N
            IF(P(I,5).LT.0.0D0.OR.K(I,1).GT.2) GOTO 170
            DO 160 J=1,4
               IF(ABS(P(I,J)).LE.1.E-6) THEN
                  P(I,J)=0.0D0
               ENDIF
  160       CONTINUE
  170    CONTINUE
         CHEC = 5D-4
         CHEC = 1D-2
         P1 =ABS(P(1,1)+P(2,1))
         P2 =ABS(P(1,2)+P(2,2))
         P3 =ABS(P(1,3)+P(2,3))
         P4 =ABS(P(1,4)+P(2,4))
         PC1 = (ABS(PYP(0,1)) - P1)/P4
         PC2 = (ABS(PYP(0,2)) - P2)/P4
         PC3 = (ABS(PYP(0,3)) - P3)/P4
         PC4 = (ABS(PYP(0,4)) - P4)/P4
         IF(DABS(PC1).GT.CHEC.
     +      OR.DABS(PC2).GT.CHEC.
     +      OR.DABS(PC3).GT.CHEC.
     +      OR.DABS(PC4).GT.CHEC) THEN
            write(6,*) ' caps: energy of final particles not correct'
     +      ,chec,LST(21)
            write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
            write(6,*) ' x1 = ',pyvar(31),' x2 = ',pyvar(32)
            write(6,*) ' xg1 = ',xg1,'xg2_bar',xg2_bar
            write(6,*) ' xg2 = ',xg2,'xg1_bar',xg1_bar

            write(6,*) ' nia1 = ',nia1,' nia2 = ',nia2

            IF(IRES(1).EQ.0) THEN
               CALL PYROBO(0,0,STHETA(2),SPHI(2),0.D0,0.D0,0.D0)
               CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),
     +  		   DBETA(2,3))
            ENDIF
            do  i=1,n
               p(i,3) = - p(i,3)
            enddo
            CALL PYROBO(0,0,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),
     + 		 CM(3)/CM(4))
            call PYLIST(1)
            LST(21) = 100
            return
         ENDIF
      ENDIF
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' in cm system '
         CALL PYLIST(1)
      ENDIF
      IF(IRES(1).EQ.0.OR.(IRES(1).EQ.1.and.IABS(K(1,2)).EQ.11)) THEN
	   IF(IRES(1).EQ.1.and.IABS(K(1,2)).EQ.11) then
	   endif
         CALL PYROBO(0,0,STHETA(2),0d0,0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0d0,SPHI(2),0.D0,0.D0,0.D0)
         CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
      ENDIF
c... change z coordinate  to CASCADE standard
      do 180 i=1,n
  180 p(i,3) = - p(i,3)
      DO 190 I=1,4
         IF(NIEL.GT.0) then
            IF(ABS(P(NIEL,I)).LE.1.E-6) THEN
               P(NIEL,I)=0.0D0
            ENDIF
         ENDIF
         P(1,I)=PBEAM(1,I)
         P(2,I)=PBEAM(2,I)
  190 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' after boost dbeta in caps'
         CALL PYLIST(1)
      ENDIF
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' caps final '
         CALL PYLIST(1)
      ENDIF
c gg veto condition
      if(iveto.eq.1) then 
         LST(21) = 90
      endif
c end gg veto
      RETURN
      END
