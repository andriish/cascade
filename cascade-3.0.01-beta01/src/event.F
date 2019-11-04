      SUBROUTINE EVENT
      IMPLICIT NONE
      Integer NDIM
      COMMON/INTDIM/ NDIM
#include "caseffic.inc"
#include "caspara.inc"
#include "caspara1.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "casludat1.inc"
#include "casgki.inc"
#include "caslujets.inc"
#include "casvalues.inc"
#include "caspara5.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspar.inc"

      Integer Iremn
      common/casflags/Iremn
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
	Integer MSTP,MSTI
	Double Precision PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

      Double Precision X(20)
      Double Precision  XGR(20)
      Integer IMIX
      COMMON /OALPHAS/ IMIX

      Double Precision XV
      Integer NDIMEN
      COMMON /XVAL/ XV(20),NDIMEN
      Double Precision dcasrn
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI

      Integer Nqqb,Nccb,Nbbb,Npr1,Npr2,Npr3
      Common/NEVOUT/Nqqb,Nccb,Nbbb,Npr1,Npr2,Npr3

      Integer LERR
      COMMON/ERR/LERR(1000)
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      Double Precision remPARJ32
      LOGICAL FIRST,first2

      Integer NPRIN,IDEBUGO,IPROO,KPAO,KINT22,MSTP128
      Integer IHFO
      Integer I,J
      Double Precision WEIGL,WEIGC,WEIGB,WEIGT
      COMMON /weiglh/weigl,weigc,weigb,weigt
      DOUBLE PRECISION WeigW1,WeigW2,WeigW3,WeigW4,WeigW5,WeigW6
      COMMON /weigW/WeigW1,WeigW2,WeigW3,WeigW4,WeigW5,WeigW6

      Double Precision WTOT,Xsec
      Integer IQUARK
      Double Precision p1,p2,p3,p4,pc1,pc2,pc3,pc4,chec
      Integer KWZ,KWZP1,KWZP2
      Common/CASWZ/KWZ,KWZP1,KWZP2
      
      Double Precision dtest
      EXTERNAL dcasrn
      DATA FIRST/.TRUE./
      DATA FIRST2/.TRUE./
      DATA NPRIN/5/
c      DATA NPRIN/100/
      SAVE remPARJ32, IHFO, IPROO, KPAO, KINT22   
      
      IDEBUGO = IDEBUG
      MSTU(22)=100
c information on resonances in event record: 
c MSTP(128)>0 destroys strategy in caupevnt      
      MSTP128 = MSTP(128)
      MSTP(128) = 0
	IF(IDEBUG.EQ.1) write(6,*) ' in event '
   10 CONTINUE
      NIN=NIN+1
c      IF(MOD(NOUT,1).EQ.0) write(6,*) ' start: Nr of events ',NOUT

      IF(FIRST) THEN
         first = .FALSE.
         remPARJ32 = PARJ(32)
         IPROO = IPRO
         KPAO = KPA
         KINT22 = KINT(2,2)
         FIRST=.FALSE.
         Nqqb = 0
         Nccb = 0
         Nbbb = 0
         Npr1 = 0
         Npr2 = 0
         Npr3 = 0
         IHFO = IHF
         IF(ILHA.GE.1.or.ILHA.eq.-1) then
 	      call caupinit
 	   endif
      ELSE
         PARJ(32)=remPARJ32
         IHF = IHFO
         IPRO = IPROO
         KPA = KPAO
         KINT(2,2) = KINT22
      ENDIF
      PT2H = 0.d0
      SHH = 0.d0
      IWEI = 0
      IGENFL = 0
      IMIX = 0
      
      if(IPRO.ne.-1) then
         DO 20 I=1,20
            XV(I)= 0.D0
   20    X(I) = 0.D0
         NDIMEN = NDIM
         CALL casgen(NDIM,XGR)
         DO 30 I=1,NDIM
            X(I) = XGR(I)
   30    XV(I) = X(I)
         IF(IDEBUG.EQ.1) write(6,*) ' events are generated ',IPRO
         IMEPHI = 1
            
         IF(IPRO.GE.10.and.IPRO.LT.20) THEN
c select light or heavy quark production
            Call Xsec1(X,Xsec)
         elseif(ipro.eq.2) THEN
            call Xsecpsi(X,Xsec)
         elseif(ipro.eq.3) then
            Call Xsecchi(X,Xsec)
         elseif(ipro.eq.4) then
            Call Xsecqg(X,Xsec)
         elseif(ipro.eq.21) then
            Call Xsecjj(X,Xsec)
         elseif(ipro.eq.102) then
            Call Xhiggs(X,Xsec)
         elseif(ipro.eq.500) then
            Call Xv1(X,Xsec)
         elseif(ipro.eq.501) then
            Call XVq(X,Xsec)
         elseif(ipro.ge.503.and.ipro.le.509) then
            Call XVQQ(X,Xsec)
         elseif(ipro.eq.510) then
            write(6,*) ' event: process not yet implemented: STOP '
         elseif(ipro.eq.511) then
            WTOT = WEIGW1+WEIGW2
            IF(dcasrn().GT.WEIGW1/WTOT) THEN
              KWZ = -24
            else
              KWZ = 24
            Endif
c            write(6,*) KWZ
            Call XVq(X,Xsec)
         elseif(ipro.ge.513.and.ipro.le.519) then
            Call XVQQ(X,Xsec)
         Endif

         IF(NFLAV.GT.3.AND.IPRO.EQ.10) THEN
            WTOT = WEIGL + WEIGC + WEIGB + WEIGT
            IPRO = 10
            Call Xsec1(X,Xsec)
            dtest = dcasrn()
            IF(dtest.GT.WEIGL/WTOT) THEN
               IPRO = 11
               IHFLA = 4
               IF(dtest.GT.(WeigL+WeigC)/Wtot) IHFLA = 5
               IF(dtest.GT.(WeigL+WeigC+WeigB)/Wtot) IHFLA = 6
               Call Xsec1(X,Xsec)
            ENDIF
         ENDIF
         if(xsec.le.0 ) goto 10
	   IF(IDEBUG.EQ.1) write(6,*) ' after event mixing etc 1st '
         IMEPHI = 0

         IDEBUG = IDEBUGO
c check whether event mixing wanted by user
         IGENFL = 1
c............
      
      Endif
      
      IWEI = 0
      IHF=0

c      call pylist(1)
            
      IF((KBEAM(1,2).eq.22.or.Iabs(KBEAM(1,2)).eq.11).and.
     &  (KBEAM(2,2).eq.22.or.Iabs(KBEAM(2,2)).eq.11)) Then
	   Call capsee
	   Else
	   if(Mult.eq.0) then
	      Call caps
c            write(6,*) ' after caps '
            
		else
		call caps_ms
	   Endif
	Endif

      IHF = IHFO
c check on caps and pyremn
	IF(IDEBUG.EQ.1) write(6,*) ' after caps ',lst(21)
      IF(LST(21).NE.0) THEN
c         write(6,*) LST(21) 
         LERR(LST(21)) = LERR(LST(21)) + 1
         If(ipro.ne.-1) then
           GOTO 10
           else
           return
         endif
      ENDIF
	IF(IDEBUG.EQ.1) write(6,*) ' after event mixing etc '
c change interacting gamma to z0 on request from experiment
c      DO 40  I=3,N
c         IF(K(I,2).EQ.22.AND.K(I,1).EQ.21) K(I,2)=23
c   40 CONTINUE
      MSTP(128) = MSTP128
c      write(6,*) 'before pyexe '
c      call pylist(1)
      IF(NFRAG.GE.1.and.ILHA.le.0.and.Iremn.eq.1) then
        if(p6frag.eq.1) then 
           CALL PYEXEC
        else
c           call pyedit(11)
c           call pyedit(12)
c           call pyedit(13)
        call pyprep(0)
        call pylist(2)
c       do II =1,N
cc       write(6,*) ' I , KC, KCC, ',II,PYCOMP(K(II,2)),KCHG(PYCOMP(K(II,2)),2)
c        enddo
           call pythia8frag
        endif
c      write(6,*) 'after pyexe '
      Endif 
	
      IF(ILHA.GE.1) then
	   call caupevnt
	endif

      CALL PYROBO(0,N,0.d0,0.d0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
	
	
c....
      IF(MSTU(24).NE.0) THEN
         WRITE(6,*) 'MSTU(24)= ',MSTU(24)
         WRITE(6,*) ' Event rejected '
         CALL PYLIST(1)
         If(ipro.ne.-1) GOTO 10
         LERR(91) = LERR(91) + 1
         return
      ENDIF

      IF(ILHA.le.0.and.Iremn.eq.1) then
         CHEC = 5.D-2
c         CHEC = 5.D-5
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
            write(6,*) ' event: energy of final particles not correct'
     +      ,chec
            write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
            write(6,*) ' Nr of events ',NOUT
            write(6,*) ' N = ',N
            Call PYLIST(1)
         ENDIF
      ENDIF

      NOUT = NOUT + 1
      IF(MOD(NOUT,1000).EQ.0) write(6,*) ' Nr of events ',NOUT
      if(ILHA.eq.1) then 
      else
c      write(6,*) ' event: IPRO = ',IPRO,KPA
      Iquark = Iabs(K(NF1,2))
      If(IQUARK.LE.3) THEN
         Nqqb = Nqqb +1
      Elseif(IQUARK.EQ.4) THEN
         Nccb = Nccb + 1
      Elseif(IQUARK.EQ.5) THEN
         Nbbb = Nbbb + 1
      Endif
      If(Irespro.eq.1) then
         Npr1 = Npr1 + 1
      Elseif(Irespro.eq.2) then
         Npr2 = Npr2 + 1
      Elseif(Irespro.eq.3) then
         Npr3 = Npr3 + 1
      Endif
c      write(6,*) ' event ',nout,nprin,ilha
      IF(NOUT.LE.NPRIN.and.ILHA.le.0) THEN
         write(6,*) ' event: IPRO = ',IPRO
         CALL PYLIST(1)
      ENDIF
      endif
      RETURN
      END
