      SUBROUTINE p_mscatt(X,WPART)
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
      Double Precision X(20)
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
	Double Precision  QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision QF,W12
      DATA QF/0.0D0/,W12/0.0D0/
      Integer NDIMEN
      Integer I,J,IN,NP,NPP,KPA2,NRN,IST,NPFIN,KI
      Integer nb2,kpf,kstart,kpao,kpao2,nacc
      Double Precision WPART,W02,QFT,YX,QT2
      Double Precision QG2MAX,QG2MIN,WQG22,SMIN
      Double Precision SHAT1,PHI,FGAM,FWEI,THETE,PR,CPHP,SPHP
c      Double Precision STHP,CTHP,YGAM,YGAM_BAR
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
      
	double precision q2m

      NDIMEN = NDIM
c      DO 10  IN=1,20
c         K(IN,1) = 0
c         K(IN,2) = 0
c   10 CONTINUE

C...  GIVE BEAM  FOUR VECTORS
c      DO 20 I=1,2
c         DO 20 J=1,5
c            K(I,J)=KBEAM(I,J)
c   20 P(I,J) = PBEAM(I,J)
      ME =P(1,5)
      MP =P(2,5)
      NP=2
      NPP=3*NP
c      N=2
c      IDEBUG = 1
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
         AM(2) = PYMASS(KPA)
      ELSEIF(IPRO.EQ.2) THEN
         AM(1) = PYMASS(Kpsi)
         AM(2) = 0.d0
      ELSE
         WRITE(6,*) 'wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
c	write(6,*) ' in p_mscatt s=',sss,N,IPRO
	Fwei1 = 1d0
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
      W02=(AM(1)+AM(2)+MP)**2
c      IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
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
c            KPA=4
c            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
            IF(dcasrn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
         AM(1)=PYMASS(KPA)
         AM(2)=PYMASS(KPA2)
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
	Q2M = 10000.
c      IF(IRES(1).EQ.0) THEN

	IF(Iabs(KE).eq.2212) then
c only p-p included
         NRN = NRN + 1
         YX = YMIN*((YMAX/YMIN)**X(NRN))
c	   write(6,*) 'p_mscatt ',yx,Ymin,ymax,x(nrn),nrn
         XG1=YX
	   XG1SF = XG1
         QG2MIN = 0.0001d0
         Q2MIN = QG2MIN
         Q2MAX=YX*(SSS-MP**2) - W12
         IF(Q2MAX.LT.Q2MIN) GOTO 100
	   If(Q2MAX.GT.Q2M) Q2MAX = Q2M
         NRN = NRN + 1
         Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         XMIN = 0.0D0
         XMIN = 0.0D0
         qg22 = 0.D0
         Wqg22 = 1.D0
         QT2 = Q2*(1.d0-yx)
         QG2MAX = SSS*YX - Q2 - 4.D0 * AM(1)**2
c            IF(ICCFM.EQ.0) QG2MAX = DMIN1(Q2+4D0*AM(1)**2,QG2MAX)
	   If(QG2MAX.GT.Q2M) QG2MAX = Q2M
         QG2MAX=DMAX1(QG2MIN,QG2MAX)
         NRN = NRN + 1
         XZOHRN = X(NRN)
c 1/qt**2
         qg22 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
         Wqg22 = qg22*DLOG(QG2MAX/QG2MIN)
c        write(6,*) ' p_mscatt Q2,qg22: ',Q2,qg22
         SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
         XMIN=(SMIN+QT2+qg22)/(YX*SSS)
C ... xg2 = E_parton /E_proton
         NRN = NRN + 1
         xg2= XMIN*((XMAX/XMIN)**X(NRN))
c	   write(6,*) ' xg2 ',nrn
	   XG2SF = XG2
      Else
	   Write(6,*) ' IRES(1)=1 for KE = ',KE,' not implemented '
         Stop
	Endif


      IF(xg2.GE.0.9999.OR.xg2.LT.XMIN) GOTO 120
      SHAT1 = SSS*XG1*xg2
ctest      SHAT1=SSS*XG1*xg2 - Q2 - qg22
ctest      IF(ISEMIH.EQ.1)  SHAT1=SSS*XG1*xg2 - Q2 - qg22/(1-xg2)
      IF(SHAT1.LT.(AM(1)+AM(2))**2) GOTO 130

      IST = 0
      YY=YX
      XG=xg2
      PHI = 99999.D0
      NRN=NRN+1
      PHI = 2.D0*PI*X(NRN)
      CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      XP=0.D0
      IRES(2)=1

C FINAL STATE PROTON
      NPFIN=NIA1+4
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      K(NPFIN,3)=2
      P(NPFIN,5)=PYMASS(KP)
      PR = P(2,3)*(1.D0 - xg2)
      P(NPFIN,4) = PR
      CPHP=1.D0
      SPHP=DSQRT(1.D0 - CPHP**2)
      IF(ISEMIH.EQ.1) THEN
         STHP = DSQRT(qg22)/PR
      ENDIF
      IF(DABS(STHP).GT.1.D0) goto 190
      CTHP=DSQRT(1.D0 - STHP**2)
      P(NPFIN,1)= PR*STHP*CPHP
      P(NPFIN,2)= PR*STHP*SPHP
      P(NPFIN,3)= PR*CTHP
c      write(6,*) ' lin 372 '
c      CALL PYLIST(1)
C MOMENTA OF PARTON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KGL
      K(NIA1+1,3)=2
      DO 40  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   40 CONTINUE
      IF(ISEMIH.EQ.0) THEN
         P(NIA1+1,4)=ABS(P(NIA1+1,3))
         P(NIA1+1,5)=0.0D0
      ELSE
         P(NIA1+1,4)=DSQRT(P(NIA1+1,1)**2+P(NIA1+1,2)**2+P(NIA1+1,3)**2
     +               - qg22)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
c here define new momenta: introduce maximum virtuality mu**2 = Q2 + 4M**2
c for the matrix element and kinematics we set Virt_max=kt**2
c     define lightcone variable xg=xg2
c         XG = xg2
c     define lightcone variable xg_bar=(-kt**2)/(1-xg)/SSS
         XG2_bar = (-qg22)/(1.d0-xg2)/SSS
c here check on virt. gamma
         QT2 = PYP(NIA1,9)
         IF(XG1.NE.1D0) THEN
            XG1_BAR = (-QT2)/(1.d0 - XG1)/SSS
         ELSE
            XG1_BAR=0D0
         ENDIF
         KT2_1 = ABS(QT2)
c         xsi_hard=(ygam+xg_bar)/(xg+ygam_bar)
         xsi_hard=(xg1+xg2_bar)/(xg2+xg1_bar)

	   xsim = xsi_hard
         ykt = sqrt(abs(qg22))
         YY=YX
         if(xsi_hard.le.0) goto 200
         P(NIA1+1,4) =  xg2_bar*P(1,4) + xg2*P(2,4)
         P(NIA1+1,3) =  xg2_bar*P(1,3) + xg2*P(2,3)
         NRN = NRN + 1
         PHI_K2 = 2*PI*X(NRN)
         P(NIA1+1,1) = DSQRT(qg22)*SIN(PHI_K2)
         P(NIA1+1,2) = DSQRT(qg22)*COS(PHI_K2)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
         KT2_2=ABS(qg22)
c now we also change the photon: set plus component to zero
c         P(NIA1,4) =  ygam*P(1,4) + ygam_bar*P(2,4)
c         P(NIA1,3) =  ygam*P(1,3) + ygam_bar*P(2,3)
         P(NIA1,4) =  xg1*P(1,4) + xg1_bar*P(2,4)
         P(NIA1,3) =  xg1*P(1,3) + xg1_bar*P(2,3)
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
c        write(6,*) ' p_mscatt xg: ',xg
         P(2,4)=dabs(p(2,3))
         DO 50   KI=1,4
            P(NPFIN,KI)=P(2,KI)-P(NIA1+1,KI)
   50    CONTINUE
      ENDIF
      NIA2 = NIA1+1
c      call PYLIST(1)
      NF1=NIA1+2
      NF2=NIA1+3
      K(NF1,1)=2

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
      NB2 = NIA2
c      write(6,*) ' kp2,kpf ',kpa2,kpf
C...   VECTOR OF GAMMA GLUON CMS SYSTEM

      DBCMS(1)=  P(NIA1,1) + P(NB2,1)
      DBCMS(2)=  P(NIA1,2) + P(NB2,2)
      DBCMS(3)=  P(NIA1,3) + P(NB2,3)
      DBCMS(4)=  P(NIA1,4) + P(NB2,4)
      DO 60  I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
   60 CONTINUE
      SHAT=DOT(DBCMS,DBCMS)
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
            KPAO = KPA
            KPAO2 = KPA2
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
c         write(6,*) ' IPRO = 13 ,NRN= ',NRN
      NRN = NRN + 1
      XY(1) = X(NRN)
      NRN = NRN + 1
      XY(2) = X(NRN)
c         write(6,*) am(1),am(2),ipro
c         write(6,*) 'p_mscatt ',XY(1),XY(2),NDIM,NDIMEN,IMIX,IWEI
c      write(6,*) ' 1st p_mscatt NRN ',NRN,N
c      write(6,*) ' 1st p_mscatt Nia1,2,nf1,2 ',Nia1,Nia2,nf1,nf2,npfin
	CALL PHASE(NP,ECM,AM,PCM,WT)
	
      IF(WT.LE.0.D0) GOTO 180
      DO 80  I=1,4
         P(NF1,I)=PCM(I,1)
         P(NF2,I)=PCM(I,2)
   80 CONTINUE
c      write(6,*) ' p_mscatt ', KPA,KPA2,K(NF1,2),K(NF2,2),NRN,X(NRN)
      PT2 = PYP(NF1,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) GOTO 191

C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
      call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      XG = xg2
      if(xg2.GT.1.d0) write(6,*) 'p_mscatt xg2>1',xg2,ipro
c      call pylist(1)
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
      XP = 1.d0/xg2
      XPINT=DLOG(XMAX/XMIN)
      WTGLU = xg2 * XP * XPINT
      WTGLU = WTGLU*Wqg22
      WMAT=WT *WTGLU*FGAM*FWEI * Fwei1
c	write(6,*) FGAM,FWEI,FWEI1
      WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
      FLUX=(2.D0*PI)**4
      FLUX=FLUX/(2.D0*xg2*xg1*sss)
      WPART = WMAT *FLUX

      PT2H = PT2
      SHH = SHAT


      NDIMC = NRN
c      write(6,*) ' 2nd p_mscatt NRN ',NRN,N
c      write(6,*) ' 2nd p_mscatt Nia1,2,nf1,2 ',Nia1,Nia2,nf1,nf2,npfin
c	call pylist(1)
      CALL PYEDIT(14)
c      write(6,*) ' 3rd p_mscatt NRN ',NRN,N
c      write(6,*) ' 3rd p_mscatt Nia1,2,nf1,2 ',Nia1,Nia2,nf1,nf2,npfin
c      write(6,*) ' end of p_mscatt ',NDIMC,NRN
c      CALL PYLIST(1)
c      write(6,*)  'p_mscatt: WPART ',WPART
      KPAO = KPA
      KPAO2 = KPA2
	
      CALL PYROBO(0,N,0.D0,0.D0, DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +DBCMS(3)/DBCMS(4))
	
c	write(6,*) ' end p_mscatt N= ',N
c	call pylist(1)
      if(wpart.ne.wpart) write(6,*) ' p_mscatt wpart = ',wpart,xsi_hard
c      write(6,*) ' p_mscatt end ,kpa,kpa1,kpa2,kpao2',kpa,kpa1,kpa2,kpao2
      RETURN
   90 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      RETURN
  100 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
         NDIMC = 9999
      ENDIF
      RETURN
  110 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: xlimit ; RETURN xmin,xmax',xmin,xmax
      ENDIF
      NDIMC = 9999
      RETURN
  120 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: xg2imit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  130 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: slimit ; RETURN ',shat1,am(1),am(2),IPRO
         write(6,*) ' sss,xg1,xg2,q2,qg22 ',sss,xg1,xg2,q2,qg22
         write(6,*) ' p_mscatt q2,y,sss',q2,yx,sss
         write(6,*) ' p_mscatt 2*pq ',2*dot1(nia1,2),sss*yx
      ENDIF
      NDIMC = 9999
      RETURN
  140 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: theta limit ; RETURN ',themi,thema,thete
      ENDIF
      NDIMC = 9999
      RETURN
  150 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: bochck limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  160 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: shat ; RETURN ',shat,shat1,kpa
c        write(6,*) ' p_mscatt: QG2MAX,qg22 ',QG2MAX,qg22
c        write(6,*) ' p_mscatt: yx,xg2,q2,s ',yx,xg2,q2,sss
c        write(6,*) ' p_mscatt 2*pq ',2*dot1(nia1,2),sss*yx
c        write(6,*) ' p_mscatt nia1,nia2 ',nia1,nia2
c        write(6,*) ' p_mscatt 2*ktgam*ktglu ',2*(P(Nia1,1)*P(nia2,1)+
c     &  P(Nia1,2)*P(NIA2,2))
      ENDIF
      NDIMC = 9999
      RETURN
  170 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: ECM limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  180 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: PHASE WT=0 ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  191 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: PTCUT limit ; RETURN ',PT2,KPA
      ENDIF
      NDIMC = 9999
      RETURN

  190 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: qg22 limit; RETURN ',DSQRT(qg22),PR,STHP
      ENDIF
      NDIMC = 9999
      RETURN
  200 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' p_mscatt: xsi_hard limit ; RETURN ',xsi_hard
      ENDIF
      NDIMC = 9999
      RETURN


      END
