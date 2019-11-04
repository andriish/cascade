      SUBROUTINE cascps(IPU1,IPU2)
C initial state backward parton shower
C with angular ordering
C
C
C
C
C              |
C       XB     |
C    (IFLB)    |
C              +-------
C              |
C    XA= XB/Z  |
C     (IFLA)   |
C
C
      Implicit None
C...GENERATES SPACELIKE PARTON SHOWERS
      Double Precision XLP,YLP,Q2LP
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /CASTOU/PARL(30),XLP,YLP,Q2LP,LLST(40)
#include "caslujets.inc"
#include "casludat1.inc"
#include "casluco.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casshwr.inc"
#include "casgki.inc"
#include "caspar.inc"
#include "casherup.inc"
#include "caslhe.inc"

      Double Precision XA,XB
      Double Precision TEVS(2)
      DOUBLE PRECISION KT2,KT2B
      DIMENSION KT2(2)
      DOUBLE PRECISION THE2(2)
      Double Precision dcasrn
      EXTERNAL dcasrn
#include "caspara5.inc"
      LOGICAL LPRINT,LSPRINT,LPORDER
      LOGICAl CCFM
      DOUBLE PRECISION DETOT,DBETA,STHETA,SPHI
      COMMON /MYBOOST/ DETOT,DBETA(2,3),STHETA(2),SPHI(2)
      DOUBLE PRECISION QG21,QG22
      COMMON/SEMIH/ QG21,QG22
      Integer IS,IFLS,KFL,JR,IPO,IFLB
      DIMENSION IS(2),IFLS(4)
      Double Precision X,SH
      COMMON /MYPROC/KFL(3,2),X(2),SH

      double precision QT2P(2),XS(2),XS_BAR(2)
      double precision k2d
      double precision kxg(2,1000),kxg_bar(2,1000)
      double precision QK(2,1000),QKT(2,1000)
      double precision QQT(2,1000),QQTP(2,1000),TH(2,1000),ZZ(2,1000)
      Double Precision PSFIN(4)
      INTEGER i_ng,NPROP(2,1000)
      Double Precision k2,kt2bs
      common /plot/ k2,kt2bs
      Double precision kx,ky
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	
      Integer i,j,it,jt,Nprint,ncall,ifla,nq,np,ntry,isemiho,ilep,ns
      Integer i1,i2,itemp,ipu1,ipu2,iflag,kn1,kd1,id1,id2,jtfin,npsfin
      Double Precision z,q0
      Double Precision pxsum,pysum,xktsum,yktsum,qt2
      Double Precision yph,yph_bar,xsi_hard(2),qt2max,qt2cut,tmax
      Double Precision xsi,qt2_glu,t_maxg,t_hard,xsi_old
      Double Precision xsim,zmax,dx,dxp,pxt,pyt,tevb,yi,qt2b
      Double Precision DMAS_MAX,xgbar,theta,q2_test,tnew
      Double Precision xsi_test
      Double Precision pp_a0
      Double precision xgc,x0,Qs
      Double Precision the_max,DOT1
      Integer MSTJO
      common /cbran/ xgc,x0,Qs
	Integer kf
	Common /cbran1/ kf
	Double Precision  BB,ktqt,Qcdlam,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlam,scal,qtd,ktd

      Double Precision PZM,PEM,PZF,DBETAZ
      Double Precision PZMB,PEMB,PZFB
      Double Precision PB_SAV(2,5),dboost,dboostl
      Double Precision PSUMF(4),cmtest,dot

	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar


      Double Precision mtest
c      data mtest/5.d0/
      data mtest/1.d0/

      Logical reshuff
c      Data reshuff/.false./
      Data reshuff/.true./ ! new default
      Integer Iresl1,Iresl2
      Double Precision en,rap_lab,pt_lab,ptot,et_lab
      common /hzgrap/rap_lab,et_lab
c      DATA LPRINT/.TRUE./
      DATA LPRINT/.FALSE./
      DATA LSPRINT/.FALSE./
c      DATA LSPRINT/.TRUE./
c      DATA LPORDER/.TRUE./
      DATA LPORDER/.FALSE./
      DATA NPRINT /1000/
      DATA NCALL /0/
      DATA IFLA,NQ/0,0/,Z,XA/2*0./
      DATA IS/2*0/
c      Data pp_a0/0.54/
      Data pp_a0/0.1/
      DATA JTFIN/2/
c      Character*5 xxxx
c      Data xxxx/'bran'/

	Double Precision bitx,bity,bitz,Etest
	double precision pshw(4)
      double precision pxff,pxs,es,pyff,pys,pzs,pzff,betax,betay,betaz
      DOUBLE PRECISION PHITIM,THETIM
      
      Double Precision STMAX1,STMAX2
      common/myscale/STMAX1,STMAX2
      
      double precision shtest(4)


      logical first
      Data first/.true./
      
      integer itimtest
      Data itimtest/1/
      
      if(first) then
         first=.false.

      endif
      do j=1,4
      psumf(j) = 0.
      end do

      IF(IORDER.LT.1.OR.IORDER.GT.3) THEN
         Write(6,*) ' wrong ordering selected IORDER = ',IORDER
         Write(6,*) ' swtich to IORDER= 1 (angular ordering) '
         IORDER=1
      ENDIF
      Do I=1,4
         PSFIN(I)=0.D0
      Enddo
      NTRY = 0
      QCDLAM = 0.20D0
      QCDLAM = PARU(112)
      if(Qg0.le.0) QG0=1.
      Qg0 = 0.3 ! new hack for new TMDs and to prevent large remnant pt
      Q0 = QG0
	iresl1=0
	if(iabs(k(1,2)).eq.11.and.ipy(41).eq.22) iresl1=1
      iresl2=0
	if(iabs(k(2,2)).eq.11.and.ipy(42).eq.22) iresl2=1
      NCALL = NCALL + 1
      ISEMIHO = ISEMIH
      IF(ICCFM.EQ.0.or.ICCFM.EQ.3) THEN
         CCFM = .FALSE.
      ELSE
         CCFM = .TRUE.
      ENDIF
c      IF(LSPRINT) write(6,*) ' New event'
      IF(QG21.LE.0.AND.QG22.LE.0)   then
c         LST(21)=55
         IF(NCALL.lt.10) write(6,*) ' cascps: Kt < 0 : no ps possible'
         GOTO 60
      endif
      do jt=1,2
         do I=1,1000
            KXG(jt,I) = 0.d0
            KXG_BAR(jt,I) = 0.d0
            QK(jt,I) =0.d0
            QKT(jt,I)=0.d0
            QQT(jt,I)=0.d0
            QQTP(jt,I)=0.d0
            TH(jt,I)=0.d0
            ZZ(jt,I)=0.d0
            NPROP(jt,I) = 0
         ENDDO
      enddo

      IF(NCALL.GE.NPRINT) LPRINT=.FALSE.
C...COMMON CONSTANTS, SET UP INITIAL VALUES
c perform the boost into the overall cms (ep cms)
C...Transform to e p cms, boost parameters in double precision.
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      if(ILEP.eq.0.and.(iresl1.eq.1.or.Iresl2.eq.1)) then
         dboostl = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
         CALL PYROBO(5,0,0.d0,0.D0,0.D0,0.D0,dboostl)
      Endif	
      CALL PYROBO(0,0,STHETA(2),0.D0,0.D0,0.D0,0.D0)
      CALL PYROBO(0,0,0.D0,SPHI(2),0.D0,0.D0,0.D0)
      CALL PYROBO(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
c      write(6,*) ' in cascps start ',dboostl,PYVAR(31),PYVAR(32)
c      call pylist(1)
      
      DO I=1,5
         PB_SAV(1,I) = P(1,I)
         PB_SAV(2,I) = P(2,I)
         P(1,I) = PBEAM(1,I)
         P(2,I) = PBEAM(2,I)
      ENDDO
c      write(6,*) ' cascps ILHA ',ILHA
      If (ILHA.ge.0) then 
        P(1,3) = -P(1,3)
        P(2,3) = -P(2,3)
      endif
c now go to lab system
      CALL PYROBO(0,0,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4))
c      write(6,*) ' in cascps  beginning '
c      call pylist(1)
      
      
      Qs =Q0
      NQ = 21
      NP = 23
      xgbar = 0.d0
      xgc = xg
      QT2_glu = PYP(NP,9)
      xg_bar = -QT2_glu/(1.d0-xg)/SSS
      QT2 = PYP(NQ,9)
      yph = xg1
      yph_bar = -QT2/(1.d0 - yph)/SSS
      
     
      
c      write(6,*) 'cascps start: ',xg,xg_bar,qt2_glu,yph,yph_bar,qt2
c      write(6,*) 'cascps start: ',PYVAR(31),yph,PYVAR(32),xg

      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)	
	xsi_hard(2) = min(xsi_hard(2),5.d8)
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
	xsi_hard(1) = min(xsi_hard(1),5.d8)	

	IF(ILEP.EQ.1) then 
         QT2MAX = (xsi_hard(2)*SSS)*xg**2
C check the following: must apply also for gg 
         IF(QG22.LE.Q0**2)   then
           LST(21)=55
           IF(LSPRINT) write(6,*) ' cascps: QG2 < Q0 : no ps possible'
c           write(6,*) ' qg21,qg22 ',qg21,qg22,q0         
           xsi_test=0.5
           GOTO 50
         endif
	ELSEIF(ILEP.EQ.2) then 
         QT2MAX = (xsi_hard(1)*SSS)*yph**2
         IF(QG21.LE.Q0**2)   then
           LST(21)=55
           IF(LSPRINT) write(6,*) ' cascps: QG2 < Q0 : no ps possible'
c           write(6,*) ' qg21,qg22 ',qg21,qg22,q0         
           xsi_test=0.5
           GOTO 50
         endif
	ELSE
         QT2MAX = max((xsi_hard(2)*SSS)*xg**2,(xsi_hard(1)*SSS)*yph**2)
	ENDIF
c	write(6,*) ' cascps ILEP = ',ILEP,sqrt(QT2MAX),
c     &   sqrt(min((xsi_hard(2)*SSS)*xg**2,(xsi_hard(1)*SSS)*yph**2))
c	write(6,*) ' cascps stmax: ',stmax1,stmax2
      NPSFIN=N

c define z for first branching z = 1
      z = 1.d0

      IF(CCFM) THEN
         QT2CUT = 1.D22
c  taking (ala PYTHIA) masses of incoming and outgoing partons + pt         
cc         QT2CUT = (P(21,5)**2 + P(23,5)**2 + P(25,5)**2 + P(26,5)**2)/2.
cc     &            + (P(25,1)**2+ P(25,2)**2 + P(26,1)**2+ P(26,2)**2)/2.
c         write(6,*) ' cascps test 21,26 = ',P(21,5),P(26,5)
c         write(6,*) ' cascps test 23,25 = ',P(23,5),P(25,5)
c         write(6,*) ' cascps test 27 = ',P(27,5)
      ELSE
         QT2CUT = max((P(21,1)**2 + P(21,2)**2), (P(23,1)**2 + P(23,2)**2))
c         QT2CUT = 1.D22
c         write(6,*) ' qt2cut ',qt2cut 
       ENDIF
ccc      QT2CUT = max((P(21,1)**2 + P(21,2)**2), (P(23,1)**2 + P(23,2)**2))

c here use scale and cut from lhe file, if they exist
      if(scale_process.gt.0) then
         QT2MAX = scale_process**2
c         if(scalup_cas**2.le.qt2cut) 
c     &    write (6,*) ' cascps qt2cut,scalup ', sqrt(qt2cut), scalup_cas
         QT2CUT = scalup_cas**2
c         write (6,*) ' cascps qt2max,scale ', sqrt(qt2max),scale_process
c         write (6,*) ' cascps qt2cut,scalup ', sqrt(qt2cut), scalup_cas
      endif
c         write (6,*) ' cascps qt2cut,qt2max ', sqrt(qt2cut), sqrt(qt2max)
c      write(6,*) ' new event '
c      call pylist(1)
c      write(6,*) 'new ', Qt2max 

cc 

      IF(LSPRINT) write(6,*) ' initial QT2MAX ',QT2MAX
      IF(Lporder) write(6,10000) sqrt(QT2MAX),sqrt(qg22)
10000 format(' max qt"    ',F8.3,' gen. kt = ',F8.3)
      TMAX = DLOG(QT2MAX/QCDLAM**2)
      
c      write(6,*) ' initial QT2MAX ',QT2MAX,sqrt(QT2MAX),QT2CUT,CCFM
c      write(6,*) 'cascps ',xg,xg_bar,qt2_glu,yph,yph_bar,qt2
c      write(6,*) ' start shower '
      IF(LPRINT) THEN
         write(6,*) ' hard scattering limits '
         write(6,*) 'p.(P+P) ',DOT1(2,27),'e.(P+P) ',DOT1(1,27)
         write(6,*) ' xsi_hard  ',xsi_hard
         write(6,*) ' QG2 = ',QG21,QG22
         write(6,*) ' tmax ',tmax
         write(6,*) ' cascps TMAX ',TMAX
         call PYLIST(1)
      ENDIF
      IF(QT2MAX.LE.pp_a0**2) THEN
         IF(LSPRINT) then
            write(6,*) ' return : tmax ',tmax
         endif
         IF(LPRINT) THEN
            write(6,*) ' cascps 1st check: '
            write(6,*) ' QCDLAM ',QCDLAM
            write(6,*) ' Q2LP ',Q2LP
         ENDIF
c         LST(21)=55
c         write(6,*) ' LST(21) = 55 ... goto 50 '
         GOTO 50
      ENDIF
      IF(LPRINT) write(6,*) 'cascps:check ok '
      NS=N
      MSTU(2)=0
      N=NS
      DO  I=N+1,N+8
         DO J=1,5
            K(I,J)=0
            P(I,J)=0.D0
         ENDDO
      ENDDO
      DO   JT=1,2
         IFLS(JT)=KFL(2,JT)
         IF(KFL(2,JT).EQ.21) IFLS(JT)=0
         IFLS(JT+2)=IFLS(JT)
c define the sudakov stuff for exchanged parton
         IF(JT.EQ.1) I2=NQ
         IF(JT.EQ.2) I2=NP
          KT2(JT) =  P(I2,1)**2 + P(I2,2)**2
         IF(JT.EQ.1) Then
            XS(JT)=yph
            XS_BAR(JT)=yph_bar
            xgc = yph
         ELSEIF(JT.EQ.2) THEN
            XS(JT)=xg
            XS_BAR(JT)=xg_bar
            xgc = xg
         ENDIF
c test hannes
         
         TEVS(JT)=TMAX
         
      ENDDO

      THE2(1)=180.D0
      THE2(2)=0.D0

      PEM=P(21,4)+P(23,4)
      PZM=P(21,3)+P(23,3)
      
     
      NPROP(1,1) = 21
      NPROP(2,1) = 23
cc new
      IS(1)=21
      IS(2)=23
cc

* loop over the legs
      JT=0
   10 CONTINUE
      JT = JT + 1
      IF(ILEP.EQ.1) JT=2
      PXSUM = 0.d0
      PYSUM = 0.d0
      XKTSUM = 0.d0
      YKTSUM = 0.d0
ccccccc      
ccc	JT=JTFIN ! test only 
ccccccc

      IFLB=IFLS(JT)
cccc      IFLA=IFLS(JT) ! hannes hack
      XSI = xsi_hard(jt)
      I_NG=1

      kf = iabs(IPY(40+JT)) 

c hannes hack for DIS: needs to be checked      
      if(kf.eq.11) then
c        write(6,*) ' cascps: KF = ',kf,JT
        goto 50
      endif
c end hannes hack for DIS
      
      IF(ILEP.EQ.0) THEN
c p p production in pp cms
         JR = JT
         I=N+1
         IPO=19+2*JR
         K(I,1)=14
         K(I,2)=IFLS(JR+2)
         IF(IFLS(JR+2).EQ.0) K(I,2)=21
         IF(IFLS(JR+2).EQ.22) K(I,1)=11
         K(I,3)=0
         K(I,4)=IPO
         K(I,5)=IPO
         IF(K(IPO,2).EQ.22) K(I,4)=0
         IF(K(IPO,2).EQ.22) K(I,5)=0
         If(JR.EQ.1) THEN
            DO J=1,5
               P(I,J) = P(NQ,J)
            ENDDO
         ELSE
            DO J=1,5
               P(I,J) = P(NP,J)
            ENDDO
         ENDIF

         K(I+1,1)=-1
         K(I+1,2)=K(IPO+1,2)
ccc         K(I+1,2)=IFLB ! hannes hack
         K(I+1,3)=I
         K(I+1,4)=0
         K(I+1,5)=0
         P(I+1,1)=0.D0
         P(I+1,2)=0.D0
         P(I+1,3)=IPO
         P(I+1,4)=IPO
         P(I+1,5)=0.D0
         P(IPO+1,1)=I
         P(IPO+1,2)=I
         K(IPO,4)=MOD(K(IPO,4),MSTU(5))+I*MSTU(5)
         K(IPO,5)=MOD(K(IPO,5),MSTU(5))+I*MSTU(5)
         IF(K(IPO,2).EQ.22) K(I,4)=0
         IF(K(IPO,2).EQ.22) K(I,5)=0
         IS(JR) = I
         I2 = I
         NPROP(JT,I_NG) = I2
         N=N+2
c      write(6,*) '  Parton branching: branching shower initiator iflb',IFLS(JT+2),K(IPO+1,2),I

      ELSE
C..LEPTOPRODUCTION EVENTS: in e p cms
         I1=NS+2*ILEP-1
         I2=NS-2*ILEP+5
         DO ITEMP=NS+1,NS+4
            DO J=1,5
               K(ITEMP,J)=0
               P(ITEMP,J)=0.D0
            ENDDO
         ENDDO
         DO J=1,5
            P(I1,J)=P(NQ,J)
            P(I2,J)=P(NP,J)
         ENDDO
         K(NS+1,1)=11
         K(NS+3,1)=14
         IF(ILEP.EQ.2) THEN
            K(NS+1,1)=14
            K(NS+3,1)=11
         ENDIF
         K(NS+2,1)=-1
         K(NS+4,1)=-1
         K(NS+1,3)=0
         K(NS+2,3)=NS+1
         K(NS+3,3)=0
         K(NS+4,3)=NS+3
         K(I1,2)=KFL(2,ILEP)
         K(I2,2)=KFL(2,3-ILEP)
         P(I2+1,3)=MAX(IPU1,IPU2)
         P(I2+1,4)=MAX(IPU1,IPU2)
         K(I2,4)=K(I2,4)-MOD(K(I2,4),MSTU(5))+MAX(IPU1,IPU2)
         K(I2,5)=K(I2,5)-MOD(K(I2,5),MSTU(5))+MAX(IPU1,IPU2)
         P(26-2*ILEP,1)=I2
         P(26-2*ILEP,2)=I2
         K(25-2*ILEP,4)=MOD(K(25-2*ILEP,4),MSTU(5))+I2*MSTU(5)
         K(25-2*ILEP,5)=MOD(K(25-2*ILEP,5),MSTU(5))+I2*MSTU(5)
         IS(2) = I2
         NPROP(JT,I_NG) = I2
         N=N+4
      ENDIF

c      write(6,*) ' start JT ',JT,' N = ',N



      IF(JT.EQ.1) THEN
         P(IS(1),3) = XS_BAR(1)*P(2,3) + XS(1)*P(1,3)
         P(IS(1),4) = XS_BAR(1)*P(2,4) + XS(1)*P(1,4)
      ELSEIF(JT.EQ.2) THEN
         P(IS(2),3) = XS(2)*P(2,3) + XS_BAR(2)*P(1,3)
         P(IS(2),4) = XS(2)*P(2,4) + XS_BAR(2)*P(1,4)
      ENDIF
      
c      write(6,*) ' check 1st ',is(jt),dot1(is(jt),is(jt))
      
      KXG(JT,I_NG) = XS(JT)
      KXG_BAR(JT,I_NG) = XS_BAR(JT)
c      write(6,*) ' 1st jt ',jt,i_ng,KXG(JT,I_NG),KXG_BAR(JT,I_NG)
   30 CONTINUE
      IF(N.GT.MSTU(4)-10) THEN
         WRITE(6,*) ' cascps: no more memory in LUJETS'
         LST(21)=51
         GOTO 50
      ENDIF
      DO   I=N+1,N+8
         DO  J=1,5
            K(I,J)=0
            P(I,J)=0.D0
         ENDDO
      ENDDO
      IF(LPRINT) THEN
         write(6,*) ' cascps JT,KT2(1),KT2(2) ',JT,KT2(1),KT2(2)
      ENDIF
c define xb as light cone variable

      XB=XS(JT)
      IFLB=IFLS(JT)
c      write(6,*) ' after 30 ',iflb,ifls(JT),JT
      KT2B=KT2(JT)
      TEVB=TEVS(JT)
      XSI_OLD = XSI
	
c here calculate z and qt' for next branching

      NTRY = NTRY+1
c      write(6,*) 'cascps ',Ntry
      IF(NTRY.GT.15000) THEN
         LST(21)=56
         IF(LSPRINT) write(6,*) ' NTRY > 15000 : RETURN'
         write(6,*) ' NTRY > 15000 : RETURN'
         GOTO 50
      ENDIF

      IF(CCFM) THEN
         T_MAXG = TEVB
      ELSE
c new upper limit: qt<qt_hard
         T_hard=log(QT2CUT/QCDLAM**2)
c new ordering criteria (one loop): qt'_i > qt'_(i-1)
c         write(6,*)  'cascps qt2cut check ',sqrt(qt2cut), sqrt(exp(Tevb))*qcdlam
         T_MAXG = MIN(TEVB,T_hard)
      ENDIF
c new
      zmax = 1.D0 - XB*Q0/SQRT(QCDLAM**2*EXP(T_MAXG))
      
c      write(6,*) ' selected t ',sqrt(exp(T_maxg))*qcdlam,' z_M = ', zmax, JT

cc      IF(XB.GE.ZMAX) THEN
cc         QT2B=0.D0
cc         IFLAG = 0
ccc         write(6,*) 'cascps: z_max =',zmax,' xb ',xb,q0,qcdlam,t_maxg,qt2cut
cc         GOTO 40
cc      ENDIF

      kx=P(IS(JT),1)
      ky=P(IS(JT),2)
      dx = xb


      tmax = T_MAXG
      xsim = xsi_old
	
      dxp = 0.d0
      pxt = 0.d0
      pyt = 0.d0
      z = 1.d0
      iflag=0
      if(dx.le.0.98) then
c      write(6,*) ' before casbran ',ifla,iflb,dx,tmax,kx,ky,pxt,pyt
      If(Iccfm.le.1) then 
          Call casbran(ifla,iflb,dx,kx,ky,tmax,xsim,dxp,pxt,pyt,z,tnew,sss,Iflag)
          elseif(Iccfm.ge.2) then
          Call casbran_new(ifla,iflb,dx,kx,ky,tmax,xsim,scalup_cas,dxp,pxt,pyt,z,tnew,sss,Iflag)
      endif 
      endif
      xa = dxp
c new t value
      qt2b = pxt**2 + pyt**2
      tevb = tnew
      
cc      if(qt2b.gt.qt2cut.and.iflag.ne.0) goto 30  ! hannes test
            
cc      write(6,*) ' after casbran ',dx,z,dx/z,sqrt(qt2b),sqrt(exp(Tevb))*qcdlam
c      write(6,*) ' after casbran: iflag = ', iflag,' ifla = ',ifla,'  iflb = ',iflb
c      if(iflb.eq.0) iflag =0
c      if(Iflag.ne.0) then 
c        write(6,*) ' after casbran: iflag = ', iflag,' ifla = ',ifla,'  iflb = ',iflb
c        write(6,*) ' after casbran: z = ',z,'  jt = ',jt, ' it = ',it
c      endif
      if(Iflag.eq.0) Then
         QT2B = 0.d0
      ELSE
         IF(LPRINT) THEN
            write(6,*) ' now we have z,xb,xa values:',z,xb,xa
            write(6,*) 'before HARD scatter in cm frame ',N,NS+4-2*
     +      MIN(1, ILEP)
         ENDIF
      ENDIF
      IF(iflag.ne.0) THEN

c
         N=N+2
c now do the branching with z and qt from the previous loop
         IF(LPRINT) THEN
            write(6,*) 'find max allowed mass: xa,xb,z ',xa,xb,z
         ENDIF
         IT=N-1
         K(IT,1)=3
         K(IT,2)=21

c         write(6,*) ' cascps check :',IFLB-IFLS(JT+2),IFLB,IFLS(JT+2)
c hannes hack
ccc         K(IT,2)=IFLA
ccc         if(IFLA.EQ.0) K(IT,2)=21
         IFLS(JT+2) = IFLA
         K(IT,2)=-IFLB+IFLS(JT+2)
         IF(IFLB-IFLS(JT+2).EQ.0) K(IT,2)=21
c hannes hack         
c         write(6,*)  ' Parton branching: branching time - like : ',IT ,K(IT,2),IFLB,IFLS(JT+2)
         P(IT,5)=PYMASS(K(IT,2))
         DO J=1,5
            K(IT+1,J)=0
            P(IT+1,J)=0.D0
         ENDDO
         K(IT+1,1)=-1
         K(IT+1,2)=K(IS(JT)+1,2)
         K(IT+1,3)=IT

C...RECONSTRUCT KINEMATICS OF BRANCHING: TIMELIKE PARTON SHOWER
         MSTU(2) = IT +1
         P(IT,1)=PXT
         P(IT,2)=PYT
         PXSUM = PXSUM + PXT
         PYSUM = PYSUM + PYT
         XKTSUM = PXT+P(IS(JT),1)
         YKTSUM = PYT+P(IS(JT),2)
         MSTU(2) = 0
         YI = (1.d0-Z)* XB/Z
         XSI =QT2B/SSS/YI**2
c         write(6,*) ' xsi ',xsi,qt2b,yi
c         if(xsi.gt.xsi_old.or.xsi.gt.xsi_hard(jt)) then
c            write(6,*) ' cascps: xsi = ',xsi,' xsi_old = ',xsi_old
c            write(6,*) ' cascps: xsi = ',xsi,' xsi_hard = ',xsi_hard(jt)
c         endif
         IF(JT.EQ.1) THEN
            P(IT,3) =YI*(P(1,3)+XSI*P(2,3))
            P(IT,4) =YI*(P(1,4)+XSI*P(2,4))
         ELSEIF(JT.EQ.2) THEN
            P(IT,3) =YI*(P(2,3)+XSI*P(1,3))
            P(IT,4) =YI*(P(2,4)+XSI*P(1,4))
         ENDIF
c	   write(6,*) ' P(it,i) before  jt =',jt,P(IT,3),P(IT,4)
c         write(6,*) ' yi = ',yi, ' xsi = ',xsi
c	   write(6,*) ' P(1(2),3) before =',P(1,3),P(2,3)
* calculate true rapidities
         EN = SQRT(P(IT,1)**2+P(IT,2)**2+P(IT,3)**2)
         If (P(IT,3).ge.EN) then
            rap_lab=1000.
         Elseif (P(IT,3).le.-EN) then
            rap_lab=-1000.
         Else
            rap_lab=0.5*Log((EN-P(IT,3))/(EN+P(IT,3)))
        Endif
         pt_lab = sqrt(P(IT,1)**2 + P(IT,2)**2)
         ptot = sqrt(P(IT,1)**2 + P(IT,2)**2 + P(IT,3)**2)
         et_lab = pt_lab/ptot*P(IT,4)

c end rap


         P(IT,5) =
     +   dsqrt(dabs(QT2B+(P(IT,3)-P(IT,4))*(P(IT,3)+P(IT,4))))
c         write(6,*) ' mtest ,p(it,5) ',mtest,P(IT,5)
         IF(ITIMSHR.EQ.1) THEN
c insert here time like showers
c  add q0 as a safetly margin to avoid xsi old = xsi_new           
            DMAS_MAX = DSQRT(DMAX1(XSI_OLD*YI**2*SSS - QT2B-Q0**2,0.d0))
            if(xsi_old.eq.xsi_hard(jt)) then
c              write(6,*) ' xs_old = xsi_hard'
              DMAS_MAX=0
            Endif
c            Dmas_max=Dmas_max*0.1
*! in radians
c            the_max = PYP(IT,13)
            MSTJO = MSTJ(48)
c impose maximum angle for 1st emission MSTJ(48) = 1
c            MSTJ(48) = 1
c            PARJ(85) = the_max
c            write(6,*) ' before tim shower max mass: ',DMAS_MAX,IT,N
c            call PYLIST(1)
            if(itimtest.eq.1) then
c boost to showering parton rest frame
		   Etest = DMAS_MAX**2 + P(IT,1)**2+P(IT,2)**2+P(IT,3)**2
               Etest = dsqrt(Etest)
               bitx = P(IT,1)/Etest
               bity = P(IT,2)/Etest
               bitz = P(IT,3)/Etest
            else            
               Etest = P(IS(1),4)+P(IS(2),4)
               bitx = (P(IS(1),1)+P(IS(2),1))/Etest
               bity = (P(IS(1),2)+P(IS(2),2))/Etest
               bitz = (P(IS(1),3)+P(IS(2),3))/Etest
            endif
            IF(DMAS_MAX.ge.parj(82)) then 
c               write(6,*) ' beta boost ',bitx,bity,bitz,Dmas_max
               CALL PYROBO(0,N,0.D0,0.D0,-bitx,-bity,-bitz)
c            write(6,*) ' before tim shower max mass: ',DMAS_MAX,IT,N
c            write(6,*) ' cascps MSTJ(45) ',mstj(45)
c            call PYLIST(1)
               CALL PYSHOW(IT,0,DMAS_MAX)
               CALL PYROBO(0,N,0.D0,0.D0,bitx,bity,bitz)
            Endif
            MSTJ(48) = MSTJO
c            do i=it,N
c              if(K(i,2).ge.3.and.k(i,2).le.10) then
c               write(6,*) ' caps: tim like: heavy flavor produced ',
c     &         k(i,2)
c              endif
c            end do 
              
            IF(N.GE.IT+2) then
c               write(6,*) ' max mass: ',DMAS_MAX,' gen mass ',P(IT+2,5),
c     &          ' Nr of timelike partons ',(N-IT+2)/2
c               write(6,*) ' shower cutoff parj(82)',parj(82)
c at it+2 is the virtual parton stored, that is showered            
               P(IT,5)=P(IT+2,5)
c boost showered system back such that pt is the same as of shower initiator             
               pxff = P(IT,1)
               es = P(IT+2,4)
               pxs = P(IT+2,1)
               betax = (-pxs*es + pxff*dsqrt(es**2-pxs**2+pxff**2))/
     &        (es**2+pxff**2)
               call pyrobo(it+2,n,0.d0,0.d0,betax,0d0,0d0)
               pyff = P(IT,2)
               es = P(IT+2,4)
               pys = P(IT+2,2)
               betay = (-pys*es + pyff*dsqrt(es**2-pys**2+pyff**2))/
     &        (es**2+pyff**2)
               call pyrobo(it+2,n,0.d0,0.d0,0d0,betay,0d0)
ccc               P(IT,5)=DSQRT(DABS(DOT1(IT+2,IT+2)))
c calculate xsi using the mass of the shower            
               XSI =(P(IT,5)**2+QT2B) /SSS/YI**2
c calculate 4-momentum of parton using yi from before to boost shower to have same pz
               IF(JT.EQ.1) THEN
                  P(IT,3) =YI*(P(1,3)+XSI*P(2,3))
                  P(IT,4) =YI*(P(1,4)+XSI*P(2,4))
               ELSEIF(JT.EQ.2) THEN
                  P(IT,3) =YI*(P(2,3)+XSI*P(1,3))
                  P(IT,4) =YI*(P(2,4)+XSI*P(1,4))
               ENDIF
c boost showered system back such that yi and xsi  are the same as of shower initiator with proper mass             
               pzff = P(IT,3)
               es = P(IT+2,4)
               pzs = P(IT+2,3)
               betaz = (-pzs*es + pzff*dsqrt(es**2-pzs**2+pzff**2))/
     &        (es**2+pzff**2)
               call pyrobo(it+2,n,0.d0,0.d0,0d0,0d0,betaz)
   
c               call PYLIST(1)
               
            endif
            IF(N.GT.MSTU(4)-10) THEN
               WRITE(6,*) ' cascps: no more memory in PYJETS'
               LST(21)=52
               RETURN
            ENDIF

c               P(IT,5) = DMAS_MAX
c               XSI =(P(IT,5)**2+QT2B) /SSS/YI**2
c               write(6,*) ' check xsi ',xsi_old,QT2B/SSS/YI**2,xsi
               
            DO I=N+1,N+8
               DO J=1,5
                  K(I,J)=0
                  P(I,J)=0.d0
               ENDDO
            ENDDO
c            write(6,*) ' after tim shower 2nd',IT
c        call PYLIST(1)
         ENDIF
c	   write(6,*) ' P(it,3) after tim ',P(IT,3),P(IT,4)
c         write(6,*) ' shower check it  ', (P(it,i),i=1,4)
c         write(6,*) ' shower check sum  ', (P(it+2,i),i=1,4)
c         IF(JT.EQ.1) THEN
c            P(IT,3) =YI*(P(1,3)+XSI*P(2,3))
c            P(IT,4) =YI*(P(1,4)+XSI*P(2,4))
c         ELSEIF(JT.EQ.2) THEN
c            P(IT,3) =YI*(P(2,3)+XSI*P(1,3))
c            P(IT,4) =YI*(P(2,4)+XSI*P(1,4))
c         ENDIF
c         write(6,*) ' masses ', dot1(it,it),dot1(it+2,it+2)
c         write(6,*) ' pt it  ', pyp(it,10),sqrt(qt2b)
c         write(6,*) ' pt check sum  ', pyp(it+2,10)

         Do I=1,4
            PSFIN(I) = PSFIN(I) + P(IT,I)
         EndDo
         if(lporder) then
            write(6,*) ' y ',yi,' xg ',xg,' xsi ',xsi
            write(6,*) ' y_bar ',yi*xsi,' xgbar ',xgbar
            write(6,*) ' new xg ',xg+yi, ' new xgbar ',xgbar+yi*xsi
            write(6,*) ' new xsi = ',(yph+xgbar+yi*xsi)/(xg+yi+yph_bar)
         endif
         IF(LSPRINT) THEN
            write(6,*) '       QT2B = ',QT2B,' Z = ',Z
            write(6,*) '       qt**2 = ',QT2B,P(IT,1)**2+P(IT,2)**2
         endif

         theta = PYP(IT,14)
         K(IT,3)=N+1
c         IF(THETA.LT.THE2(JT).AND.IORDER.EQ.1.AND.JT.EQ.2) Then
cc         IF(xsi.gt.xsi_old.AND.IORDER.EQ.1.AND.JT.EQ.2) Then
cc            write(6,*) ' cascps: angular ordering violated',JT
cc            write(6,*) ' old angle: ',THE2(JT),' new angle: ',THETA
cc            write(6,*) ' xsi old ',xsi_old,' xsi ',xsi
cc            write(6,*) ' mass ',P(IT,5),' max mass ',dmas_max
cc            write(6,*) ' y = ',yi,' yi*xsi = ',yi*xsi,' z = ',z
cc            write(6,*) ' xb = ',xb,' QT2B = ',QT2B,' s = ',sss
cc            write(6,*) ' z*QT" =',z*sqrt(QT2B),' Z =',Z
cc            write(6,*) ' t ',tevb
cc            elseif(jt.eq.2) then
ccc            write(6,*) 'JT = ',jt,' xsi old ',xsi_old,' xsi ',xsi
ccc         ELSEIF(THETA.GT.THE2(JT).AND.IORDER.EQ.1.AND.JT.EQ.1) Then
cc         ELSEIF(xsi.gt.xsi_old.AND.IORDER.EQ.1.AND.JT.EQ.1) Then
cc            write(6,*) ' cascps: angular ordering violated',JT
cc            write(6,*) ' old angle: ',THE2(JT),' new angle: ',THETA
cc            write(6,*) ' xsi old ',xsi_old,' xsi ',xsi
cc            write(6,*) ' mass ',P(IT,5),' max mass ',dmas_max
cc            write(6,*) ' y = ',yi,' yi*xsi = ',yi*xsi,' z = ',z
cc            write(6,*) ' xb = ',xb,' QT2B = ',QT2B,' s = ',sss
cc            write(6,*) ' z*QT" =',z*sqrt(QT2B),' Z =',Z
cc            write(6,*) ' t ',tevb
cc            elseif(jt.eq.1) then
ccc            write(6,*) 'JT = ',jt,' xsi old ',xsi_old,' xsi ',xsi
cc         endif
         THE2(JT)=THETA
         IF(LPRINT) then
            write(6,*) ' p(it,3) ',P(IT,3)
            write(6,*) ' p(it,4) ',P(IT,4)
            write(6,*) 'YI,XSI,QT2B,Z',YI,XSI,QT2B,Z
            write(6,*) 'XA,XB,XS(JT)',XA,XB,XS(JT)
         endif
C...RECONSTRUCT KINEMATICS OF BRANCHING: SPACELIKE PARTON
         K(N+1,1)=14 
c hannes hack
         K(N+1,2)=IFLA
         IF(IFLA.EQ.0) K(N+1,2)=21
c hannes hack
ccc         write(6,*) ' Parton branching: branching space - like ',N+1,K(N+1,2)
         K(N+1,3)=0
c
         P(N+1,1)= P(IT,1)+P(IS(JT),1)
         P(N+1,2)= P(IT,2)+P(IS(JT),2)
         P(N+1,3)= P(IT,3)+P(IS(JT),3)
         P(N+1,4)= P(IT,4)+P(IS(JT),4)
c         P(N+1,1)= -P(IT,1)+P(IS(JT),1)
c         P(N+1,2)= -P(IT,2)+P(IS(JT),2)
c         P(N+1,3)= -P(IT,3)+P(IS(JT),3)
c         P(N+1,4)= -P(IT,4)+P(IS(JT),4)
c         write(6,*) ' check N+1 ',dot1(n+1,n+1)
         KT2B = P(N+1,1)**2+P(N+1,2)**2
         IF(LPRINT) write(6,*) ' checking kt"s old new ',kt2(JT),kt2b
         k2d =-P(N+1,1)**2-P(N+1,2)**2+
     +   (P(N+1,4)-P(N+1,3))*(P(N+1,4)+P(N+1,3))
         P(N+1,5) = SIGN(1.d0,k2d)*dSQRT(dabs(k2d))
c         write(6,*) ' check ',N+1,k2d,dSQRT(dabs(k2d))
         KT2(JT) = KT2B
         k2 = dabs(P(N+1,5))
         kt2bs = SQRT(kt2b)
c define xg as lightcone variable
c         write(6,*) ' jt ',jt,I_NG,KXG(JT,I_NG),KXG_BAR(JT,I_NG)
         I_NG=I_NG + 1
         NPROP(JT,I_NG) = N+1
         KXG(JT,I_NG) = KXG(JT,I_NG-1) + YI
c         write(6,*) ' jt new',jt,I_NG,KXG(JT,I_NG-1)/Z,KXG(JT,I_NG)
         KXG_BAR(JT,I_NG) = KXG_BAR(JT,I_NG-1) + YI*XSI
         QK(JT,I_NG) =P(N+1,5)
         QKT(JT,I_NG)=SQRT(KT2B)
         QQT(JT,I_NG)=sqrt(QT2B)
         QQTP(JT,I_NG)=sqrt(QT2B)
         TH(JT,I_NG)=theta
         ZZ(JT,I_NG)=z
         IF(Lporder) then
            IF(CCFM) THEN
               write(6,10100) sqrt(QT2B),xg,xgbar,sqrt(QT2B),SQRT(KT2B)
     +         ,P(N+1,5),Z,theta
            ELSE
            ENDIF
10100 format(' qt" =',F7.3,' xg = ',F9.7,' xg_bar =',
     +        F9.7,' qt =',F6.3,
     +       ' kt    =',F6.3,F9.3,' z = ',F5.3,' th =',F8.3)
         endif
         DO J=1,5
            K(N+2,J)=0
            P(N+2,J)=0.D0
         ENDDO
         K(N+2,1)=-1
         K(N+2,2)=K(IS(JT)+1,2)
         K(N+2,3)=N+1
ccc         write(6,*) ' before color flow K(N+2) ',K(N+2,2),N+2,IS(JT)+1

C...DEFINE COLOUR FLOW OF BRANCHING
         K(IS(JT),1)=14
         K(IS(JT),3)=N+1
c         write(6,*) ' 14  line nr ',IS(JT)
         ID1=IT
         KN1=ISIGN(500+IABS(K(N+1,2)),2*K(N+1,2)+1)
         KD1=ISIGN(500+IABS(K(ID1,2)),2*K(ID1,2)+1)
         IF(K(N+1,2).EQ.21) KN1=500
         IF(K(ID1,2).EQ.21) KD1=500

         IF((KN1.GE.501.AND.KD1.GE.501).OR.(KN1.LT.0.AND. KD1.EQ.500)
     +   .OR.(KN1.EQ.500.AND.KD1.EQ.500.AND. dcasrn().GT.0.5)
     +   .OR.(KN1.EQ.500.AND.KD1.LT.0)) ID1=IS(JT)
         ID2=IT+IS(JT)-ID1
         P(N+2,3)=ID1
         P(N+2,4)=ID2
         P(ID1+1,1)=N+1
         P(ID1+1,2)=ID2
         P(ID2+1,1)=ID1
         P(ID2+1,2)=N+1
cc         write(6,*) ' color flow branch ',ID1,ID2
         K(N+1,4)=K(N+1,4)-MOD(K(N+1,4),MSTU(5))+ID1
         K(N+1,5)=K(N+1,5)-MOD(K(N+1,5),MSTU(5))+ID2
         K(ID1,4)=MOD(K(ID1,4),MSTU(5))+(N+1)*MSTU(5)
         K(ID1,5)=MOD(K(ID1,5),MSTU(5))+ID2*MSTU(5)
         K(ID2,4)=MOD(K(ID2,4),MSTU(5))+ID1*MSTU(5)
         K(ID2,5)=MOD(K(ID2,5),MSTU(5))+(N+1)*MSTU(5)
         N=N+2

         MSTU(1)=0
         IFLS(JT+2)=IFLS(JT)
         IFLS(JT)=IFLA
c         write(6,*) ' end IFLS(JT) ',ifls(jt),ifls(jt+2)
         XS(JT)=XA
         XS_BAR(JT) = KXG_BAR(JT,I_NG) 
c         write(6,*) ' test xs ',xs(JT),KXG(jt,I_NG)
         KT2(JT)=KT2B
         TEVS(JT)=TEVB
c       write(6,*) ' end JT ',JT,' N = ',N
      ENDIF
c      write(6,*) ' after iflag if ' 
      if((XS(1)+XS_BAR(2)).gt.1.or.(XS(2)+XS_BAR(1)).gt.1) then 
         IF(LPRINT) then
            write(6,*) ' cascps: xs = ',XS(1)+XS_BAR(2),' xs = ',XS(2)+XS_BAR(1),JT
            write(6,*) ' cascps: xs(1) = ',XS(1),XS_BAR(2),' xs(2) = ',XS(2),XS_BAR(1),JT
            write(6,*) ' cascps: yi, xsi = ',YI,XSI,Z, XB,QT2B
         endif 
       endif

C...SAVE QUANTITIES, LOOP BACK
   40 CONTINUE
      IF(LPRINT) THEN
         write(6,*) ' save quantities and loop back ' 
         call PYLIST(1)
      ENDIF
      IS(JT)=N-1
      QT2P(JT)=QT2B
      IF(IFLAG.EQ.1) THEN
      ELSE
         IF(JT.EQ.1) IPU1=N-1
         IF(JT.EQ.2) IPU2=N-1
      ENDIF
      IF(LPRINT) THEN
         write(6,*) 'check end: N=',N
         write(6,*) 'check end: xa,xb,z ',xa,xb,z
      ENDIF
      IF(LSPRINT) then
         write(6,*) ' check end',QT2P(1),QT2P(2),
     +   MAX(IABS(1-ILEP)*QT2P(1),MIN(1,2-ILEP)*QT2P(2)),
     +   N,NS+2
      ENDIF
c         write(6,*) ' save quantities and loop back ' 
c         call PYLIST(2)
      IF(IFLAG.EQ.1)
     +   GOTO 30
c      write(6,*) ' end iflag = 0 '
c      write(6,*) ' before reshuff '
c	call pylist(1) 
      xsi_test = 1.d0
      IF(reshuff) then
C
         XGBAR = 0.d0
         DO J=I_NG,1,-1
            IF(J.GT.0) THEN
               I=NPROP(JT,J)
               KXG_BAR(JT,J) = XGBAR
               IF(JT.EQ.1) THEN
                  P(I,4) = KXG_BAR(JT,J)*P(2,4) + KXG(JT,J)*P(1,4)
                  P(I,3) = KXG_BAR(JT,J)*P(2,3) + KXG(JT,J)*P(1,3)
               ELSE
                  P(I,4) = KXG_BAR(JT,J)*P(1,4) + KXG(JT,J)*P(2,4)
                  P(I,3) = KXG_BAR(JT,J)*P(1,3) + KXG(JT,J)*P(2,3)
               ENDIF
               XGBAR=KXG_BAR(JT,J) + XGBAR
               q2_test =-P(I,1)**2-P(I,2)**2+ (P(I,4)-P(I,3))*(P(I,4)+
     +         P(I, 3))
               P(I,5) = SIGN(1.d0,q2_test)*dSQRT(dabs(q2_test))
               QK(JT,J) = P(I,5)
c		   if(jt.eq.1) then
c                  xsi_test=(kxg(JT,j)+xgbar)/(kxg_bar(JT,j)+xg)
c		   elseif(jt.eq.2) then
                  xsi_test=(kxg_bar(JT,j)+yph)/(kxg(JT,j)+yph_bar)
c		   endif
c		   write(6,*) ' xsi_test ', xsi_test,jt
            ENDIF
         ENDDO


      endif

c      write(6,*) ' after reshuff ',JT
c	call pylist(2) 
      IF(lporder) then
c         write(6,*) ' after reshuffling propagator '
         DO I=2,I_NG
            write(6,10100) QQTP(JT,I),kxg(JT,I),kxg_bar(JT,I),QQT(JT,I),
     +      QKT(JT,I),QK(JT,I),ZZ(JT,I),th(JT,I)
         ENDDO
      ENDIF
      k2 = dabs(QK(JT,I_NG))
      kt2bs = QKT(JT,I_NG)

      IF(JT.NE.JTFIN) GOTO 10

      IF(reshuff.and.ilep.eq.0.and.iresl1.eq.0.and.iresl2.eq.0) then
         PZF=P(NPROP(1,1),3)+P(NPROP(2,1),3)
c         write(6,*) ' NPROP(1,1),NPROP(2,1) ',NPROP(1,1),NPROP(2,1)
         DBETAZ=(-PEM*PZM+
     +   PZF*DSQRT(PZF**2+PEM**2-PZM**2))/(PZF**2+PEM**2)
c         write(6,*) ' before boost ',npsfin
c         call pylist(1)
         CALL PYROBO(25,npsfin,0.D0,0.D0,0.D0,0.D0,DBETAZ)
ccc         CALL PYROBO(21,npsfin,0.D0,0.D0,0.D0,0.D0,DBETAZ)
         PZF=P(NPROP(1,1),4)+P(NPROP(2,1),4)
c         write(6,*) ' PZF check final',PZF,P(21,4)+P(23,4),P(25,4)+P(26,4) 
         
c         call pylist(1)
c	   write(6,*) ' here in reshuffling... 1'
      elseif(reshuff.and.ilep.eq.0.and.
     +	(iresl1.eq.1.or.iresl2.eq.1)) then
c	   write(6,*) ' here in reshuffling... 2 '
      elseif(reshuff.and.ilep.ne.0) then
         Do i=1,4
            P(NPROP(2,I_NG),I) = PSFIN(I) + P(25,I) + P(26,I) - P(21,I)
         Enddo
         P(NPROP(2,I_NG),5) = -DSQRT(Dabs(DOT1(NPROP(2,I_NG),NPROP(2,I_
     +   NG))))
c	   write(6,*) ' here in reshuffling... 3 '
      endif
      do i=21,N
       Do j=1,4
         if(k(i,1).eq.3) then 
           PSUMF(J) = P(i,j) + psumF(j)
         endif
       end do
      end do
      cmtest = dot(psumf,psumf)

C...STORE USER INFORMATION
      K(21,1)=14
      IF(IRES(1).EQ.0) K(21,1)=11
      K(23,1)=14
      K(21,3)=NS+1
      K(23,3)=NS+3
      DO JT=1,2
         KFL(1,JT)=IFLS(JT)
         IF(IFLS(JT).EQ.0) KFL(1,JT)=21         
      ENDDO
      PYVAR(31)=min(0.999d0,XS(1)+XS_BAR(2))
      PYVAR(32)=min(0.999d0,XS(2)+XS_BAR(1))
c      write(6,*) ' cascps 2nd ',PYVAR(31),PYVAR(32),XS(1),XS_BAR(2),XS(2),XS_BAR(1)
      DO I=NS+1,N
         DO J=1,5
            V(I,J)=0.d0
         ENDDO
      ENDDO
   50 CONTINUE
      if(lporder) then
         write(6,*) ' qt_sum_x = ',pxsum,' qt_sum_y = ',pysum
         write(6,*) ' kt_sum_x = ',xktsum,' kt_sum_y = ',yktsum
         write(6,*) ' kt_x = ',P(NP,1),' kt_y = ',P(NP,2)
c      call PYLIST(1)
      endif
c      write(6,*) ' end of cascps '
c      call pylist(1)
c now go back from  lab system
      CALL PYROBO(0,0,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
c      write(6,*) ' after boost 1'

      DO I = 1,4
         P(1,I) = PB_SAV(1,I)
         P(2,I) = PB_SAV(2,I)
      ENDDO

      dboost=DBETA(2,3)
      if(ILEP.eq.0.and.Iresl1.eq.0.and.Iresl2.eq.0) then
         dboost = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
      Endif
      if(dabs(dboost).gt.0.999999) then
         LST(21)=53
c         return
      Endif
      CALL PYROBO(0,0,0.D0,0.D0,-DBETA(2,1),-DBETA(2,2),-dboost)
      CALL PYROBO(0,0,0.D0,-SPHI(2),0.D0,0.D0,0.D0)
      CALL PYROBO(0,0,-STHETA(2),0.D0,0.D0,0.D0,0.D0)
c      write(6,*) ' after boost 2',dboost,DBETA(2,3),PYVAR(31),PYVAR(32)
      if(ILEP.eq.0.and.(Iresl1.eq.1.or.Iresl2.eq.1)) then
c	write(6,*) ' cascps: bef pyvar ',pyvar(31),pyvar(32),yx1,yx2    
c	write(6,*) ' cascps: ipy(41),(42)',IPY(41),IPY(42)     
	   if(Iresl1.eq.1) pyvar(31) = pyvar(31)/yx1
	   if(Iresl2.eq.1) pyvar(32) = pyvar(32)/yx2
         dboostl = (PYVAR(31)-PYVAR(32))/(PYVAR(31)+PYVAR(32))
         CALL PYROBO(5,0,0.d0,0.D0,0.D0,0.D0,-dboostl)
c      write(6,*) ' cascps: end after dboostl boost ',dboostl
c  	write(6,*) ' cascps 1 end: pyvar ',pyvar(31),pyvar(32),yx1,yx2    
c	call pylist(1)
      Endif
   60 CONTINUE
c	write(6,*) ' cascps 2 end: pyvar ',pyvar(31),pyvar(32),yx1,yx2    
      if(xsi_test.le.0)  THEN
c         LST(21)=58
c	   write(6,*) ' xsi_test < 0 ',xsi_test
      Endif
      ISEMIH = ISEMIHO
      RETURN
      END
