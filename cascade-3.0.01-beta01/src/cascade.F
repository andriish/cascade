      SUBROUTINE CASCADE
      IMPLICIT None
#include "caslujets.inc"
#include "casluco.inc"
#include "caspart.inc"
#include "caspara.inc"
#include "caspara1.inc"
#include "caseffic.inc"
#include "casgki.inc"
#include "casludat1.inc"
#include "casludat3.inc"
#include "casshwr.inc"
#include "casprkt.inc"
#include "casvalues.inc"
#include "casscal.inc"
#include "caspara5.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspsipol.inc"
#include "caspar.inc"
c common for bases/spring 5.1
      Integer Mxdim
      PARAMETER (MXDIM = 50 )
      Double precision xlb,xub
      Integer ndimb,nwildb,ig,ncallb
      COMMON /BPARM1/ XLB(MXDIM),XUB(MXDIM),NDIMB,NWILDB,
     +   IG(MXDIM),NCALLB
      Double precision xlmc,xumc
      Integer ndimmc,ncallmc
      Common /MCINT1/XLMC(MXDIM),XUMC(MXDIM),NDIMMC,NCALLMC
      Double Precision acc1b,acc2b
      Integer itmx1b,itmx2b
      COMMON /BPARM2/ ACC1B,ACC2B,ITMX1B,ITMX2B

	Integer kf
	Common /cbran1/ kf
      Double Precision ACC1,ACC2
      Integer iint,ncb
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer ndim
      COMMON/INTDIM/ NDIM
      Double Precision ANS,ERROR
      Integer Ndimc
      COMMON /DIMEN/ NDIMC
      double Precision XD(20)
      Double Precision XGR(20)
      CHARACTER CHIN*100
      CHARACTER * 16 CNAM
      CHARACTER * 7 CINT
      Character*12 CPSI,CCHI
      CHARACTER*20 VERSQQ
      Character*6 DCPSI 
      Integer Lerr
      COMMON/ERR/LERR(1000)
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Integer IGLU
	Common/CAGLUON/Iglu
      double precision DFUN,FXN1,FXNB
      EXTERNAL DFUN,FXN1,FXNB
      Integer I,J,KPCH,KPAO,IT1,IT2
      Double Precision SQRTS,DDUM,CTIME
      Double Precision DOT
	Integer MSTP,MSTI
	Double Precision PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      Integer kc
      Double precision asb, asbtest
	Double Precision mz,QCDlam
      Common/casalphas/mz,qcdlam
      external asb
#include "casvers.inc"

      IGENFL = 0
      DO 10  I=1,1000
   10 LERR(I) = 0
      ERROR = 0.0
      CHIN='PARU(11)=0.010;'
      CALL PYGIVE(CHIN)
      ALPH=7.299D-3
      ALPHS=0.3D0
      PI=4.D0*DATAN(1.D0)
      IWEI = 0
C...  GIVE ELECTRON FOUR VECTOR
      K(1,1)=21
      K(1,2)=KE
      P(1,1) = 0.0D0
      P(1,2) = 0.0D0
      P(1,3) = DBLE(PLEPIN)
      P(1,5) = PYMASS(KE)
c      P(1,5) = 0
      P(1,4) = DSQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
C...  GIVE PROTON FOUR VECTOR
      K(2,1) = 21
      K(2,2) = KP
      P(2,1) = 0.0D0
      P(2,2) = 0.0D0
      P(2,3) = DBLE(PPIN)
      P(2,5) = PYMASS(KP)
c      P(2,5) = 0
      P(2,4) = DSQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
      N=2
c      CALL PYLIST(1)
C... CALCULATE CMS ENERGY
      DO 20 I=1,4
         CM(I) =P(1,I)+P(2,I)
   20 CONTINUE

C BOOST TO EP CMS
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4))
      DO 30 I=1,2
         DO 30 J=1,5
            PBEAM(I,J)=P(I,J)
   30 KBEAM(I,J)=K(I,J)
      KINT(1,2) = 22
      KINT(2,2)= KP
      KEB = KPH
      SSS = DOT(CM,CM)
      SQRTS=DSQRT(SSS)
      IF(IPRO.EQ.10.OR.IPRO.EQ.11.OR.IPRO.EQ.10.OR.IPRO.EQ.15.OR.
     & IPRO.EQ.2.or.IPRO.eq.3.OR.ipro.eq.4.or.
     & IPRO.EQ.21.or.
     & IPRO.eq.102.or.
     & IPRO.ge.500.and.IPRO.le.515.or.
     & IPRO.le.-1) THEN
      ELSE
         WRITE(6,*) ' wrong subprocess selected ', IPRO
         WRITE(6,*) ' PROGRAM STOPPED '
         STOP
      ENDIF
c deal withinut/output LHE files      
      CALL PYNAME(KINT(2,2),CNAM)
      CINT = CNAM(1:7)
      KPCH = 213
      CALL PYNAME(KPCH,CNAM)
      WRITE(6,10000)

      WRITE(6,10100)

      WRITE(6,10200)

      WRITE(6,10300) VERSQQ
      
      IF(IPRO.lt.0) then
         if(IPRO.eq.-1) then
         write(6,*) '*     Input of LHE file                           '
     +   //'  *'
         ILHA=-2
         call cas_lheinit

         else
         write(6,*) '*     Output of LHE file                          '
     +   //'  *'
         IPRO=Iabs(IPro)
         ILHA=-1
         endif
	Endif
	kf = 2212
	If(Ipro.ne.-1) then
      if(Iabs(K(1,2)).eq.2212) then
c proton proton
         write(6,*) '*       proton proton version selected            '
     +   //'  *'
         KEB=21
	   IRES(1)=1
	   IRES(2)=1
         QMI = 0.
         YMI = 0.
         YMA = 1. 
         
         IF(IPRO.EQ.2) THEN
            if(ihfla.eq.3) then
               kpsi = 333
            elseif(ihfla.eq.4) then
               kpsi = 443
            elseif(ihfla.eq.5) then
               kpsi = 553
            else
               kpsi = 443
            endif

            if(i23s.eq.0) then
            elseif(i23s.eq.2) then
               kpsi = 100000 + kpsi
               write(6,*) '*          2S state selected                        *'                
            elseif(i23s.eq.3) then
               kpsi = 200000 + kpsi
               write(6,*) '*          3S state selected                        *'                
            else
               write(6,*) ' cascade: wrong values for i23s selected: i23s = ',i23s
               write(6,*) ' cascade: possible choices are i23s =0,2,3'
               write(6,*) ' cascade: program stop '
               stop
            endif   
            CALL PYNAME(Kpsi,CNAM)
            Cpsi = CNAM(1:12)
            IF(INTER.EQ.0) WRITE(6,11001) Cpsi
            if(Ipsipol.eq.1) then
               write(6,*) '*      including lepton polarisation ',
     &         Ipsipol,'  *' 
               write(6,*) '*      only leptonic decay calculated  ',
     &         '             *' 
     	         write(6,*) '*    No final state PS ',
     &         ' with lepton polarisation    *'
     	         write(6,*) '*      Final state PS ',
     &         'from gluon switched off       *'
               if(IFPS.eq.2) IFPS=0
               if(IFPS.eq.3) IFPS=1
c               write(6,*) ' original decay table '
c               call pylist(12)
               KC = pycomp(kpsi)
c           switch off all decay channels               
               do i=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1 
                  mdme(i,1) = 0
               end do
               mdme(mdcy(kc,2),1) = 1
               mdme(mdcy(kc,2)+1,1) = 1
c               write(6,*) ' modified decay table ',mdcy(kc,2),mdcy(kc,3)
c               call pylist(12)
            endif
         ELSEIF(IPRO.EQ.3) THEN
            if(ihfla.eq.3) then
               kchi = 10331
            elseif(ihfla.eq.4) then
               kchi = 10441
c               kchi = 445
            elseif(ihfla.eq.5) then
               kchi = 10551
            else
               kchi = 10441
            endif
            
            if(i23s.eq.0) then
            elseif(i23s.eq.2.and.ihfla.eq.5) then
               if(ihfla.eq.5) kchi = 100000+kchi
c               if(ihfla.eq.4) kchi = 100000+kchi
               write(6,*) '*          2P state selected                        *',kchi              
            elseif(i23s.eq.3.and.ihfla.eq.5) then
               if(ihfla.eq.5) kchi = 200000+kchi
c               if(ihfla.eq.4) kchi = 100000+kchi
               write(6,*) '*          3P state selected                        *',kchi              
c            elseif(i23s.eq.3) then
c               kpsi = 300000 + kchi
c               write(6,*) '*          3P state selected                        *'                
            else
               write(6,*) ' cascade: wrong values for i23s selected: i23s = ',i23s
               if(ihfla.eq.4) write(6,*) ' cascade: possible choices are i23s =0 for ihfla = 4'
               if(ihfla.eq.5) write(6,*) ' cascade: possible choices are i23s =0,2 for ihfla = 5'
               write(6,*) ' cascade: program stop '
               stop
            endif   

            CALL PYNAME(Kchi,CNAM)
            Cchi = CNAM(1:12)
         
            IF(INTER.EQ.0) WRITE(6,11002) Cchi
         ENDIF
         IF(IPRO.EQ.4) then
           write(6,11515)
         endif
         IF(IPRO.EQ.10) Then
            IHFLA = 3
c            IF(IRPA.ne.0) WRITE(6,11300)
            IF(IRPA.NE.0) THEN
               IF(IHFLA.LT.4) write(6,11205)
               IF(IHFLA.GE.4) write(6,11305)
            ENDIF
            IF(IRPB.NE.0) write(6,11405)
            IF(IRPC.NE.0) THEN
               IF(IHFLA.LT.4) write(6,11505)
               IF(IHFLA.GE.4) write(6,11605)
            ENDIF
            IF(IRPD.NE.0) write(6,11705)
            IF(IRPE.NE.0) THEN
               IF(IHFLA.LT.4) write(6,11805)
               IF(IHFLA.GE.4) write(6,11905)
            ENDIF
            IF(IRPF.NE.0) THEN
               IF(IHFLA.LT.4) write(6,12005)
               IF(IHFLA.GE.4) write(6,12105)
            ENDIF
         Endif
         
         IF(IPRO.EQ.11) THEN
            IRPA = 1
            WRITE(6,11700)
            WRITE(6,11800) MAX(4,IHFLA)
            KPA=4
         ENDIF
c         IF(IPRO.EQ.15) WRITE(6,11304)
         IF(IPRO.EQ.21)  then
            if(ihfla.eq.3) then
               kpsi = 333
            elseif(ihfla.eq.4) then
               kpsi = 443
            elseif(ihfla.eq.5) then
               kpsi = 553
            else
               kpsi = 443
            endif
            CALL PYNAME(Kpsi,CNAM)
            Dcpsi = CNAM(1:6)
          WRITE(6,*) '*            gluon + gluon  --> ',Dcpsi,' + ',Dcpsi,'     *'
         ENDIF


         IF(IPRO.EQ.102) WRITE(6,11302)
         IF(IPRO.EQ.500) WRITE(6,14500)
         IF(IPRO.EQ.501) WRITE(6,14501)
         IF(IPRO.EQ.503) WRITE(6,14503)
         IF(IPRO.EQ.504) WRITE(6,14504)
         IF(IPRO.EQ.510) WRITE(6,14510)
         IF(IPRO.EQ.511) WRITE(6,14511)
         IF(IPRO.EQ.513) WRITE(6,14513)
         IF(IPRO.EQ.514) WRITE(6,14514)
         IF(IPRO.EQ.504.OR.IPRO.EQ.514) THEN
            WRITE(6,11811) MAX(3,IHFLA)
            KPA=MAX(4,IHFLA)
         ENDIF
c	   write(6,11801) Mult
      elseif(Iabs(K(2,2)).eq.22.or.
     & (Iabs(K(1,2)).eq.11.and.Iabs(K(2,2)).eq.11)) then
	if(Iabs(K(2,2)).eq.22) then 
c gamma gamma 
         write(6,*) '*     gamma gamma version selected                '
     +   //'  *'
         else
         write(6,*) '*     e+ e-   version selected                    '
     +   //'  *'
         endif
         KEB=22
	   kf=keb
         KINT(1,2) = 22
         KINT(2,2) = 22 
         IRPA = 1 
         IF(IPRO.EQ.10) WRITE(6,11301)
         IF(IPRO.EQ.11) THEN
            WRITE(6,11701)
            WRITE(6,11800) MAX(4,IHFLA)
            KPA=4
         ENDIF
         IF(IPRO.EQ.102) WRITE(6,11302)
        write(6,*) '*     beam 1:  direct/resolved (0/1): ',Ires(1),
     +   '            *'
        write(6,*) '*     beam 2:  direct/resolved (0/1): ',Ires(2),
     +   '            *'
	   
      else
	   IRES(2)=1
         IRPA = 1
         IF(INTER.EQ.0) WRITE(6,10400)
         IF(INTER.EQ.2) WRITE(6,10500)
         IF(IPRO.EQ.10) THEN
            IF(NFLAV.LE.3) THEN
               IF(INTER.EQ.0) WRITE(6,10600) CINT
               IF(INTER.EQ.2) WRITE(6,10700) CINT
            ELSE
               IF(INTER.EQ.0) WRITE(6,10800) CINT,NFLAV
               IF(INTER.EQ.2) WRITE(6,10900) CINT,NFLAV
            ENDIF
         ENDIF
         if(ihfla.eq.3) then
         kpsi = 333
         elseif(ihfla.eq.4) then
         kpsi = 443
         elseif(ihfla.eq.5) then
         kpsi = 553
         else
         kpsi = 443
         endif

         CALL PYNAME(Kpsi,CNAM)
         Cpsi = CNAM(1:7)
         
         IF(IPRO.EQ.2) THEN
            IF(INTER.EQ.0) WRITE(6,11000) CINT,Cpsi
            IF(INTER.EQ.2) WRITE(6,11100) CINT
            if(Ipsipol.eq.1) then
               write(6,*) '*      including lepton polarisation ',
     &         Ipsipol,'  *' 
               write(6,*) '*      only leptonic decay calculated  ',
     &         '             *' 
     	         write(6,*) '*    No final state PS ',
     &         ' with lepton polarisation    *'
     	         write(6,*) '*      Final state PS ',
     &         'from gluon switched off       *'
               if(IFPS.eq.2) IFPS=0
               if(IFPS.eq.3) IFPS=1
c               write(6,*) ' original decay table '
c               call pylist(12)
               KC = pycomp(kpsi)
c           switch off all decay channels               
               do i=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1 
                  mdme(i,1) = 0
               end do
               mdme(mdcy(kc,2),1) = 1
               mdme(mdcy(kc,2)+1,1) = 1
c               write(6,*) ' modified decay table ',mdcy(kc,2),mdcy(kc,3)
c               call pylist(12)
            endif
         ENDIF
        IF(IPRO.EQ.4) then
           write(6,11516)
        endif
        IF(IPRO.EQ.10) WRITE(6,11200)

         IF(ISEMIH.EQ.1) THEN
            IF(IPRO.EQ.2) then
               write(6,11500)
            else
               write(6,11400)
            endif
         ENDIF
         IF(IPRO.EQ.11) THEN
            IF(INTER.EQ.0) THEN
               WRITE(6,11600) CINT
               WRITE(6,11800) MAX(4,IHFLA)
            ENDIF
            IF(INTER.EQ.2) WRITE(6,11900) CINT
            KPA=4
         ENDIF
         IF(IPRO.EQ.11) WRITE(6,12000)

         IHF = 0
         IF(THEMA.NE.180.0D0) THEN
            WRITE(6,12100) THEMA

         ELSE
            WRITE(6,*) '*  no cut on max angle of scattered electron   '
     +      //'     *'
         ENDIF
         IF(THEMI.NE.0.0D0) THEN
            WRITE(6,12200) THEMI

         ELSE
            WRITE(6,*) '*  no cut on min angle of scattered electron   '
     +      //'     *'
         ENDIF
         IF(QMI.GT.0.0) THEN
            WRITE(6,12300) QMI

         ELSE
            WRITE(6,*) '*  Q^2 _min according to kinematics            '
     +      //'     *'
         ENDIF
         IF(QMA.LT.10D8) THEN
            WRITE(6,12400) QMA

         ELSE
            WRITE(6,*) '*  Q^2 _max according to kinematics            '
     +      //'     *'
         ENDIF
         IF(YMI.GT.0.0) THEN
            WRITE(6,12500) YMI

         ELSE
            WRITE(6,*) '*  y_min according to kinematics               '
     +      //'     *'
         ENDIF
         IF(YMA.LT.1.0) THEN
            WRITE(6,12600) YMA

         ELSE
            WRITE(6,*) '*  y_max according to kinematics               '
     +      //'     *'
         ENDIF
      endif
      IF(IPRO.eq.3.or.IPRO.eq.102) then
      else
      IF(PT2CUT(IPRO).GT.0D0) write(6,17100) DSQRT(PT2CUT(IPRO)),IPRO
      endif 
      IF(Ipro.eq.10.or.IPRO.eq.11) then
         IF(PT2CUT(IPRO).GT.0D0) then
            PT2CUT(10) = PT2CUT(IPRO)
            PT2CUT(11) = PT2CUT(10)
c            write(6,11888)
         Endif
      endif
      endif
      
      IPY(8) = NFLAV

      WRITE(6,*) '*****************************************************'
      IF(IFPS.EQ.0) THEN
         IPY(13) = 0
         IPY(14) = 0
         MSTJ(41) = 0 ! final state PS from PYTHIA for resonance
      ENDIF
      IF(IFPS.NE.0) THEN
         WRITE(6,*) '##################################################'
     +   //'###'
         WRITE(6,*) '#     parton shower selection:                    '
     +   //'  #'
         IF(IFPS.EQ.0) THEN
            IPY(13) = 0
            IPY(14) = 0
         ELSEIF(IFPS.EQ.1) THEN
            IPY(13) = 0 
            IPY(14) = 1		
            MSTJ(41) = 0 ! final state PS from PYTHIA for resonance
            if(iglu.eq.1.or.iglu.ge.6) then
		   else
               IPY(14) = 0
               WRITE(6,*) '############################################'
     +         //'#########'
               WRITE(6,*) '#     no initial state PS possible for IGLU>'
     +         //'1       #'
               WRITE(6,*) '#        initial state PS switched off      '
     +         //'        #'
               WRITE(6,*) '#     see manual.                           '
     +         //'        #'
            endif
         ELSEIF(IFPS.EQ.2) THEN
            IPY(13) = 1
            IPY(14) = 0
         ELSEIF(IFPS.EQ.3) THEN
            IPY(13) = 1
            IPY(14) = 1
            if(ires(1).eq.1.and.iabs(K(1,2)).eq.11) then
               IPY(14) = 0
               WRITE(6,*) '############################################'
     +         //'#########'
               WRITE(6,*) '#     no initial state PS possible for resol'
     +         //'ved     #'
               WRITE(6,*) '#     initial state PS switched off         '
     +         //'        #'
               WRITE(6,*) '#     see manual.                           '
     +         //'        #'
            endif
            if(iglu.eq.1.or.iglu.ge.6) then
		   else
               IPY(14) = 0
               WRITE(6,*) '############################################'
     +         //'#########'
               WRITE(6,*) '#     no initial state PS possible for IGLU>'
     +         //'1 or 6    #'
               WRITE(6,*) '#        initial state PS switched off      '
     +         //'        #'
               WRITE(6,*) '#     see manual.                           '
     +         //'        #'
            endif
         ENDIF
         if(IPY(14).EQ.1.and.Iccfm.eq.1) Then
            WRITE(6,*) '#     CCFM inital state parton showers '
     +      //'including    #'
            WRITE(6,*) '#     angular ordering                         '
     +      //'     #'
            elseif(IPY(14).EQ.1.and.Iccfm.ge.2) Then
            WRITE(6,*) '#     New inital state parton showers '
c     +      //'including    #'
c            WRITE(6,*) '#     angular ordering                         '
     +      //'              #'
            write(6,4678) IordasPS
4678        format(' #     Order Splitting fct in shower = ',I1,'             # ') 
            asbtest =   asb(mz)   
            write(6,4679) asbtest*2.*pi,qcdlam
4679        format(' #     alphas(Mz) = ',F6.4,' with QCDlam = ',F6.4,'      # ')      
         endif
         if(IPY(14).EQ.1.and.Iccfm.eq.0) Then
            WRITE(6,*) '#     DGLAP inital state parton showers '
     +      //'including   #'
            WRITE(6,*) '#     angular ordering                         '
     +      //'     #'
         endif
         if(itimshr.eq.1)
     +   WRITE(6,*) '#     timelike inital state partons may shower    '
     +   //'  #'
         if(itimshr.eq.0)
     +   WRITE(6,*) '#     no timelike inital state partons            '
     +   //'  #'
         IF(IPY(13).EQ.1) then
         WRITE(6,*) '#     final state parton shower  '
     +   //'                   #'
         if(Ifinal.eq.1) WRITE(6,*) '#     FPS scale: 2(m^2_1t+m^2_2t)'
     +   //'                   #'
         if(Ifinal.eq.2) WRITE(6,*) '#     FPS scale: shat            '
     +   //'                   #'
         if(Ifinal.eq.3) WRITE(6,*) '#     FPS scale: 2(m^2_1+m^2_2)  '
     +   //'                   #'
         if(Ifinal.eq.4) WRITE(6,*) '#     FPS scale: SCALUP          '
     +   //'                   #'
         WRITE(6,20001) SCALFAF
20001 format(' #     FPS scale factor   ',F5.2,'                      #')       
         endif
         WRITE(6,*) '##################################################'
     +   //'###'
         WRITE(6,*) ' '
      ENDIF
	if(ILHA.le.0) then
	write(6,*)  '#####################################################'
	write(6,*)  '#     proton remnant parameters                     #'
	write(6,*)  '#     treated inside CASCADE                        #'
	write(6,12701) IRspl
12701	format(' #     energy sharing IREM = ',I2,'                      #')
	write(6,*)  '#####################################################'
	elseif(ILHA.eq.1) then
	write(6,*)  '##################################################'
	write(6,*)  '#             output of events in                #' 
	write(6,*)  '#          Les Houges Accord Format              #'
	write(6,*)  '##################################################'
	elseif(ILHA.gt.1) then
         MSTP(61)=0 ! initial state PS
         MSTP(81)=0 ! multiple interactions
         MSTP(71)=0 
         IPY(13) = 0		
         IPY(14) = 0		
         if(ifps.eq.1) then 
c            MSTP(71)=0 
            ifps = 0
            IPY(14) = 1		
            endif
         if(ifps.eq.2) then 
c            MSTP(71)=1 
            ifps = 0
            IPY(14) = 0		
            IPY(13) = 1		
            endif
         if(ifps.eq.3) then
c            MSTP(71)=1 
            IFPS = 1
            IPY(14) = 1
            IPY(13) = 1		
         Endif
	write(6,*)  '##################################################'
	write(6,*)  '#              proton remnant and                #' 
	write(6,*)  '#                final state PS                  #'
	write(6,*)  '#               treated in PYTHIA                #'
	write(6,*)  '##################################################'
	
	endif
      
      IF(IPRO.ne.-1) then
      
      WRITE(6,12700) SQRTS
c      IF(IRUNAEM.EQ.0) THEN
c         WRITE(6,*) ' alpha_em fixed; alpha_em = ',ALPH
c      ELSE
c         WRITE(6,*) ' running alpha_em selected '
c      ENDIF
      IF(IRUNA.EQ.0) THEN
         MSTU(111) = 0
         PARU(111) = ALPHS
         WRITE(6,*) ' alpha_s fixed; alpha_s = ',ALPHS
      ENDIF
      IF(IQ2.EQ.1) WRITE(6,*) ' scale for alpha_s: 4*m_q**2 '
      IF(IQ2.EQ.2) WRITE(6,*) ' scale for alpha_s: shat '
      IF(IQ2.EQ.3) WRITE(6,*) ' scale for alpha_s: m_q**2 + p_t **'
     +//'2 '
      IF(IQ2.EQ.4) WRITE(6,*) ' scale for alpha_s: Q2 '
      IF(IQ2.EQ.5) WRITE(6,*) ' scale for alpha_s: Q2 + p_t **2'
      IF(IQ2.EQ.6) WRITE(6,*) ' scale for alpha_s: kt**2 '
c      IF(IQ2.EQ.7) WRITE(6,*) ' scale for alpha_s:'
c     +//'two*(sL*rt*ru)/((sL*rt*ru)**2+rt**2+ru**2)'
	IF(IQ2.EQ.7) WRITE(6,*) ' scale for alpha_s: max(stmax1,stmax2) '

      WRITE(6,*) ' scale for alpha_s is multiplied by: ',SCALFA
c      WRITE(6,13300) PARU(112),MSTU(112)
c      PYPAR(21) = PARU(112)



      KPAO = KPA


      IF(IPRO.EQ.10.OR.IPRO.EQ.11) NDIM=8
      NDIMC = 9999
   40 CALL dcasrnv(XGR,20)
      DO 50 I=1,20
   50 XD(I)=XGR(I)
      DDUM=DFUN(20,XD)
*         WRITE(6,*) ' NDIM dynamically calculated ',NDIM,NDIMC
      IF(NDIMC.GT.20) goto 40
      NDIM = NDIMC
      NIN=0
      NOUT=0
      IF(iint.eq.0) then
********************************************************************
*     Initialization of BASES/SPRING V5.1
********************************************************************
*===> Initialization of BASES by calling BSINIT
c
         CALL BSINIT
         NDIMB=NDIM
         NWILDB = NDIM
         NCALLB = NCB
         ITMX1B = 30
         ACC1B = DBLE(ACC1)
         ITMX2B = 100
         ACC2B = DBLE(ACC2)
         DO 60 I=1,NDIMB
           IG(I) = 1
           XLB(I)=1.0D-12
   60    XUB(I)=1.D0
********************************************************************
*     Nimerical Integration by BASES V5.1
********************************************************************

         CALL BASES( FXNB, AVGI, SD, CTIME, IT1, IT2 )
         ANS = AVGI
         ERROR = SD
         SD=ERROR
         AVGI=ANS
      elseif(iint.eq.1) then
         NDIMMC=NDIM
         DO I=1,NDIMMC
           IG(I) = 1
           XLMC(I)=1.0D-12
           XUMC(I)=1.D0
         End do
         NCALLMC = NCB
         Call MCINT(FXNB,AVGI,SD)
         ANS = AVGI
         ERROR = SD
         SD=ERROR
         AVGI=ANS
      else
      write(6,*)' CASCADE: IINT = ',IINT,' not implemented: STOP '
      endif
      
      endif

      KPA = KPAO
      RETURN

10000 FORMAT(' *****************************************************')
10100 FORMAT(' *                                                   *')
10200 FORMAT(' *     You are using the CASCADE MC generator        *')
10300 FORMAT(' *           version  ',A20,'           *')
10400 FORMAT(' *    neutral current interaction selected           *')
10500 FORMAT(' *    charged current interaction selected           *')
10600 FORMAT(' *    gamma + gluon_',A6,' --> q q_bar(light quarks) *')
10700 FORMAT(' *      W + gluon_',A6,' --> q q_bar(light quarks)   *')
10800 FORMAT(' *    gamma + gluon_',A6,' --> q q_bar max flav =',I2,
     +'  *')
10900 FORMAT(' *      W + gluon_',A6,' --> q q_bar max flav =',I2,
     +'    *')
11000 FORMAT(' *    gamma + gluon_',A6,' --> ',A7,' + gluon',
     +'       *')
11001 FORMAT(' *        gluon + gluon   --> ',A7,' + gluon',
     +'        *')
11002 FORMAT(' *        gluon + gluon   --> ',A10,'        ',
     +'     *')
11205 FORMAT(' *               g + g     --> q + q_bar             *')
11305 FORMAT(' *               g + g     --> Q + Q_bar             *')
11405 FORMAT(' *               g + g     --> g + g                 *')
11505 FORMAT(' *               g + q     --> g + q                 *')
11515 FORMAT(' *              g* + q*    --> q                     *')
11516 FORMAT(' *           gamma* + q*   --> q                     *')
11605 FORMAT(' *               g + Q     --> g + Q                 *')
11705 FORMAT(' *               q + q_bar --> g + g                 *')
11805 FORMAT(' *               q + q_bar --> q + q_bar             *')
11905 FORMAT(' *               q + q_bar --> Q + Q_bar             *')
12005 FORMAT(' *               q + q     --> q + q                 *')
12105 FORMAT(' *               q + q     --> Q + Q                 *')
11100 FORMAT(' *      W + gluon_',A6,' --> q q_bar max flav =',I2,
     +'    *')
11200 FORMAT(' *    EPA + gamma* gluon --> q q_bar  used           *')
11300 FORMAT(' *    gluon* gluon* --> q q_bar  used                  *')
11301 FORMAT(' *     gamma gamma --> q q_bar  used                   *')
11302 FORMAT(' *       g g --> h0                                    *')
14503 FORMAT(' *       g g --> Z q qbar                            *')
14504 FORMAT(' *       g g --> Z Q Qbar                            *')
14513 FORMAT(' *       g g --> W q_i q_j                           *')
14514 FORMAT(' *       g g --> W q Q                               *')
14500 FORMAT(' *       q_i qbar_i --> Z                            *')
14510 FORMAT(' *       q_i q_j --> W                               *')
14501 FORMAT(' *       q g --> Z q                                 *')
14511 FORMAT(' *       q g --> W q                                 *')
11404 FORMAT(' *       q g --> q g                                 *')
11400 FORMAT(' *  semihard approach for BGF Catani et al           *')
11500 FORMAT(' *  semihard approach for BGF Baranov/Zotov et al    *')
11600 FORMAT(' *     gamma + gluon_',A4,
     +   ' --> Q Q_bar                *')
11700 FORMAT(' *       gluon* + gluon*',
     +   ' --> Q Q_bar                 *')
11701 FORMAT(' *       gamma  + gamma',
     +   ' --> Q Q_bar                 *')
11800 FORMAT(' *       heavy flavor produced is :',
     +    I3,'               *')
11811 FORMAT(' *       heavest flavor produced is :',
     +    I3,'             *')
11801	   format(' *       multiple scattering selected: ',I2,
     +'            *')
11900 FORMAT(' *      W + gluon_',A4,
     +   ' --> c c_bar                   *')
12000 FORMAT(' *    EPA + gamma gluon --> Q Q_bar (massive)   used   *')
12100 FORMAT(' *    scattered electron theta_max = ',F7.3,'          *')
12200 FORMAT(' *    scattered electron theta_min = ',F7.3,'          *')
12300 FORMAT(' *    Q^2 _min = ',F8.3,'                            *')
12400 FORMAT(' *    Q^2 _max = ',F8.3,'                            *')
12500 FORMAT(' *    y _min = ',F7.3,'                               *')
12600 FORMAT(' *    y _max = ',F7.3,'                               *')
12700 FORMAT('  cm energy   ',F10.3,' GeV')
17100 FORMAT(' *  cut on p_t = ',f10.2,' for IPRO = ',I5,
     +       '         *')
11888 FORMAT(' *  cut on pt applied to process 10 and 11           *')
      END
