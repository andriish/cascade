      SUBROUTINE CAREMN(IPU1,IPU2)
      IMPLICIT None
C...ADDS ON TARGET REMNANTS (ONE OR TWO FROM EACH SIDE) AND
C...INCLUDES PRIMORDIAL KT.
      DOUBLE PRECISION DBETAX,DBETAZ,DROBO(5)
      Double Precision X,Y,Q2LP
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /CASTOU/PARL(30),X,Y,Q2LP,LLST(40)
#include "caslujets.inc"
#include "casludat1.inc"
#include "casgki.inc"
#include "casluco.inc"
#include "caspara1.inc"
#include "casparas.inc"
#include "casprkt.inc"
#include "caspara5.inc"
#include "casherup.inc"
#include "caspsipol.inc"
#include "caspar.inc"

      DOUBLE PRECISION ROBO,SINTH,PHIPT
      Double Precision XPY,SH
      Integer KFL
      COMMON /MYPROC/KFL(3,2),XPY(2),SH
      Double Precision PMS,PSYS
      Integer KFLCH,KFLSP,IS,ISN
      DIMENSION KFLCH(2),KFLSP(2),PMS(0:6),IS(2),ROBO(5)

c.hju
      DIMENSION PSYS(0:2,5),ISN(2)

      Double Precision GGPT,GPPT
      COMMON /CASKT/ GGPT,GPPT
      Double Precision prkt

      Double Precision PTSPL,PHISPL,CHI(2),zz,pmst
      Integer I,J,JT,IFLS,IPU,IQ,IPU1,IPU2,ILEP,IP,NS,NTRY
      Integer IL,IR,IG,IMAX,IMIN
      Double Precision PEI,PE,PZI,PZ,SHS,SHR,PMMIN0,PMMIN1
      Double Precision PMMIN2,PMMIN,PPB,PNB,PMTB,PMTR,PMTL,PW1
      Double Precision SQLAM,SQSGN,RKR,RKL,BER,BEL,PZF,PT2,RQP
      Double Precision PEM,PZM,PDEV
      double precision chmin,chmax,wtz,wtmax
	Double Precision px,py
	Common/carespt/px(2),py(2)
      
      Double Precision cmspart(4),dot
      
      Double Precision dcasrn
	Integer kfmy, icallz
      Integer Iremn
      common/casflags/Iremn
      Integer K3remn
	Double Precision done
	data done/.999999999d0/
      External dcasrn
      DATA IPU,IQ/0,0/,PEI,PE,PZI,PZ,SHS/5*0./
      
c      Iremn = 1 
c      Iremn = 0  ! switch off remnant for hadronisaton etc.    
c      write(6,*) ' caremn ',Iremn  
      K3remn = 3
      If(Iremn.eq.0) K3remn = 21 
      
      GGPT = -9999.d0
      GPPT = -9999.d0
      IF(IPU1.EQ.0.AND.IPU2.EQ.0) then
        Write(6,*) ' CAREMN: IPU1/2 RETURN'
        RETURN
      Endif
c      write(6,*) ' caremn ',IPU1,IPU2,Isemih
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      IF(ILEP.EQ.1) IQ=21
      IF(ILEP.EQ.2) IQ=23
      IF(ILEP.eq.1.and.ILEPTO.eq.0) ILEP=0

      IF(ILEP.EQ.1.OR.(IRES(1).EQ.0.and.ILEPTO.EQ.1)) Q2 = Q2LP
      PYVAR(33) =  0
      PYVAR(34) =  0
c      write(6,*) ' CAREMN ',ILEP,PYVAR(31),PYVAR(32)
      IP=MAX(IPU1,IPU2)
      NS=N
C...DEFINE INITIAL PARTONS, INCLUDING PRIMORDIAL KT
   10 DO 30  I=3,4
         IF(I.EQ.3) IPU=IPU1
         IF(I.EQ.4) IPU=IPU2
         K(I,1)=21
         K(I,3)=I-2
         DO 20 J=1,5
   20    P(I,J)=0.D0
         IF(IPU.NE.0) THEN
            K(I,2)=K(IPU,2)
            P(I,5)=P(IPU,5)
            PTSPL=0
            P(I,5)=0.D0
            P(I,1)=PTSPL*COS(PHISPL)
            P(I,2)=PTSPL*SIN(PHISPL)
            PMS(I-2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
         ELSE
            K(I,2)=K(IQ,2)
            IF(ILEPTO.EQ.0) THEN
               P(I,5)= 0.0D0
               pms(i-2)=0.0D0

            ELSEIF(ILEPTO.EQ.1) THEN
               P(I,5)=-SQRT(Q2)
               PMS(I-2)= - Q2
c. shs = w**2 but calculated for photoprod.
               SHS= (- Q2 + PARL(22)+PYVAR(7-I)**2)
            ENDIF
         ENDIF
   30 CONTINUE
C...KINEMATICS CONSTRUCTION FOR INITIAL PARTONS
      IF(ILEP.EQ.0) SHS=PYVAR(31)*PYVAR(32)*PYVAR(2)
     + +(P(3,1)+P(4,1))**2+(P(3,2)+P(4,2))**2
      SHR=DSQRT(DMAX1(0.D0,SHS))
      IF(ILEP.EQ.0) THEN
         IF((SHS-PMS(1)-PMS(2))**2-4.D0*PMS(1)*PMS(2).LE.0.) then
	   IPY(48) = 1
	   RETURN
	   ENDIF
         P(3,4)=0.5D0*(SHR+(PMS(1)-PMS(2))/SHR)
         P(3,3)=DSQRT(DMAX1(0.D0,P(3,4)**2-PMS(1)))
         P(4,4)=SHR-P(3,4)
         P(4,3)=-P(3,3)
	   
      ELSEIF(ILEP.EQ.1) THEN
         P(3,4)=P(IQ,4)
         P(3,3)=P(IQ,3)
         P(4,4)=P(IP,4)
         P(4,3)=P(IP,3)
      ELSEIF(ILEP.EQ.2) THEN
         P(3,4)=P(IP,4)
         P(3,3)=P(IP,3)
         P(4,4)=P(IQ,4)
         P(4,3)=P(IQ,3)
      ENDIF
C...TRANSFORM PARTONS TO OVERALL CM-FRAME (NOT FOR LEPTOPRODUCTION)
      IF(ILEP.EQ.0) THEN
         MSTU(1)=3
         MSTU(2)=4
         DROBO(3)=(P(3,1)+P(4,1))/SHR
         DROBO(4)=(P(3,2)+P(4,2))/SHR
         CALL PYROBO(MSTU(1),MSTU(2),0.D0,0.D0,-DROBO(3),-DROBO(4),0.D0)
         ROBO(2)=PYANGL(P(3,1),P(3,2))
         CALL PYROBO(MSTU(1),MSTU(2),0.D0,-ROBO(2),0.D0,0.D0,0.D0)
         ROBO(1)=PYANGL(P(3,3),P(3,1))
         CALL PYROBO(MSTU(1),MSTU(2),-ROBO(1),0.D0,0.D0,0.D0,0.D0)
         MSTU(2)=MAX(IPY(47),IPU1,IPU2)
         CALL PYROBO(MSTU(1),0,ROBO(1),ROBO(2),DROBO(3), DROBO(4),0D0)
         DROBO(5)=MAX(-done,MIN(done,(PYVAR(31)-PYVAR(32))/
     +   (PYVAR(31)+PYVAR(32))))
         CALL PYROBO(MSTU(1),0,0.D0,0.D0,0.D0,0.D0,DROBO(5))
         MSTU(1)=0
         MSTU(2)=0

      ENDIF
C...CHECK INVARIANT MASS OF REMNANT SYSTEM:
C...HADRONIC EVENTS OR LEPTOPRODUCTION
      IF(ILEP.LE.0) THEN
c here is changed; now according to PYTHIA 
         PSYS(0,4)=P(3,4)+P(4,4)
         PSYS(0,3)=P(3,3)+P(4,3)
         PMS(0)=DMAX1(0.D0,PSYS(0,4)**2 - PSYS(0,3)**2)
         PMMIN0=DSQRT(PMS(0))
         PMMIN1=P(1,5)**2 + P(3,1)**2 + P(3,2)**2
         PMMIN1=DSQRT(PMMIN1)
         PMMIN2=P(2,5)**2 + P(4,1)**2 + P(4,2)**2
         PMMIN2=DSQRT(ABS(PMMIN2))
         PSYS(1,4) = 0.5*PYVAR(1)*(1. - PYVAR(30+1))
         PSYS(2,4) = 0.5*PYVAR(1)*(1. - PYVAR(30+2))
         IF(IRES(1).EQ.0.AND.(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13)
     +   ) PSYS(1,3)=PSYS(1,4)
         IF(PMMIN0+PMMIN1+PMMIN2.GT.PYVAR(1).OR.PMMIN1.GT.PSYS(1,4)
     +   .OR. PMMIN2.GT.PSYS(2,4)) THEN
c            write(6,*) ' caremn PYVAR(31)=',PYVAR(31),' PYVAR(32)=',PYVAR(32)
c            write(6,*) ' caremn PMMIN0,PMMIN1,PMMIN2 ',PMMIN0,PMMIN1,
c     +      PMMIN2
c            write(6,*) ' caremn PYVAR(1) ',PYVAR(1)
c            write(6,*) ' caremn PMMIN1=',PMMIN1,' GT PSYS(1,4)=',
c     +      PSYS(1,4)
c            write(6,*) ' caremn PMMIN2=',PMMIN2,' GT PSYS(2,4)=',
c     +      PSYS(2,4)
            IPY(48)=1
            RETURN
         ENDIF
      ELSE
         PEI=P(IQ,4)+P(IP,4)
         PZI=P(IQ,3)+P(IP,3)
         PMS(ILEP)=DMAX1(0.D0, PEI**2-PZI**2+P(5-ILEP,1)**2+
     +   P(5-ILEP,2)**2)
         PMMIN=P(3-ILEP,5)+PYMASS(K(5-ILEP,2))+DSQRT(PMS(ILEP))
         PMMIN=PMMIN+PYPAR(12)

         IF(SHR.LE.PMMIN) THEN
c            write(6,*) ' caremn : shr < PMmin ',shr,pmmin
            IPY(48)=1
            RETURN
         ENDIF
      ENDIF
C...SUBDIVIDE REMNANT IF NECESSARY, STORE FIRST PARTON
      NTRY = 0
   40 I=NS-1
      DO 70  JT=1,2
         ISN(JT)=0
         IF((ILEP.EQ.JT).OR. (IPY(40+JT)
     +   .EQ.22.AND.IRES(JT).EQ.0)) GOTO 70
         IF(JT.EQ.1) IPU=IPU1
         IF(JT.EQ.2) IPU=IPU2
c         CALL caspli(IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT))
         CALL pyspli(IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT))
c         write(6,*) ' caremn ',JT,IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT)
         I=I+2
         IS(JT)=I
         ISN(JT)=2
         K(I,1)=K3remn
         K(I,2)=KFLSP(JT)
         K(I,3)=JT
         P(I,5)=PYMASS(K(I,2))
c cheat with masses for ILHA 
         if(ilha.gt.0) P(I,5)= PYMASS( KFLCH(JT))

C...FIRST PARTON COLOUR CONNECTIONS AND TRANSVERSE MASS
         K(I+1,1)=-1
         K(I+1,3)=I
         K(I+1,2)=1000
         IF(IPY(34).GE.1) K(I+1,2)=1000+JT
         DO 50 J=1,5
   50    P(I+1,J)=0.D0
         IF(IPU.NE.0.AND.KFLSP(JT).EQ.21) THEN
            P(I+1,3)=IPU
            P(I+1,4)=IPU
            P(IPU+1,1)=I
            P(IPU+1,2)=I
            K(I,4)=IPU+IPU*MSTU(5)
            K(I,5)=IPU+IPU*MSTU(5)
            K(IPU,4)=MOD(K(IPU,4),MSTU(5))+I*MSTU(5)
            K(IPU,5)=MOD(K(IPU,5),MSTU(5))+I*MSTU(5)
         ELSEIF(IPU.NE.0) THEN
            IFLS=(3-ISIGN(1,KFLSP(JT)*(1102-IABS(KFLSP(JT)))))/2
            P(I+1,IFLS+2)=IPU
            P(IPU+1,3-IFLS)=I
            K(I,IFLS+3)=IPU
            K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
         ENDIF
         IF(KFLCH(JT).EQ.0) THEN
            P(I,1)=-P(JT+2,1)
            P(I,2)=-P(JT+2,2)
            IF(ISEMIH.EQ.1) THEN
               P(I,1) = - P(IPU,1)
               P(I,2) = - P(IPU,2)
c               write(6,*) ' caremn',P(I,1),P(I,2),IPU
            ENDIF
            PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
         ELSE

C...WHEN EXTRA REMNANT PARTON OR HADRON: FIND RELATIVE PT, STORE
C...PRIMORDIAL KT SPLIT SHARED BETWEEN REMNANTS
            CALL CAPRIKT(PRKTREM,PTSPL,PHISPL)
            IF((IABS(K(JT,2)).EQ.11.OR.IABS(K(JT,2)).EQ.13).AND.
     +      IRES(JT).NE.0) THEN
               CHI(JT)=(1.d0 - yy)/(1.d0 - PYVAR(30+JT))
            ELSEIF((IABS(K(JT,2)).EQ.11.OR.IABS(K(JT,2)).EQ.13).AND.
     +      IRES(JT).EQ.0)  THEN
               CHI(JT)=YY
            ELSE
C...RELATIVE DISTRIBUTION OF ENERGY; EXTRA PARTON COLOUR CONNECTION
               If(IRspl.ge.1.and.IRspl.le.3) then
123               CALL caremh(KFLSP(JT),KFLCH(JT),CHI(JT))
c Oct 2016: check 0.65 limit for other processes
                  if((chi(jt).gt.0.65).and.ipro.ne.-1) goto 123
		   else
C...Using chosen fragmentation function in JETSET/PYTHIA
                  PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
			pmst=pms(jt+4)
                  icallz = 0
124               If(IRspl.eq.4) then 
                    icallz = icallz + 1
			  CALL PYZDIS(1,0,pmst,ZZ)
			elseif(IRspl.eq.5) then   
                    CALL PYZDIS(1000,0,pmst,ZZ)
			else
			  write(6,*) ' CAREMN IRspl not implemented ', IRspl
			  stop
			endif
                  if(icallz.gt.300) then
                     write(6,*) ' caremn icallz limit '
                     IPY(48)=1
                     RETURN
                  endif
c Oct 2016: check 0.65 limit for other processes
                  if((zz.gt.0.65).and.ipro.ne.-1) goto 124
 			CHI(JT)=ZZ
		   endif
            ENDIF
            P(I,1)=-P(JT+2,1)*(1.D0-CHI(JT))
     +       +PTSPL*COS(PHISPL)
            P(I,2)=-P(JT+2,2)*(1.D0-CHI(JT))
     +       +PTSPL*SIN(PHISPL)
            IF(ISEMIH.EQ.1) THEN
               P(I,1)=-P(IPU,1)*(1.d0-CHI(JT))+PTSPL*COS(PHISPL)
               P(I,2)=-P(IPU,2)*(1.d0-CHI(JT))+PTSPL*SIN(PHISPL)
c               write(6,*) ' carmen1 ',P(I,1),P(I,2),IPU,ilha
c               if(ilha.gt.0) then 
               if(ilha.ne.0) then 
c                 P(I,1)=0.
c                 P(I,2)=0.
               endif
            ENDIF

            PMS(JT+2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
            I=I+2
            ISN(JT)=4
            DO 60 J=1,5
               K(I,J)=0
               K(I+1,J)=0
               P(I,J)=0.D0
   60       P(I+1,J)=0.D0
            K(I,1)=1
            K(I,2)=KFLCH(JT)
c            write(6,*) ' caremn ',KFLCH(JT),IPY(34)
            K(I,3)=JT
            P(I,5)=PYMASS(K(I,2))
            P(I,1)=-P(JT+2,1)*CHI(JT)-PTSPL*COS(PHISPL)
            P(I,2)=-P(JT+2,2)*CHI(JT)-PTSPL*SIN(PHISPL)
            IF(ISEMIH.EQ.1) THEN
               P(I,1)=-P(IPU,1)*CHI(JT)-PTSPL*COS(PHISPL)
               P(I,2)=-P(IPU,2)*CHI(JT)-PTSPL*SIN(PHISPL)
c               write(6,*) ' caremn2 ',P(I,1),P(I,2),P(IPU,1),ilha,IPU
c                if(ilha.gt.0) then 
                if(ilha.ne.0) then 
c                  P(I,1)=-P(IPU,1)
c                  P(I,2)=-P(IPU,2)
               endif
           ENDIF
            PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
C...end of update
            PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/
     +               (1.D0-CHI(JT))
            K(I+1,1)=-1
            K(I+1,3)=I
            K(I+1,2)=1000
            IF(IPY(34).GE.1) K(I+1,2)=1000+JT
            IF((IABS(KFLCH(JT)).GE.1.AND.IABS(KFLCH(JT)).LE.8).OR.
     +      IABS(KFLCH(JT)).EQ.21.OR.pycomp(IABS(KFLCH(JT))).EQ.90)
     +      THEN
               IFLS=(3-ISIGN(1,KFLCH(JT)*(1102-IABS(KFLCH(JT)))))/2
               P(I+1,IFLS+2)=IPU
               P(IPU+1,3-IFLS)=I
               K(I,1)=K3remn
               K(I,IFLS+3)=IPU
               K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
            ELSE
               IF(IPY(34).GE.1) THEN
                  K(I,1)=1
                  K(I,3)=JT
               ENDIF
               K(I,1)=K3remn
            ENDIF
         ENDIF
   70 CONTINUE
      IF(ILEP.GE.1.AND.
     +   SHR.LE.(DSQRT(PMS(1))+DSQRT(DABS(PMS(2))))) THEN
         NTRY = NTRY + 1
         IF(NTRY.LT.500) GOTO 40
         IPY(48) = 2
c         write(6,*) ' casremn IPY(48)=2 ',IPY(48)
         RETURN
      ENDIF


      IF(ILEP.EQ.0.AND.(PMS(1).GT.PSYS(1,4)**2.OR.
     +   PMS(2).GT.PSYS(2,4)**2)) THEN
         NTRY = NTRY + 1
         IF(NTRY.LT.500) GOTO 40
         IPY(48) = 2
         RETURN
      ENDIF
      N=I+1
C...RECONSTRUCT KINEMATICS OF REMNANTS
      DO 80  JT=1,2
         IF((ILEPTO.EQ.1.AND.IRES(JT).EQ.0).OR.(IPY(40+JT)
     +   .EQ.22.AND.IRES(JT).EQ.0)) GOTO 80
         IF(ILEP.EQ.1) THEN
            PE=0.5D0*(SHR+(PMS(JT)-PMS(3-JT))/SHR)
            PZ=DSQRT(PE**2-PMS(JT))
         ELSE
c new ala PYTHIA  
            PE=PSYS(JT,4)
            PZ=DSQRT(PE**2-PMS(JT))
         ENDIF
         IF(KFLCH(JT).EQ.0) THEN
            P(IS(JT),4)=PE
            P(IS(JT),3)=PZ*(-1)**(JT-1)
         ELSE
            PW1=CHI(JT)*(PE+PZ)
            P(IS(JT)+2,4)=0.5D0*(PW1+PMS(JT+4)/PW1)
            P(IS(JT)+2,3)=0.5D0*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)
            P(IS(JT),4)=PE-P(IS(JT)+2,4)
            P(IS(JT),3)=PZ*(-1)**(JT-1)-P(IS(JT)+2,3)
         ENDIF
         PSYS(JT,3)=PZ*(-1)**(JT-1)
   80 CONTINUE


C...HADRONIC EVENTS: BOOST REMNANTS TO CORRECT LONGITUDINAL FRAME
      IF(ILEP.LE.0) THEN
ccc      CALL PYLIST(1)
C new a la PYTHIA
C... check if longitudinal boost needed - if so pick two systems.
         PDEV = ABS(PSYS(0,4)+PSYS(1,4)+PSYS(2,4) -
     +    PYVAR(1) + ABS(PSYS(0,3)+PSYS(1,3)+PSYS(2,3)))
c         write(6,*) ' CAREMN: PDEV return ',PDEV,PYVAR(1),ISN(1),ISN(2)
c          call pylist(1)
         IF(PDEV.LE.1.E-7*PYVAR(1).OR. (K(1,2).EQ.22.AND.IRES(1).EQ.0))
     +      then
c            write(6,*) ' CAREMN: PDEV return ',PDEV,PYVAR(1),ISN(1),ISN(2)
c            CALL PYLIST(1)
            RETURN
         Endif
         IF(ISN(1).EQ.0) THEN
            IR=0
            IL=2
         ELSEIF(ISN(2).EQ.0) THEN
            IR=1
            IL=0
         ELSEIF((1.-PYVAR(31)).GT.0.05.AND.(1.-PYVAR(32)).GT.0.05) THEN
            IR=1
            IL=2
         ELSEIF((1.-PYVAR(31)).GT.0.05) THEN
            IR=1
            IL=0
          ELSEIF((1.-PYVAR(32)).GT.0.05) THEN
           IR=0
            IL=2
         ELSEIF(PMS(1)/PSYS(1,4)**2.GT.PMS(2)/PSYS(2,4)**2) THEN
            IR=1
            IL=0
         ELSE
            IR=0
            IL=2
         ENDIF
         IG=3-IR-IL

C...Construct longitudinal boost
         IF((IG.EQ.1.AND.ISN(1).EQ.0).OR.(IG.EQ.2.AND.ISN(2).EQ.0))
     +   THEN
            PPB=PYVAR(1)
            PNB=PYVAR(1)
         ELSE
            PPB=PYVAR(1)-PSYS(IG,4)-PSYS(IG,3)
            PNB=PYVAR(1)-PSYS(IG,4)+PSYS(IG,3)
         ENDIF

         PMTB=PPB*PNB
         PMTR=PMS(IR)
         PMTL=PMS(IL)
         SQLAM=DSQRT(DMAX1(0.D0,(PMTB-PMTR-PMTL)**2-4.D0*PMTR*PMTL))
         SQSGN=DSIGN(1.D0,PSYS(IR,3)*PSYS(IL,4)-PSYS(IL,3)*PSYS(IR,4))
         RKR=(PMTB+PMTR-PMTL+SQLAM*SQSGN)/(2.D0*(PSYS(IR,4)+PSYS(IR,3))*
     +   PNB)
         RKL=(PMTB+PMTL-PMTR+SQLAM*SQSGN)/(2.D0*(PSYS(IL,4)-PSYS(IL,3))*
     +   PPB)
         BER=(RKR**2-1.D0)/(RKR**2+1.D0)
         BEL=-(RKL**2-1.D0)/(RKL**2+1.D0)

C...Perform longitudinal boost
         MSTU(1)=0
         MSTU(2)=0
         IF(IR.EQ.1) THEN
           CALL PYROBO(IS(1),IS(1)+ISN(1)-2,0.D0,0.D0,0D0,0D0,BER)
         ELSE
           CALL PYROBO(3,NS,0.D0,0.D0,0D0,0D0,BER)
         ENDIF
         IF(IL.EQ.2) THEN
           CALL PYROBO(IS(2),IS(2)+ISN(2)-2,0.D0,0.D0,0D0,0D0,BEL)
         ELSE
           CALL PYROBO(3,NS,0.D0,0.D0,0D0,0D0,BEL)
         ENDIF
c        write(6,*) ' caremn after boost',IS(1),ISN(1),IS(2),ISN(2),IR,IL
c        call PYLIST(1)
C...LEPTOPRODUCTION EVENTS: BOOST COLLIDING SUBSYSTEM
      ELSE

         IMIN=21
         if(ipsiel1.ge.1) imin=14
         MSTU(1) = IMIN
         IMAX=MAX(IP,IPY(47))
         MSTU(2) = IMAX
         PZF=PZ*(-1)**(ILEP-1)
         PT2=P(5-ILEP,1)**2+P(5-ILEP,2)**2
         PHIPT=PYANGL(P(5-ILEP,1),P(5-ILEP,2))
         CALL PYROBO(IMIN,IMAX,0.D0,-PHIPT,0.D0,0.D0,0.D0)
         RQP=P(IQ,3)*(PT2+PEI**2)-P(IQ,4)*PEI*PZI
         SINTH=P(IQ,4)*DSQRT(PT2*(PT2+PEI**2)/(RQP**2+PT2* P(IQ,4)**2*
     +   PZI**2))*DSIGN(1.D0,-RQP)
         CALL PYROBO(IMIN,IMAX,DASIN(SINTH),0.D0,0.D0,0.D0,0.D0)
         DBETAX=(-PEI*PZI*SINTH+ DSQRT(PT2*(PT2+PEI**2-
     +      (PZI*SINTH)**2)))/ (PT2+PEI**2)
         CALL PYROBO(IMIN,IMAX,0.D0,0.D0,DBETAX,0.D0,0.D0)
         CALL PYROBO(IMIN,IMAX,0.D0,PHIPT,0.D0,0.D0,0.D0)
         PEM=P(IQ,4)+P(IP,4)
         PZM=P(IQ,3)+P(IP,3)
         DBETAZ=(-PEM*PZM+
     +   PZF*DSQRT(PZF**2+PEM**2-PZM**2))/(PZF**2+PEM**2)
         CALL PYROBO(IMIN,IMAX,0.D0,0.D0,0.D0,0.D0,DBETAZ)
         MSTU(1) = 0
         MSTU(2) = 0
         CALL PYROBO(3,4,DASIN(SINTH),0.D0,DBETAX,0.D0,0.D0)
         CALL PYROBO(3,4,0.D0,PHIPT,0.D0,0.D0,DBETAZ)
         PEM=P(IQ,4)+P(IP,4)
         PZM=P(IQ,3)+P(IP,3)
      ENDIF
c      write(6,*) ' end caremn '
c      call PYLIST(1)
      RETURN
      END
