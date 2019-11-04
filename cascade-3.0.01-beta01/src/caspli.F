      SUBROUTINE CASPLI(KF,KFLIN,KFLCH,KFLSP)
      Implicit None
C...In case of a hadron remnant which is more complicated than just a
C...quark or a diquark, split it into two (partons or hadron + parton).
      Integer KF,KFLIN,KFLCH,KFLSP,KFL,KFA,KFS,KFLR,IAGR,KFDUMP,NAGR
	Integer J,ID1,ID2,KSP
      DIMENSION KFL(3)
#include "caspara5.inc"
	Double Precision dcasrn
	Double Precision RAGR
C...Preliminaries. Parton composition.
      KFA=IABS(KF)
      KFS=ISIGN(1,KF)
      KFL(1)=MOD(KFA/1000,10)
      KFL(2)=MOD(KFA/100,10)
      KFL(3)=MOD(KFA/10,10)
      IF(KFA.EQ.22) THEN
         KFL(2)=INT(1.5+dcasrn())
         KFL(3)=KFL(2)

      ELSEIF(KFA.EQ.113.AND.dcasrn().GT.0.5) THEN
         KFL(2)=2
         KFL(3)=2
      ELSEIF(KFA.EQ.223.AND.dcasrn().GT.0.5) THEN
         KFL(2)=1
         KFL(3)=1
      ENDIF
      IF(KFLIN.NE.21.AND.KFLIN.NE.22.AND.KFLIN.NE.23) THEN
         KFLR=KFLIN*KFS
      ELSE
         KFLR=KFLIN
      ENDIF
      KFLCH=0

C...Subdivide lepton.
      IF(KFA.GE.11.AND.KFA.LE.18) THEN
c     in case of resolved photon
         IF(KFLIN.EQ.21) THEN
            KFLSP = KFLIN*KFS
            KFLCH = KFA
         ELSEIF(IABS(KFLIN).LT.6) THEN
            KFLSP = KFLIN
            KFLCH = KFA
         ELSEIF(KFLIN.EQ.22) THEN
            KFLSP = KFLIN*KFS
            KFLCH = KFA
         ENDIF
C...Subdivide photon.
      ELSEIF(KFA.EQ.22) THEN
         IF(KFLR.NE.21.AND.KFLR.LE.6) THEN
            KFLSP=-KFLR
            KFLCH = 0
         ELSE
            KFLSP = KFLIN
            KFLCH = 0
         ENDIF
c.. changed for pomeron
      ELSEIF(KFA.EQ.100) THEN
c         write(6,*) ' caspli: here for pomeron splitting KFLIN: ',KFLIN
         IF(KFLIN.EQ.21) THEN
            KFLCH = 0
            KFLSP = 21
         ELSE
            KFLCH = 0
            KFLSP = -KFLIN
         ENDIF
c          write(6,*) 'caspli: KFLCH ',KFLCH,' KFLSP = ',KFLSP
C...Subdivide meson.
      ELSEIF(KFL(1).EQ.0) THEN
         KFL(2)=KFL(2)*(-1)**KFL(2)
         KFL(3)=-KFL(3)*(-1)**IABS(KFL(2))
         IF(KFLR.EQ.KFL(2)) THEN
            KFLSP=KFL(3)
         ELSEIF(KFLR.EQ.KFL(3)) THEN
            KFLSP=KFL(2)
         ELSEIF(KFLR.EQ.21.AND.dcasrn().GT.0.5) THEN
            KFLSP=KFL(2)
            KFLCH=KFL(3)
         ELSEIF(KFLR.EQ.21) THEN
            KFLSP=KFL(3)
            KFLCH=KFL(2)
         ELSEIF(KFLR*KFL(2).GT.0) THEN
            CALL PYKFDI(-KFLR,KFL(2),KFDUMP,KFLCH)
            KFLSP=KFL(3)
         ELSE
            CALL PYKFDI(-KFLR,KFL(3),KFDUMP,KFLCH)
            KFLSP=KFL(2)
         ENDIF

C...Subdivide baryon.
      ELSE
         NAGR=0
         DO 10 J=1,3
   10    IF(KFLR.EQ.KFL(J)) NAGR=NAGR+1
         IF(NAGR.GE.1) THEN
            RAGR=0.00001+(NAGR-0.00002)*dcasrn()
            IAGR=0
            DO 20 J=1,3
               IF(KFLR.EQ.KFL(J)) RAGR=RAGR-1.
   20       IF(IAGR.EQ.0.AND.RAGR.LE.0.) IAGR=J
         ELSE
            IAGR=INT(1.00001+2.99998*dcasrn())
         ENDIF
         ID1=1
         IF(IAGR.EQ.1) ID1=2
         IF(IAGR.EQ.1.AND.KFL(3).GT.KFL(2)) ID1=3
         ID2=6-IAGR-ID1
         KSP=3
         IF(MOD(KFA,10).EQ.2.AND.KFL(1).EQ.KFL(2)) THEN
            IF(IAGR.NE.3.AND.dcasrn().GT.0.25) KSP=1
         ELSEIF(MOD(KFA,10).EQ.2.AND.KFL(2).GE.KFL(3)) THEN
            IF(IAGR.NE.1.AND.dcasrn().GT.0.25) KSP=1
         ELSEIF(MOD(KFA,10).EQ.2) THEN
            IF(IAGR.EQ.1) KSP=1
            IF(IAGR.NE.1.AND.dcasrn().GT.0.75) KSP=1
         ENDIF
         KFLSP=1000*KFL(ID1)+100*KFL(ID2)+KSP
         IF(KFLR.EQ.21) THEN
            KFLCH=KFL(IAGR)
         ELSEIF(NAGR.EQ.0.AND.KFLR.GT.0) THEN
            CALL PYKFDI(-KFLR,KFL(IAGR),KFDUMP,KFLCH)
         ELSEIF(NAGR.EQ.0) THEN
            CALL PYKFDI(10000+KFLSP,-KFLR,KFDUMP,KFLCH)
            KFLSP=KFL(IAGR)
         ENDIF
c	   write(6,*) ' caspli: ID1,ID2,KSP,IAGR: ',ID1,ID2,KSP,IAGR
c	   write(6,*) ' caspli: KFLSP,KFLCH,KFS: ',KFLSP,KFLCH,KFS
      ENDIF

C...Add on correct sign for result.
      KFLCH=KFLCH*KFS
      KFLSP=KFLSP*KFS

      RETURN
      END
