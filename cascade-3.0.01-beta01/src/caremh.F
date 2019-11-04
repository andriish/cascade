      SUBROUTINE CAREMH(IFLR,K2,CHI)
	Implicit None
#include "casprkt.inc"

C...Gives flavour and energy-momentum fraction Z for the particle
C...to be produced out of the target remnant when that is not a
C...simple diquark.
      Double Precision dcasrn
      Double Precision X,Y,Q2LP
      DOUBLE PRECISION PARL
	Integer LST
      COMMON /CASTOU/PARL(30),X,Y,Q2LP,LST(40)
	Integer IFLR,KSP,K2,IDUM,K2A,PYCOMP,KC2
	Double Precision AMSP,AMK2,Z,FC,A
	Double Precision CHI,ZZ,TM2,PT
	Double precision PYMASS
      EXTERNAL dcasrn

      KSP=IFLR
C...Split energy-momentum of target remnant according to functions P(z)
C...z=E-pz fraction for qq (q) forming jet-system with struck Q (Qbar)
C...1-z=E-pz fraction for qQbar (qqQ) hadron
C...mq=mass of (light) parton remnant q (qq) in jet system
C...mQ=mass of produced (heavy flavour) hadron
      AMSP=PYMASS(KSP)
      AMK2=PYMASS(K2)
      KC2=IABS(PYCOMP(K2))
	
      IF(IRspl.EQ.1) THEN
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=1/3 since a=1 fixed
      Z=1.-SQRT(dcasrn())
      IF(KC2.GE.301.AND.KC2.LE.400) Z=1.-Z
      ELSEIF(IRspl.EQ.2) THEN
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=mq/(mq+mQ) --> a=a(mq,mQ)
         A=(AMSP+AMK2)/AMSP - 2.
         Z=dcasrn()**(1./(A+1.))
      ELSEIF(IRspl.EQ.3) THEN
C...Using Peterson fragmentation function
C...P(z)=N/(z(1-1/z-c/(1-z))**2)  where c=(mq/mQ)**2  (FC=-c)
         FC=-(AMSP/AMK2)**2
   30    Z=dcasrn()
         IF(-4.*FC*Z*(1.-Z)**2.LT.dcasrn()
     +	   *((1.-Z)**2-FC*Z)**2) GOTO 30
	Else
	  write(6,*) ' caremh: parameter IRspl not valid: ',IRspl
	  stop
      ENDIF
      LST(27)=1
      K2A=IABS(K2)
      IF((K2A.GE.1.AND.K2A.LE.8).OR.K2A.EQ.21.OR.PYCOMP(K2A).EQ.90)
     +LST(27)=2
      CHI = Z
c      write(6,*) ' caremh : ',Z,IFLR,K2,AMSP,AMK2
      RETURN
      END
